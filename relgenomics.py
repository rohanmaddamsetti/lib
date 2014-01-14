#!/usr/bin/python

'''relgenomics.py by Rohan Maddamsetti and David Groendyk
This module contains code that compares evolutionary rates between the LTEE
and wild E. coli strains. More functionality may be added as needed in the
future.
'''

from os import listdir, rename, getcwd, remove
from os.path import basename, exists
from sys import stdout, exit
import subprocess
from Bio import AlignIO
from Bio.Seq import Seq
from pprint import pprint
from codeml import Codeml, read
import csv
import pickle

def backtranslate(prot_aln_file, dna_lookup_file, outfile):
    '''Input: a protein alignment file produced by SATe, and a DNA
    sequence file in FASTA format. These DNA sequences MUST correspond
    to the sequences in the protein alignment. If not, a fatal error
    is raised. A filename for the output file is required.
    Output: a DNA alignment in FASTA format in the file named by outfile.
    '''
    prot_aln = open(prot_aln_file, 'r')
    dna_lookup = open(dna_lookup_file, 'r')
    outhandle = open(outfile, 'w')
    prot_dict = {}
    dna_dict = {}

    for line in prot_aln:
        line = line.strip()
        if len(line):
            if line.startswith('>'):
                seqid = line[1:]
            else:
                prot_dict[seqid] = line
 
    for line in dna_lookup:
        line = line.strip()
        if len(line):
            if line.startswith('>'):
                seqid = line[1:]
                dna_dict[seqid] = ""
            else:
                dna_dict[seqid] = dna_dict[seqid] + line

    for seqid, protseq in prot_dict.iteritems():
        matching_dnaseq = dna_dict[seqid]
        ##Verify that the dnaseq matches the protseq.
        ##NOTA BENE: This might fail for non-standard translations.
        tr_dnaseq = str(Seq(matching_dnaseq).translate(table="Bacterial"))
        nogap_protseq = protseq.replace('-','')
        ##non-canonical amino acids X (unknown), U (selenocysteine),
        ## and O (pyrrolysine) are converted to match internal stops.
        ##this is to prevent fatal errors being raised when valid
        ##internal 'stop' codons are used to recode for non-canonical amino
        ## acids.
        nogap_protseq = nogap_protseq.replace('X', '*')
        nogap_protseq = nogap_protseq.replace('U', '*')
        nogap_protseq = nogap_protseq.replace('O', '*')
        if tr_dnaseq != nogap_protseq:
                        raise TypeError("FATAL ERROR: The protein sequence " + 
                                        nogap_protseq + " doesn't match " +
                                        "the translated dna sequence: " + tr_dnaseq + "\n")
            
        outhandle.write('>' + seqid + '\n')
        ## Print each DNA triplet per amino acid.
        ## If there's a gap char, print three gap chars.
        pos = 0
        for residue in protseq:
            if residue == '-':
                outhandle.write('---')
            else:
                codon = matching_dnaseq[3*pos:3*pos+3]
                outhandle.write(codon)
                pos = pos + 1
        outhandle.write('\n')

    return None
            
def align_using_sate(file_of_fasta_sequences, datatype="Protein"):
    '''Input: a FASTA file with the FASTA sequences of gene sequences.
    Output: one file of the aligned genes, one file that contains the gene tree,
    one file of the output, one file of errors encountered, and one file
    that contains the score number of the gene tree.
    '''

    filename, sep, ext = basename(file_of_fasta_sequences).partition('.')
    process_args = ['/Users/Rohan/Desktop/Projects/lib/external/sate/sate-core/run_sate.py', '--datatype=' + datatype, '--iter-limit=1',
                    '--num-cpus=3', '--aligner=prank', '--job='+filename,
                    '--input=' + file_of_fasta_sequences, '--output-directory=' + getcwd()]
    process = subprocess.Popen(process_args, shell=False, stdout=subprocess.PIPE)
    print "the process arguments were: ", process_args
    my_stdout, my_stderr = process.communicate()
    if my_stderr is not None:
        print my_stdout
    return None

def convert_fasta_to_phylip(fasta_alignment, phylip_file_name):
    '''Input: a FASTA file of aligned protein sequences.
    Output: a PHYLIP file of the aligned protein sequences.
    If the length of the alignment is not divisible by three, one or two
    gap characters are added to the end of the file .
    '''
    alignment = AlignIO.read(fasta_alignment, "fasta")
    mod_alignment = len(alignment[0].seq) % 3
    while mod_alignment != 0:
        for i, record in enumerate(alignment):
            alignment[i].seq = record.seq + '-'
        mod_alignment = len(alignment[0].seq) % 3
    handle = open(phylip_file_name, "w")
    AlignIO.write(alignment, handle, "phylip")

def correct_phylip_format(old_phylip, new_phylip):
    '''Input: a phylip format alignment file of rpoB genes.
    Output: a phylip format alignment file of rpoB genes with the correct interleaved
    file format.
    '''
    first_line = 1
    for line in old_phylip:
        if first_line:
            new_phylip.write(line.rstrip("\n") + ' I' + "\n")
            first_line = 0
        else:
            new_phylip.write(line)    

def calculate_average_dN(codeml_results):
    '''Input: a dictionary of results from codeml.
    Output: the average dN for all branches.
    '''
    list_of_dN = []
    for gene, match in codeml_results['pairwise'].iteritems():
        for pair, value in match.iteritems():
            list_of_dN.append(value['dN'])
    try:
        mean_dN = sum(list_of_dN) / len(list_of_dN)
    except ZeroDivisionError:
        mean_dN = 10000
    return mean_dN

def calculate_average_omega(codeml_results):
    '''Input: a dictionary of results from codeml.
    Output: the average dN/dS ratio for all branches.
    '''
    list_of_dN = []
    list_of_dS = []
    for gene, match in codeml_results['pairwise'].iteritems():
        for pair, value in match.iteritems():
            list_of_dN.append(value['dN'])
            list_of_dS.append(value['dS'])
    try:
        omega = sum(list_of_dN) / sum(list_of_dS)
    except ZeroDivisionError:
        omega = 10000
    return omega

def read_panortholog_families (panorthologs_path, dir,dna=True):
    '''Input: the path to the hatcher panortholog file, and a directory
    containing the proteomes indexed in the hatcher panortholog file.
    The dna parameter is True if the fasta sequences are DNA sequences,
    otherwise, it will look up protein sequences.
    Output: A dictionary of locus_ids to a string containing fasta sequences, called families.
    '''

    if dna:
        suffix = ".nuc"
    else:
        suffix = ".proteins"
    ##First, map all sequence ids to locus tags in the REL606 proteome.
    seq_ids_to_tags = {}
    rel606_proteome = open("/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/hatcher_results/proteins/E-coli-B-REL606.proteins")
    for line in rel606_proteome:
        if '>' in line:
            words = line.split()
            seqid = words[0][1:]
            locustag = words[1]
            seq_ids_to_tags[seqid] = locustag
    rel606_proteome.close()
    
    families = {}
    hatcher_panorthologs = open(panorthologs_path)
    for entry in hatcher_panorthologs:
        entry = entry.strip()
        splitline = entry.split("\t")
        current_fasta_seqs = ""
        rel606seq_id = ""
        for elt in splitline:
            genome_name, sep, sequence_id = elt.partition('$')
            if sep:
                seqfile = open(dir + genome_name + suffix)
                found = 0
                if 'REL606' in genome_name:
                    rel606seq_id = sequence_id
                for line in seqfile:
                    if not found:
                        if sequence_id in line:
                            current_fasta_seqs = current_fasta_seqs + '>' + sequence_id + "\n";
                            found = 1
                    else:
                        if '>' in line:
                            break
                        else:
                            current_fasta_seqs = current_fasta_seqs + line
        current_locustag = seq_ids_to_tags[rel606seq_id]
        families[current_locustag] = current_fasta_seqs
    return(families)
                
def ltee_locustag_subset(families, ltee_hits):
    '''Input: a dictionary of locus tags to strings, each containing a set of
    FASTA sequences corresponding to a panortholog family. Also, a file containing
    locus tags of relevant genes mutated in the LTEE.
    Output: the dictionary, subsetted on those locus tags relevant to the LTEE.
    '''
    ltee_loci = open(ltee_hits)
    tags = []
    for line in ltee_loci:
        locustag = line.split().pop(0)
        tags.append(locustag)
    #remove the header element from tags.
    tags.pop(0)
    filtered_fams = {}
    for k, v in families.iteritems():
        if k in tags:
            filtered_fams[k] = v
    return filtered_fams
    
def use_codeml(locustag):

    dN = None
    dNdS = None
    print "current working directory is", getcwd()
    alignmentfile = locustag + '.phy'
    treefile = locustag + '.tre'
    print "alignment file and tree file are respectively: ", alignmentfile, treefile
    ##Catch errors caused by internal stop codons in codeml input.
    ##These should be due SOLELY due either nonsense mutations,
    ##or to non-canonical amino acids in the protein sequence.
    try:
        test_codeml = Codeml(alignment=alignmentfile, tree=treefile,
                             working_dir=getcwd(), out_file='codeml_results.txt')
        test_codeml.set_option('seqtype', 1)
        test_codeml.set_option('model', 1)
        test_codeml.set_option('fix_kappa', 1)
        test_codeml.set_option('runmode', -2)
        # In test_codeml.run(), if a control function is specified as an argument,
        # the function test_codeml.write_ctl_file() is not called.  However,
        # test_codeml.write_ctl_file() must be called in order to write
        # 'self._rel_out_file' to use in test_codeml.run(), otherwise the code breaks.  
        #test_codeml.read_ctl_file('/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/parameters/codeml.ctl')
        test_codeml.write_ctl_file()
        test_codeml.run(verbose=True)
        read_codeml = read('codeml_results.txt')
        dN = calculate_average_dN(read_codeml)
        dNdS = calculate_average_omega(read_codeml)
        print "dN/dS is : ", dNdS
    except ZeroDivisionError:
        ##TODO: do smarter, more informative error handling from the Codeml module.
        dNdS = -1
    return (dNdS, dN)

def write_panorthologs(families, family_dir):
    for locus, data in families.iteritems():
        locus_outfile = open(family_dir + locus + ".fasta", "w")
        locus_outfile.write(data)
        locus_outfile.close()
    return None        

def main():
    panorthologs = "/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/hatcher_results/lerat-point7-reciprocal/lerat-point7-reciprocal.panorthologs"
    genedir = "/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/hatcher_results/nuc/"
    protdir = "/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/hatcher_results/proteins/"
    ##load the_families if it has been pickled, otherwise, create it anew (mainly for debugging).
    the_prot_families = {}
    if exists("prot_families.p"):
        prot_family_pickle = open("prot_families.p", "r")
        the_prot_families = pickle.load(prot_family_pickle)
    else:
        the_prot_families = read_panortholog_families(panorthologs, protdir, dna=False)
        prot_family_pickle = open("prot_families.p", "w")
        pickle.dump(the_prot_families, prot_family_pickle)
    ##pprint(the_prot_families)
    genehits = "/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/non_mutator_40K_diffs/genehits.tab"
    #the_prot_families = ltee_locustag_subset(the_prot_families, genehits)
    #Write the relevant panortholog families to file.
    dna_family_dir = "/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/panortholog_families/"
    prot_family_dir = "/Users/Rohan/Desktop/Projects/ltee_comparative_genomics/data/protein_panortholog_families/"
    #write_panorthologs(the_prot_families, prot_family_dir)
    ## omegadict indexes average dN/dS by locus id.
    omegadict = {}
    ## ditto for dNdict.
    dNdict = {}
    ##load the omegadict if it has been pickled, otherwise create it anew (mainly for debugging).
    if exists("omega.p"):
        omega_pickle = open("omega.p", "r")
        omegadict = pickle.load(omega_pickle)
    if exists("dN.p"):
        dN_pickle = open("dN/p", "r")
        dNdict = pickle.load(dN_pickle)
    for prot_panortholog_input_file in listdir(prot_family_dir):
        locustag, dot, ext = prot_panortholog_input_file.partition('.')
            
        ## there's no need to recreate existing .phy and .tre files.
        if not exists(locustag + '.phy'):
            print "current file: ", prot_panortholog_input_file
            align_using_sate(prot_family_dir + prot_panortholog_input_file, "protein")
            cur_satealn_filename = locustag + '.marker001.' + locustag + '.aln'
            dna_lookup_filename = dna_family_dir + locustag + '.fasta'
            back_translated_aln_filename = 'btrans_' + locustag + '.fasta'
            temp_phylip_filename = 'temp_' + locustag + '.phy'
            backtranslate(cur_satealn_filename, dna_lookup_filename, back_translated_aln_filename)
            convert_fasta_to_phylip(back_translated_aln_filename, temp_phylip_filename)
            old_phylip_file = open(temp_phylip_filename, "r")
            correct_phylip_file = open(locustag + '.phy', "w")
            correct_phylip_format(old_phylip_file, correct_phylip_file)
            correct_phylip_file.close()
            remove(temp_phylip_filename)
            omegadict[locustag], dNdict[locustag] = use_codeml(locustag)
            
        else:
            ## the .phy and .tre files exist.
            ## run codeml on loci which don't have a defined dN/dS value.
            if locustag not in omegadict.keys() or omegadict[locustag] == -1:
                omegadict[locustag], dNdict[locustag] = use_codeml(locustag)
                
        omega_pickle = open("omega.p", "w")
        pickle.dump(omegadict, omega_pickle)
        dN_pickle = open("dN.p", "w")
        pickle.dump(dNdict, dN_pickle)
    
    ## Now, write the output file.
    with open("/Users/Rohan/Desktop/gene_rates.csv", "wb") as outputfile:
        writer = csv.writer(outputfile)
        ltee_results = open(genehits)
        ltee_dict = {}
        for i, line in enumerate(ltee_results):
            if i == 0:
                next
            else:
                data = line.split()
                cur_locus = data[0]
                nonsynonymous = data[1]
                synonymous = data[2]
                ltee_dict[cur_locus] = (nonsynonymous, synonymous)
        header = ['locus_tag', 'nonsynonymous', 'synonymous', 'dNdS', 'dN']
        writer.writerow(header)
        for locus in the_prot_families.keys():
            row = [locus, 0, 0, None, None]
            try:
                (non, synon) = ltee_dict[locus]
                row = [locus, non, synon, None, None]
            except KeyError:
                row = [locus, 0, 0, None, None]
            try:
                row[-2] = omegadict[locus]
                row[-1] = dNdict[locus]
            except KeyError:
                row[-2] = -1
                row[-1] = -1
            writer.writerow(row)
            
if __name__ == '__main__':
    main()

##TODO:

## 1) I need to break down the separate parts of the pipeline into functions,
## simply to make debugging easier, and to make the separate bits faster.
## I should make this module act more like a utility, where each different type of call
## does a different part of the analysis.
## 2) refactor code that calculates dN/dS (return dN and dS separately).
## 3) write a 'pickler' function to load/save the two relevant data structures if needed.    
## 4) I need to choose the appropriate model in PAML to use for calculating dNdS for a gene family.
## 5) I need to divide up the analyses based on metadata on these E. coli.
## Rich commented that I will find 4 categories of E. coli:
## 1) Unknown, 2) Pathogens 3) Commensals 4) Environmental isolates.
## I should focus on pathogens and commensals: perhaps I can find subcategories
## 6) My x-axis should be the number of lines that have one or more mutations (1-6).
## 7) It is as yet unclear how I should treat genes which have non-SNP mutations in the LTEE.
## 8) I get these 'refused to clean temp files' type warnings from Sate. The program seems to run okay,
## but I'd like to get rid of these messages in the future.
## 9) this module should be able to do a sensible cleanup upon exit. Perhaps bin temp files in different ways,
## and offer an option to remove these files upon exit?
## 10) Catch errors returned by codeml.
## 11) Investigate dN/dS outliers : for instance, ECB__00188 has an apparent dN/dS of 73.

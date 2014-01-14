#!/usr/bin/python

#HMMERanalysis.py by Rohan Maddamsetti.
##This module contains useful wrappers for running HMMER,
##and a function to calculate a threshold to cut off HMMER results.

from operator import itemgetter
import subprocess
import cPickle
from os.path import isfile
from Bio import AlignIO

#constant used for chucking subprocess output.
IGNORE = open("/dev/null")

#TODO:
# I'm not sure how/where to handle the file cleanup yet.
def cleanup():
    pass

##Input: a seqfile string, and a seqdb string input for phmmer, and the name
##of the output file. This function outputs a parseable table.
##Output: the output file in the working directory.
def run_phmmer(seqfile, seqdb, output):
    phmmer_proc = subprocess.Popen(["phmmer", "--tblout", output,
                                    seqfile, seqdb], stdout=IGNORE,
                                   stderr=IGNORE)
    phmmer_proc.wait()
    return None

##Input: <hmmfile> output filename, <msafile> input filename
##Ouput: the output file in the working directory.
def run_hmmbuild(hmmfile, msafile):
    #make an empty file with the name contained in hmmfile.
    out = open(hmmfile, 'w')
    out.write('')
    out.close()
    args = 'hmmbuild ' +  hmmfile + " " + msafile
#    print args
    hmmbuild_proc = subprocess.Popen(args, shell=True)
    hmmbuild_proc.wait()
    return None

##Input: <hmmfile> query filename, <seqdb> search set filename, and the name
##of the output file. This function outputs a parseable table.
##Output: the output file in the working directory.
def run_hmmsearch(hmmfile, seqdb, output):
    hmmsearch_proc = subprocess.Popen(["hmmsearch", "--tblout", output,
                                   hmmfile, seqdb], stdout=IGNORE,
                                  stderr=IGNORE)
    hmmsearch_proc.wait()
    return None

# The following commands produce the necessary output for
# the function 'sort_hmmer_scores.'
# hmmsearch -option <f> <hmmfile> <seqdb>
# hmmsearch --tblout <gene>.txt  <gene>-hmm.sto allproteomes.txt
##Input: the name of a table of HMMER results.
# Results need to be ranked by Score NOT E-value.
##Output: a sorted list of tuples, each pair being a name and a score.
def sort_hmmer_scores(hmmer_table):
    hmmer_output = open(hmmer_table)
    score_dict = {}
    for line in hmmer_output:
        if line.startswith("#"):
            continue
        line_data = line.split()
        #print line_data
        name = line_data[0]
        score = float(line_data[5])
        score_dict[name] = score

    #Sort the dictionary based on its values (the scores)
    sorted_pairs = sorted(score_dict.items(), key=itemgetter(1))
    sorted_pairs.reverse() #put into descending order.
    return sorted_pairs

##Input: The sorted list of tuples made by sort_hmmer_scores,
##and a score cutoff.
##Output: The name of the threshold, based a score threshold.
##Returns None if all hits pass.
def find_score_threshold(sorted_pairs, score_cutoff):
    for items in sorted_pairs:
        cur_name, cur_score = items
        if cur_score < score_cutoff:
            return cur_name
    return None

##Input: The sorted list of tuples made by sort_hmmer_scores.
##Output: The name of the threshold, based on my double delta heuristic.
def find_inclusion_threshold(sorted_pairs):
    names = []
    scores = []
    for items in sorted_pairs:
        names.append(items[0])
        scores.append(items[1])

    #print names
    #print scores

    delta = []
    delta_names = []
    for i in range(len(scores)-1):
        delta.append(float(scores[i]) - float(scores[i+1]))
        delta_names.append("("+names[i]+")" + " - " + "(" + names[i+1]+")")

    for i in delta: #assert that the scores are truly in descending order.
        assert i >= 0
    #print delta
    #print delta_names

    double_delta = []
    double_delta_names = []
    for i in range(len(delta)-1):
        double_delta.append(delta[i] - delta[i+1])
        #This next line is atrocious and should be changed, but it works.
        double_delta_names.append(delta_names[i].partition('-')[2][2:-1])
    #print double_delta
    #print double_delta_names

    max_index = double_delta.index(max(double_delta))
    return double_delta_names[max_index]

#Memory-intensive! And it takes a long time to load 8GB into memory!
#Input: a list of UniProt flatfile format protein entries.
#a hash table, with sequence values keyed to Uniprot ID.
def make_proteome_hash(seqdb_file):
    seqdb = open(seqdb_file)
    proteome_hash = {}
    curr_id, curr_seq = ("","")
    for line in seqdb:
        if line.startswith("ID"):
            curr_id = line.split()[1]
        elif line.startswith("SQ"):
            line = seqdb.next()
            # until the break,
            while not line.startswith('//'):
                curr_seq += "".join(line.split())
                line = seqdb.next()
            proteome_hash[curr_id] = curr_seq
            curr_id, curr_seq = ("", "")
    return proteome_hash

def main():

    #print prot_hash
    phmmer_outfile = "phmmer.out"
    run_phmmer("ecoli_alkB.txt","allproteomes.txt", phmmer_outfile)
    phmmer_results = sort_hmmer_scores(phmmer_outfile)
    phmmer_threshold = find_inclusion_threshold(phmmer_results)
    print phmmer_threshold
##     ##truncate the results.
##     threshold_index = [k for k,v in phmmer_results].index(phmmer_threshold)
##     phmmer_tophits = phmmer_results[:threshold_index]

if __name__ == "__main__":
    main()

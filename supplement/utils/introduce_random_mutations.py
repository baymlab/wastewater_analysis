#!/usr/bin/python
import random
import sys
import argparse
from csv import writer

def main():
    parser = argparse.ArgumentParser(description="Introduce random mutations to original sequence")
    parser.add_argument('--file', dest = 'file', default = "data/SARS-CoV-2-NC_045513.fasta",required=False, type=str, help="file path to the sequence to be mutated")
    parser.add_argument('--nm', dest = 'nm', required=True, type=int, help="number of mutations to be randomly made")
    args = parser.parse_args()

    sequence = open (args.file, 'r').read()
    orseq = sequence
    sequence = "".join(sequence)

    start = 0
    end = len(sequence)-1
    prev_rndm_indx = []
    for i in range (0, args.nm):

        randomIndex = random.randint(start, end)

        while (randomIndex not in prev_rndm_indx) and (sequence[randomIndex] == "\n"):
            randomIndex = random.randint(start, end)

        prev_rndm_indx.append(randomIndex)
        
    
    
        curr_nucleotide_base = sequence[randomIndex]
        nucleotide_bases  = ['A', 'C', 'G', 'T']
        nucleotide_bases.remove(curr_nucleotide_base)
        mutation = random.choice(nucleotide_bases)
        sequence_list = list(sequence)
        sequence_list[randomIndex] = mutation 
        sequence = "".join(sequence_list)
        
    #print("Edit distance matches requirement = ", editdistance.eval(sequence, orseq) == args.nm )
    outdir = "Sequence_Similarity_Experiment/sequences/"

    f = open(outdir + "sequence_with_" + str(args.nm) + "_mutations" + ".fasta", "w")
    f.write(">hCoV-19/org_sequence_"+ str(args.nm)+"_rm/"+ str(args.nm) +"\n"+ sequence)
    f.close()
    
    return

if __name__ == "__main__":
    sys.exit(main())







    






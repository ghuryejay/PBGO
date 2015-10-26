from Bio import SeqIO
import argparse
import operator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--pacbio_reads",help="input pacbio long reads")
    parser.add_argument("-k","--kmer_size",help="kmer_size")
    parser.add_argument("-o","--output_file",help="output file")

    args = parser.parse_args()

    k = int(args.kmer_size)

    handle = open(args.pacbio_reads,"rU")

    kmer_map = {}

    for record in SeqIO.parse(handle,"fasta"):
        id = record.id
        s = record.seq.lower()
       # print s
        for i in range(0,len(s)-k+1):
            kmer = s[i:i+k]
            if kmer in kmer_map:
                kmer_map[str(kmer)] += 1
            else:
                kmer_map[str(kmer)] = 1

    #sorted_x = sorted(kmer_map.items(), key=operator.itemgetter(1))

    ofile = open(args.output_file,"w")
    for elem in kmer_map:
        ofile.write(elem+"\t"+str(kmer_map[elem])+"\n")

if __name__ == '__main__':
    main()

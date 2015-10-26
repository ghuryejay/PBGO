from Bio import SeqIO
import argparse
import networkx as nx

revcompl = lambda x: ''.join(
    [{'A': 'T',
     'C': 'G',
      'G': 'C',
      'T': 'A',
      'N': 'N',
      'R': 'N',
      'M': 'N',
      'Y': 'N',
      'S': 'N',
      'W': 'N',
      'K': 'N',
      'a': 't',
      'c': 'g',
      'g': 'c',
      't': 'a',
      'n': 'n',
      ' ': '',
      }[B] for B in x][::-1])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--pacbio_reads",help="input pacbio long reads")
    parser.add_argument("-k","--kmer_size",help="kmer_size")
    parser.add_argument("-c","--kmer_count",help="file containing list of kmers")
    parser.add_argument("-f","--filter",help="kmer cutoff filter")
    parser.add_argument("-o","--output_file",help="graph to be output")

    args = parser.parse_args()

    kmer_map = {}

    #load all kmercouts
    kmer_file = open(args.kmer_count,"r")
    lines = kmer_file.readlines()
    for line in lines:
        attrs = line.split("\t")
        attrs[-1] = attrs[-1][:-1]
        kmer_map[attrs[0]] = int(attrs[-1])


    print kmer_map
    
    fil = int(args.filter)
    k = int(args.kmer_size)

    '''
    For each read, examine each kmer, if that kmer is valid to be added,
    add it go graph. If its not in order of previous kmer, find out how
    difference between start of current kmer and previous kmer

    '''
    G = nx.DiGraph()
    handle = open(args.pacbio_reads,"rU")
    for record in SeqIO.parse(handle,"fasta"):
         s = record.seq.lower()
         skip_count = 1
         prev_kmer = s[0:k]
         if prev_kmer in kmer_map and kmer_map[prev_kmer] < fil:
            prev_kmer = ""
         for i in range(1,len(s) - k + 1):
            curr_kmer = s[i:i+k]
            if curr_kmer in kmer_map and kmer_map[curr_kmer] >= fil:
                if prev_kmer == "":
                  prev_kmer = curr_kmer
                  continue
                else:
                  G.add_edge(prev_kmer,curr_kmer,weight=skip_count)
                  G.add_edge(revcompl(prev_kmer),revcompl(curr_kmer),weight=skip_count)
                  skip_count = 1
                  prev_kmer = curr_kmer
            else:
                skip_count += 1 


    nx.write_gml(G,args.output_file)
    '''
    nx.draw(G)
    plt.show()
    '''


if __name__ == '__main__':
    main()

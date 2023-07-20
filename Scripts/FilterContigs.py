#!/usr/bin/env python3

import os
import argparse
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("contigFile", help="assign file")

    parser.add_argument("summaryFile", help="summary file")

    parser.add_argument("outputDir", help="output directory")

    parser.add_argument("filtContFile", help="filtered contig file")

    args = parser.parse_args()

    contig_comp = {}
#    import ipdb; ipdb.set_trace()
    bFirst = True
    with open(args.summaryFile) as f:
        for line in f:
            if not bFirst:
                toks = line.rstrip().split("\t")
                contig_comp[toks[0]] = int(toks[2].strip('%'))
            else:
                bFirst = False
    
    keep = []
    for header,seq in sfp(open(args.contigFile)):
        bFilt = False
        if header in contig_comp and header[-1] == 'c':
            if contig_comp[header] > 90:
                with open("%s/circ_%s.fa"%(args.outputDir,header),"w") as handle:
                    handle.write(">%s\n%s\n"%(header,seq))
                    bFilt = True
        if not bFilt:
            keep.append(">%s\n%s\n"%(header,seq))
    with open(args.filtContFile,'w') as f:
        for line in keep:
            f.write(line)
if __name__ == "__main__":
    main(sys.argv[1:])


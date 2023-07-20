#!/usr/bin/env python3

import os
import argparse
import sys

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("inputFile", help="assign file")

    args = parser.parse_args()
#ctg10000l       1314    1695    .       6.8     +       Prodigal_v2.6.3 CDS     0       ID=347_3;pa
#rtial=00;start_type=TTG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.661;conf=82.58;score=6.77;c
#score=11.35;sscore=-4.58;rscore=8.88;uscore=-1.87;tscore=-11.60;

    with open(args.inputFile) as f:
        for line in f:
            toks = line.rstrip().split("\t")
            
            process = toks[9]
            ptoks = process.split(';')
            idx = (ptoks[0].split('_'))[-1]
            
            print('%s\t%s\t%s\t%s' % (toks[0],toks[1],toks[2],toks[0] + '_' + idx))
if __name__ == "__main__":
    main(sys.argv[1:])


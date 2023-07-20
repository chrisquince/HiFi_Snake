#!/usr/bin/env python3
import numpy as np
from scipy.special import gammaln
import itertools
import networkx as nx

from collections import defaultdict
from collections import Counter

import sys
import argparse
import re 



def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("gfaFile", help="unitig graph file")
    
   # parser.add_argument("unitig_len", help="unitig len file")
    
    parser.add_argument('paf_files', metavar='N', type=str, nargs='+',
                    help='gaf format files')
        
    args = parser.parse_args()

    #import ipdb;ipdb.set_trace()

    directedUnitigBiGraph = nx.DiGraph()

    overlapLengths = defaultdict(dict) 
    
    unitigLens = {}

    unitigs = []
    
    
    with open(args.gfaFile) as f:
        for line in f:
        
            line = line.rstrip()
        
            toks=line.split('\t')
    
            if toks[0] == 'L':
    
                #L	s115.ctg091013l	+	s115.ctg057300l	-	4149M	L1:i:7043
    
            
                nodeOutName = toks[1]
                
                nodeStart = toks[2]
                
                nodeInName = toks[3]
                
                nodeEnd = toks[4]
                
                nodeNameOut = nodeOutName + nodeStart

                nodeNameEnd = nodeInName + nodeEnd
                
                oL = int(toks[5][:-1])
                
                
                directedUnitigBiGraph.add_edge(nodeNameOut,nodeNameEnd,overlap=oL)
                
                reverseStart = '+'
                if nodeStart == '+':
                    reverseStart = '-'
                
                
                reverseEnd = '+'
                if nodeEnd == '+':
                    reverseEnd = '-'
                
                
                rNodeNameOut = nodeOutName + reverseStart

                rNodeNameEnd = nodeInName + reverseEnd

                
                directedUnitigBiGraph.add_edge(rNodeNameEnd,rNodeNameOut,overlap=oL)  
                
                      
                
                overlapLengths[nodeNameOut][nodeNameEnd] = oL
                overlapLengths[rNodeNameEnd][rNodeNameOut] = oL
            
            elif toks[0] == 'S':    
        
                unitigName = toks[1]
        
                unitigLens[unitigName] = len(toks[2])
        
                unitigs.append(unitigName)
        
                nodePlusName  = unitigName + '+'
                nodeMinusName = unitigName + '-'
            
                directedUnitigBiGraph.add_node(nodePlusName)
                directedUnitigBiGraph.add_node(nodeMinusName)
    
    
    unDirectedUnitigBiGraph = directedUnitigBiGraph.to_undirected()
    
    NC = nx.number_connected_components(unDirectedUnitigBiGraph)
 

    adjLens = {}

    unitigPos = {}
#import ipdb; ipdb.set_trace()
    
    for unitig in unitigs:
    
        unitigp = unitig + '+'
        
        unitign = unitig + '-'
        
        length = unitigLens[unitig]
        
        if len(overlapLengths[unitigp].values()) > 0:
            pOut = min(overlapLengths[unitigp].values())
        else:
            pOut = 0
        
        if len(overlapLengths[unitign].values()) > 0:
            nOut = min(overlapLengths[unitign].values())
        else:
            nOut = 0
        
        adjL = length - pOut - nOut 
        
        #print(unitig + ',' + str(adjL))
        
        #assert adjL > 0
        
        if adjL < 0:
            adjL = 0
            
            adjLens[unitig] = adjL
            unitigPos[unitig] = (int(length/2),int(length/2)) #defined as exclusive of ends
        else:
            adjLens[unitig] = adjL
            unitigPos[unitig] = (nOut,length - pOut) #defined as exclusive of ends
        
        
        #print(unitig + ',' + str(adjL))

    #import ipdb; ipdb.set_trace()     
    
         

#m54313U_200724_221921/67234/ccs	5271	0	5271	+	>s195.ctg008942l	122576	51833	57104	5263	5271	255	NM:i:8	dv:f:0.00151774	id:f:0.998482	cg:Z:133M1M1993M1M488M1M2186M1M36M1M143M1M14M1M160M1M110M	0.998482



#1	string	Query sequence name
#2	int	Query sequence length
#3	int	Query start (0-based; BED-like; closed)
#4	int	Query end (0-based; BED-like; open)
#5	char	Relative strand: "+" or "-"
#6	string	Target sequence name
#7	int	Target sequence length
#8	int	Target start on original strand (0-based)
#9	int	Target end on original strand (0-based)
#10	int	Number of residue matches
#11	int	Alignment block length
#12	int	Mapping quality (0-255; 255 for missing)

    
    unitigHitsAll = {}

    contigFreq = defaultdict(Counter)
    gsamples = []
    contigs = set()
    reads = {}
    
    read_hits_all = defaultdict(list)
    for gfile in args.paf_files:
        with open(gfile) as f:
            
            unitigHits = {}
            
            for unitig, lengthu in unitigLens.items():
                unitigHits[unitig] = np.zeros(lengthu,dtype=np.int)
        
            read_hits = set()
            
            gstub = gfile[:-4]
            gsamples.append(gstub)
            
            for line in f:
        
                line = line.rstrip()
            
                toks = line.split('\t')


                bFilter = True
                
                read_id = toks[0]
                read_len = int(toks[1])
                    
                nMatch = int(toks[9])
                nAlign = int(toks[10])
                
                         
                contig = toks[5]
                
                qMatch = int(toks[3]) - int(toks[2])
                
                fm = qMatch/read_len
                
                pid = nMatch/nAlign

                qid = fm*pid
                
                bFilter = qid < 0.80
                
                
                #print(str(bFilter) + ',' + str(pid) + ',' + str(fm))
                if not bFilter:
                
                    sStart = int(toks[7])
                    sEnd   = int(toks[8])
                    if contig in unitigPos:
                        (nStart,nEnd)  = unitigPos[contig]  
                
                        if sEnd >= nStart and sStart <= nEnd: 
                
                            read_hits_all[read_id].append((contig,nMatch,nAlign,pid))
                
                if not bFilter and read_id not in read_hits:
                    
                    read_hits.add(read_id)
                    
                    sStart = int(toks[7])
                    sEnd   = int(toks[8])
    
                    if contig in unitigHits:
                        unitigHits[contig][sStart:sEnd] += np.ones(sEnd - sStart,dtype=np.int)
                    
                        contigs.add(contig)
            
            unitigHitsAll[gstub] = unitigHits     

    nR = 0
    nA = 0
    for (read_id,hits) in read_hits_all.items():
        
        lHits = len(hits)
        nR += 1
        if lHits > 1:
        
        
            diff = hits[0][1] - hits[1][1]
            fd = diff/hits[0][2]
        
            if diff < 3:
                nA += 1
                       
                #print(read_id + ',' + str(diff) + ',' + str(fd))
        
    
           

    gString = ','.join(gsamples)
    
    print('Freq,' + gString)
    clist = list(unitigLens.keys())
    clist.sort()
    for contig in clist:
        covs = []
    
        for g in gsamples:
            uHits = unitigHitsAll[g][contig]
        
            (uStart, uEnd) = unitigPos[contig]
            
            if adjLens[contig] > 0:
            
                covEst = np.mean(uHits[uStart:uEnd])
            else:
                umid = int(0.5*unitigLens[contig])

                covEst = uHits [umid] 
                
            covs.append(str(round(covEst,1)))
    
        cString =  ','.join(covs)
        print ("%s,%s"  % (contig, cString))
        
        
    
if __name__ == "__main__":
    main(sys.argv[1:])     

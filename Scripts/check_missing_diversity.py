from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from collections import defaultdict,Counter
import numpy as np
import glob

# mdbg
ROOT = "/home/sebr/seb/Project/temp/HIFI_annotation/plankton_mdbg"

contig_to_cog = defaultdict(list)
cog_to_contig = defaultdict(list)
for line in open("%s/annotation/SCG_cluster_0.99.tsv"%ROOT):
    cog = "_".join(line.rstrip().split("\t")[0:2])
    for orf in line.rstrip().split("\t")[2:]:
    	contig = "_".join(orf.split("_")[:-1])
        contig_to_cog[contig].append(cog)
        cog_to_orfs[cog].append(contig)


cogs_diversity = Counter([el.split("_")[0] for el in cog_to_orfs])
initial_diversity = np.median(cogs_diversity.values())

mags = glob.glob("%s/MAGs_gaetan/mags/*.fa"%ROOT)
mag_contigs = {header for file in mags for header,seq in sfp(open(file))}
cog_magged = {cog for contig in mag_contigs if contig in contig_to_cog for cog in contig_to_cog[contig]}
cog_un_magged = set(cog_to_orfs.keys())-cog_magged
un_magged_div = np.median(Counter([el.split("_")[0] for el in cog_un_magged]).values())

# hifiasm
ROOT = "/mnt/gpfs/seb/Project/temp/HIFI_annotation/Plankton_ctg"

contig_to_cog = defaultdict(list)
cog_to_contig = defaultdict(list)
for line in open("%s/annotation/SCG_cluster_0.99.tsv"%ROOT):
    cog = "_".join(line.rstrip().split("\t")[0:2])
    for orf in line.rstrip().split("\t")[2:]:
    	contig = "_".join(orf.split("_")[:-1])
        contig_to_cog[contig].append(cog)
        cog_to_orfs[cog].append(contig)


cogs_diversity = Counter([el.split("_")[0] for el in cog_to_orfs])
initial_diversity = np.median(cogs_diversity.values())

mags = list(glob.glob("%s/MAGs/mags/*.fa"%ROOT))
for line in open("%s/binning/metabat2/metabat2_MAG_list.txt"%ROOT):
	mags.append("%s/binning/metabat2/bins/Bin_%s.fa"%(ROOT,line.rstrip()))

import pandas
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool

path_to_gtf = "/home/ubuntu/gencode.v37.annotation.gtf"

gtf = open(path_to_gtf, 'r').read().split("\n")

#### Remove first 4 lines 
gtf = gtf[5:]


######## Generate BED12 file from GTF file
bed12 = {}

#### Create fields with needed information
keys = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
    "blockSizes", "blockStarts"]

for key in keys:
    bed12[key]=[]

#### Append required entries:
list_of_transcriptids = []

for line in gtf:
    line=line.split("\t")
    if len(line) >=2 and line[2] == "transcript":
        transcript_id = line[8].split(";")[1].split(" ")[2][1:-1]
        if transcript_id not in bed12["name"]:
            bed12["chrom"].append(line[0])  # chromosome
            bed12["chromStart"].append(int(line[3])-1) # chromosome Start
            bed12["chromEnd"].append(int(line[4])) # chromosome End
            bed12["name"].append(transcript_id) # Transcript ID
            bed12["score"].append(0) # Score
            bed12["strand"].append(line[6]) # Strand
            bed12["thickStart"].append(int(line[3])-1) # thick Start = chromosome Start
            bed12["thickEnd"].append(int(line[4])) # thick End = chromosome End
            bed12["itemRgb"].append(0) # RGB 
        list_of_transcriptids.append(transcript_id)  


abridged_gtf = pandas.DataFrame(gtf)
abridged_gtf =  abridged_gtf[0].str.split('\t', expand= True)
exon_frame = abridged_gtf[abridged_gtf[2] == 'exon']

for i in range(len(exon_frame[8])):
    exon_frame.iat[i,8] = exon_frame.iat[i,8].split(";")[1].split(" ")[2][1:-1]

def exon_info(transcriptid):
    exon_transcript_id = exon_frame[exon_frame[8] == transcriptid]
    exon_count = len(exon_transcript_id)
    chromStart = bed12["chromStart"][bed12["name"].index(transcriptid)]
    exon_transcript_id = exon_transcript_id.astype({3: 'int', 4: 'int'})
    blockstarts = list(exon_transcript_id[3] - chromStart - 1)
    blocksizes = list(exon_transcript_id[4] - exon_transcript_id[3] + 1)
    blocksizes = [x for _,x in sorted(zip(blockstarts, blocksizes))]
    blockstarts = sorted(blockstarts)
    blocksizes = ','.join([str(item) for item in blocksizes])
    blockstarts = ','.join([str(item) for item in blockstarts])
    return exon_count, blocksizes, blockstarts


if __name__ == '__main__':
    pool=Pool()
    bed12["blockCount"] = [0] * len(bed12["name"])
    bed12["blockSizes"] = [0] * len(bed12["name"])
    bed12["blockStarts"] = [0] * len(bed12["name"])
    x=pool.map(exon_info, tqdm(list_of_transcriptids))
    for tr_index in tqdm(range(len(list_of_transcriptids))):
        bed12["blockCount"][bed12["name"].index(list_of_transcriptids[tr_index])] =  x[tr_index][0] # number of exons
        bed12["blockSizes"][bed12["name"].index(list_of_transcriptids[tr_index])] = x[tr_index][1] # size of exons
        bed12["blockStarts"][bed12["name"].index(list_of_transcriptids[tr_index])] = x[tr_index][2] # start of exons relative to TSS


#### Check if 0s are absent from blockCounts, blockSizes and blockStarts
if 0 in bed12["blockCount"]:
    print("There are zeros in exon counts. Please check again.")

if 0 in bed12["blockSizes"]:
    print("There are zeros in exon sizes. Please check again.")

if 0 in bed12["blockStarts"]:
    print("There are zeros in exon start positions. Please check again.")

#### Write to BED file
df = pandas.DataFrame(bed12)
df.to_csv("/home/ubuntu/gencode.v37.annotation.bed", index=False, sep="\t", header=False, mode='a')





# summarizeOntDeletions
Simple script for cleaning up amplicon ONT-seq data for further analysis of deletions sizes. Takes amplicon ONT-seq datafiles as FASTQ.GZ and a sample definition file & outputs a simple QC summary, amplicon length summary, and alignments to the target genome. 

The output of this preprocessing script is intended for use in further analysis or plotting (e.g., GraphPad Prism).
<br/><br/>

## description of main steps
FASTQ files are processed by 3 rounds of cutadapt to filter:
1. end-to-end intact amplicons
2. amplicons specific to the target region
3. maximum allowed length.

Reads passing all filtering steps are aligned to the human genome (GRCh38) using minimap2, summarized using deepTools, and the top aligning regions are reported. The majority (>95%) of reads passing filters should align to the target genomic region.
<br/><br/>

## main outputs
A summary table the output subdirectory `countLengths` for each input FASTQ file, specifying columns:
1. read lengths
2. read count
3. fraction of total
4. fraction of reads greater than read length (1)

A QC summary of reads passing each filtering step & genomic region alignment frequency is given in the `processingSummary.txt` file.
<br/><br/>

## setup conda environment
A conda environment with the required packages can be created with:

```
conda env create -f environment.yml
```

## usage 
```
conda activate summarizeOntDels
cd <outputDirectory>

./summarizeOntAmplicons.sh \
<pathToNgsData> \
<pathToSampleInfoFile> \
<pathToMiniMap2GenomeIndex>
```

## preparation of the sample definition file
An example sample definition file (`sampleInfoExample.txt`) is provided & corresponds to the ONT-seq data from K562 edited at HBB with CRISPR/Cas9 (`./example/input/barcode02.fastq.gz`).

The sample definition file is a 7-column text file specifying:

#### sample information
1.	sample name for output
2.	input FASTQ.GZ file root name (e.g., `barcode02.fastq.gz` becomes `barcode02`)

#### cutadapt round 1 (see below figure)
Sequences for filtering via full length PCR primers at both ends.

3.	forward strand, cutadapt pattern: `[ONTstubFwd][PCR primer F]…[rc_PCR primer R][rc_ONTstubRev]`
4.	reverse strand, cutadapt pattern: `[ONTstubRev][PCR primer R]…[rc_PCR primer F][rc_ONTstubFwd]`

#### cutadapt round 2 (see below figure)
Sequences for filtering genomic loci specificity via first and last 30bps of the gDNA insert:

5.	forward strand, cutadapt pattern: `[fwd_first30bps]..[fwd_last30bps]`
6.	reverse strand, cutadapt pattern: `[rev_first30bps]…[prev_last30bps]`

#### maximum amplicon length
7.	amplicon length expected for unedited cells
<br/><br/>

#### Expected structure of full-length ONT-seq read and resulting cutadapt filtering sequences:
![figure 1a](https://github.com/cornlab/summarizeOntDeletions/blob/main/misc/figure1.PNG?raw=true)
<br/><br/>

## example for sample definition setup
For the datafile in example `barcode02.fastq.gz`, the HBB locus was PCR amplified using the following 2 primers. **Lowercase sequence is ONT-seq PCR stub sequence; uppercase is specific to target gene.**

-	PCR F primer: 	`5’-tttctgttggtgctgatattgcTCAAGCTACAAAAAGCCGCC`
-	PCR R primer: 	`5’-acttgcctgtcgctctatcttcCCTTGAAGCCAGGATGATGGT`
<br/><br/>

The cutadapt **round 1** filtering sequences are:
- ONT-seq forward strand read: [PCR F primer]…[revComp PCR R primer]
```
tttctgttggtgctgatattgcTCAAGCTACAAAAAGCCGCC...ACCATCATCCTGGCTTCAAGGgaagatagagcgacaggcaagt
```
- ONT-seq reverse strand read: [PCR R primer]…[revComp PCR F primer]
```
acttgcctgtcgctctatcttcCCTTGAAGCCAGGATGATGGT...GGCGGCTTTTTGTAGCTTGAgcaatatcagcaccaacagaaa
```
The cutadapt **round 2** filtering sequences are simply the 30bps immediately adjacent to these primer binding regions for the respective forward and reverse strand ONT-seq reads.


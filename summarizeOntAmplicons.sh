#!/bin/bash

# 24 January 2024
# Charles D Yeh
#
#
# limited processing of ONT-seq data by filtering using cutadapt with custom config
# ....cutadapt 1	simple filtering & trimming for primer pairs (full primer)
# ....cutadapt 2	simple filtering for target gDNA sequence (30bp each end)
# ....cutadapt 3	simple filtering for maximum amplicon size (ref genome size +20bps ($maxLengthExtraBp))
#
# simple count of read lengths thereafter, with no consideration of alignment to target site
# minimap2 align reads to hg38; default ONT-seq configuration for long read genomic DNA
#
# required packages (conda env with: cutadapt, minimap2, samtools, deepTools, bioawk)

# USAGE
# activate conda environment with appropriate dependencies.
# runwith "./summarizeOntAmplicons.sh <pathToNgsData> <pathToSampleInfoFile> <pathToAlignerGenome>"


# load input params
inputDir=$1
sampleInfo=$2
genome=$3

maxLengthExtraBp=20
maxThreads=30


# print basic information to console
printf "\n"
date
printf "\n"

conda --version
printf "conda environment: \t$CONDA_DEFAULT_ENV\n"
printf "current directory: \t$PWD\n"
printf "NGS data directory: \t./${inputDir}\n"
printf "sample information: \t./${sampleInfo}\n"
printf "aligner genome: \t./${genome}\n\n"


# directory substructure
fastqTrimmed="cutadapt/1_trimPrimer"
fastqFiltered="cutadapt/2_filterGenomicTarget"
fastqLengthCutoff="cutadapt/3_lengthCutoff"
cutadaptReports="cutadapt/reports"

outputLengths="countLengths"
outputAlignFail="minimap2_failFilter2_genomicTarget"
outputAlignPass="minimap2_passFilter3_maxLength"


# make directories
mkdir -p ${fastqTrimmed}
mkdir -p ${fastqFiltered}
mkdir -p ${fastqLengthCutoff}
mkdir -p ${cutadaptReports}

mkdir -p ${outputLengths}
mkdir -p ${outputAlignFail}
mkdir -p ${outputAlignPass}


# main process
# FASTQ filtering for read integrity, genomic target, and maxLength
# ...countReadLength summarization
# ....minimap2 alignment to given genome
{ 
	read
	while IFS="	" read -r outputName inputFastqName ca1_fwdStrand ca1_revStrand ca2_fwdStrand ca2_revStrand maxLength;
	do
		
		#
		# 3 steps of cutadapt filtering
		
		cutadapt \
		--report=minimal \
		--quality-cutoff 0 \
		--minimum-length 100 \
		--overlap 200 \
		-e 0.2 \
		-j 12 \
		--times 1 \
		--action=trim \
		-g fwdStrand=${ca1_fwdStrand} \
		-g revStrand=${ca1_revStrand} \
		-o ./${fastqTrimmed}/${outputName}.trimmed.fastq.gz \
		--untrimmed-output ./${fastqTrimmed}/${outputName}.untrimmed.fastq.gz \
		${inputDir}/${inputFastqName}.fastq.gz \
		> ./${cutadaptReports}/${outputName}.cutadapt1.trim.log \
		\
		&& cutadapt \
		--report=minimal \
		--quality-cutoff 0 \
		--minimum-length 100 \
		--overlap 200 \
		-e 0.2 \
		-j 12 \
		--times 1 \
		--action=lowercase \
		-g fwdStrand=${ca2_fwdStrand} \
		-g revStrand=${ca2_revStrand} \
		-o ./${fastqFiltered}/${outputName}.passTargetFilter.fastq.gz \
		--untrimmed-output ./${fastqFiltered}/${outputName}.failTargetFilter.fastq.gz \
		./${fastqTrimmed}/${outputName}.trimmed.fastq.gz \
		> ./${cutadaptReports}/${outputName}.cutadapt2.filter.log \
		\
		&& cutadapt \
		--report=minimal \
		--quality-cutoff 0 \
		--minimum-length 100 \
		--maximum-length $((maxLength+maxLengthExtraBp)) \
		-j 12 \
		--times 1 \
		-o ./${fastqLengthCutoff}/${outputName}.passMaxLength.fastq.gz \
		--too-long-output ./${fastqLengthCutoff}/${outputName}.failMaxLength.fastq.gz \
		./${fastqFiltered}/${outputName}.passTargetFilter.fastq.gz \
		> ./${cutadaptReports}/${outputName}.cutadapt3.maxLength.log \
		
		#
		# bioawk: generate summary counts of readLengths passing all 3 cutadapt filtering steps
		# awk: convert above bioawk readLength to (frac total reads) & (runCumSum of that)
		# output colunms: (ampliconSize, readCount, fracReads, runCumSum)
		
		bioawk -c fastx '{print length($seq)}' < ./${fastqLengthCutoff}/${outputName}.passMaxLength.fastq.gz \
		| sort -n | uniq -c > ${outputLengths}/${outputName}.lengths.temp \
		\
		&& printf "readLength\tcount\tfracTotal\tfracReadsGreater\n" \
		> ./${outputLengths}/${outputName}.lengths.txt \
		\
		&& awk 'BEGIN{OFS="\t"} (NR==FNR){sum+=$1; next} {rowFrac=($1/sum); fracRunCumSum+=($1/sum); print $2,$1,rowFrac, 1-fracRunCumSum}' \
		${outputLengths}/${outputName}.lengths.temp \
		${outputLengths}/${outputName}.lengths.temp \
		>> ./${outputLengths}/${outputName}.lengths.txt \
		\
		&& rm ${outputLengths}/${outputName}.lengths.temp
		
		#
		# minimap2 alignment of reads passing all 3 cutadapt filter phases to given genome
		# and calculate coverage (bamCoverage)
		
		minimap2 -ax map-ont -t${maxThreads} --rmq=yes \
		${genome} \
		./${fastqLengthCutoff}/${outputName}.passMaxLength.fastq.gz \
		2> ${outputAlignPass}/${outputName}.passMaxLength.mm2.log \
		\
		| samtools view -S -b - \
		| samtools sort -@ $((maxThreads/3))- -o ${outputAlignPass}/${outputName}.passMaxLength.bam \
		\
		&& samtools index ${outputAlignPass}/${outputName}.passMaxLength.bam \
		\
		&& bamCoverage --binSize 10000 \
		--numberOfProcessors 12 \
		--normalizeUsing CPM \
		--outFileFormat bedgraph \
		--bam ${outputAlignPass}/${outputName}.passMaxLength.bam \
		--outFileName ${outputAlignPass}/${outputName}.passMaxLength.bg
		
		#
		# minimap2 alignment of reads meeting ONLY these criteria
		# ...pass cutadapt 1	amplicon integrity
		# ...fail cutadapt 2	genomic target region
		# ...NULL cutadapt 3	max length NOT assessed
		
		minimap2 -ax map-ont -t${maxThreads} --rmq=yes \
		${genome} \
		./${fastqFiltered}/${outputName}.failTargetFilter.fastq.gz \
		2> ${outputAlignFail}/${outputName}.failTargetFilter.mm2.log \
		\
		| samtools view -S -b - \
		| samtools sort -@ $((maxThreads/3)) - -o ${outputAlignFail}/${outputName}.failTargetFilter.bam \
		\
		&& samtools index ${outputAlignFail}/${outputName}.failTargetFilter.bam \
		
	done } < ${sampleInfo} 


# summarize cutadapt reports into table of reads passing each level of filtering & top alignment region information
printf "outputName\tinputReads\treadsPassCutadapt1\treadsPassCutadapt2\treadsPassCutadapt3\ttopLoci1\ttopLoci1_cpm\ttopLoci2\ttopLoci2_cpm\ttopLoci3\ttopLoci3_cpm\n" \
> ./processingSummary.txt \

{ 
	read
	while IFS="	" read -r outputName inputFastqName ca1_fwdStrand ca1_revStrand ca2_fwdStrand ca2_revStrand maxLength;
	do 
	
		paste <(printf ${outputName}) <(awk 'BEGIN{OFS="\t"}(NR==2){print $2,$7}' ./${cutadaptReports}/${outputName}.cutadapt1.trim.log) \
		<(awk '(NR==2){print $7}' ./${cutadaptReports}/${outputName}.cutadapt2.filter.log) \
		<(awk '(NR==2){print $7}' ./${cutadaptReports}/${outputName}.cutadapt3.maxLength.log) \
		<(cat ${outputAlignPass}/${outputName}.passMaxLength.bg | sort -nrk 4 | head -3 | awk 'BEGIN{OFS="\t"} (NR==1){print $1":"$2"-"$3,$4}') \
		<(cat ${outputAlignPass}/${outputName}.passMaxLength.bg | sort -nrk 4 | head -3 | awk 'BEGIN{OFS="\t"} (NR==2){print $1":"$2"-"$3,$4}') \
		<(cat ${outputAlignPass}/${outputName}.passMaxLength.bg | sort -nrk 4 | head -3 | awk 'BEGIN{OFS="\t"} (NR==3){print $1":"$2"-"$3,$4}') \
		>> ./processingSummary.txt
		
	done } < ${sampleInfo}


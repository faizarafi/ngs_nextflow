#!/usr/bin/env nextflow 

process fastp_process {

//for fastp process, make directory fastp_res for outpus
//input val for sample_id and paths for the pair of reads (raw fastq)
//output val for sample_id and path to processed fastq
//in script, load fastp module and run fastp with given inputs and listed parameters to output pair of reads

publishDir "fastp_res"
input:
tuple val(sample_id), path(reads)
output: 
tuple val(sample_id), path("${sample_id}.{1,2}.fastq.gz"), emit: fastq

script:
"""
module load fastp/intel/0.20.1

fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}.1.fastq.gz -O ${sample_id}.2.fastq.gz \
--length_required 76 \
--n_base_limit 50 \
--detect_adapter_for_pe \
-h ${sample_id}.fastp.html -j ${sample_id}.fastp.json

"""
}

process align_process {

//for align process, make directory align_res for outputs
//input val for sample_id and paths for the pair of read (processed fastq), path to indexDir and val ref
//output val for sample_id and path to bam file
//in script, load modeules bwa and samtools, run bwa-mem with given inputs and header to output sam file
// run samtools view to convert the sam file to bam file 

publishDir "align_res"
input:
tuple val(sample_id), path(reads)
path (indexDir)
val (ref)
output:
tuple val(sample_id), path("${sample_id}.bam"), emit:bam

script:
"""
module load bwa/intel/0.7.17
module load samtools/intel/1.14

bwa mem \
-M \
-t $SLURM_CPUS_PER_TASK \
-R "@RG\\tID:${sample_id}.id\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}.lb" \
${indexDir}/${ref} ${reads[0]} ${reads[1]}  > ${sample_id}.sam \

samtools view -b -h ${sample_id}.sam > ${sample_id}.bam

"""
}

process sort_process{

//for sort_process, make directory sorted.bam_res for outputs
//input a val for sample id and path to the bam files 
//output a sorted.bam file 
//in script, load picard module and run picard SortSam with given input and defined output

publishDir "sorted.bam_res"
input: 
tuple val(sample_id), path("${sample_id}.bam")
output: 
path("${sample_id}_sorted.bam")

script:
"""
module load picard/2.17.11

java -Xmx15g -XX:ParallelGCThreads=1 -jar ${params.picard} SortSam \
INPUT=${sample_id}.bam \
OUTPUT=${sample_id}_sorted.bam \
SORT_ORDER=coordinate \
MAX_RECORDS_IN_RAM=10000000 \
VALIDATION_STRINGENCY=LENIENT

"""
}

workflow {

//workflow occurs as folllows:

//define file channel to get pair of raw fastq files
file_ch=Channel.fromFilePairs("$params.inputDir/*_{1,2}.filt.100k.fastq.gz")

//do fastp_process with inputs from file_ch channel
fastp_process(file_ch)

//define channel fastq_ch with outputs (processed fastq) from fastp_process
fastq_ch= fastp_process.out.fastq

//define channels index_ch and ref_ch with path to reference file and bwa-index directory
index_ch=Channel.fromPath(params.indexDir)
ref_ch=Channel.of (params.ref)

//do align_process with inputs from fastq, index and ref channels
align_process(fastq_ch, index_ch, ref_ch)

//define channel bam_ch with output (aligned bam files) from align_process
bam_ch = align_process.out.bam

//do sort_process with input from channel bam_ch
sort_process(bam_ch)

}

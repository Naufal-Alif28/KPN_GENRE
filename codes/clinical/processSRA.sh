#!/usr/bin/env bash

#SBATCH --job-name="nalif_01"
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

module load bioinformatics/samtools/1.21
module load bioinformatics/bowtie2/2.5.3
module load bioinformatics/miniforge3/24.3.0-0

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/mamba.sh
mamba activate kpenv

#set threads/core count
CORE=16
ADDCORE=15

#make output directories
mkdir -p "sra_fastq"
mkdir -p "reads_annotate"
#paired
mkdir -p "fastq_trimmed_paired"
mkdir -p "bam_merge_paired"
mkdir -p "maps_paired"
#single
mkdir -p "fastq_trimmed_single"
mkdir -p "bam_sort_single"
mkdir -p "maps_single"

#SRA list formatting
dos2unix "sra_list.txt"
dos2unix "sra_list_paired.txt"
dos2unix "sra_list_single.txt"

#paired reads processing
echo "START PROCESSING PAIRED"
while IFS= read -r sra; do
    echo "start ${sra}"

    #extract FASTQ
    echo "start download SRA: ${sra}"
    fasterq-dump "cache/sra/${sra}.sra" --split-3 -e ${CORE} -O "sra_fastq" 
    rm "cache/sra/${sra}.sra" #clear SRA cache after split to FASTQ
    echo "end download SRA: ${sra}"

    #processFASTQ (paired)
    echo "start process FASTQ: ${sra}"
    #trimming
    sickle pe -t sanger -q 25 -f "sra_fastq/${sra}_1.fastq" -r "sra_fastq/${sra}_2.fastq" -o "fastq_trimmed_paired/${sra}_1_trim.fastq" -p "fastq_trimmed_paired/${sra}_2_trim.fastq" -s "fastq_trimmed_paired/${sra}_orphan_trim.fastq"
    #paired reads
    #mapping to KPN genome
    bowtie2 -q --fr -x "bt2index/KPN" \
        -1 "fastq_trimmed_paired/${sra}_1_trim.fastq" \
        -2 "fastq_trimmed_paired/${sra}_2_trim.fastq" \
        -S "${sra}_paired.sam" -p ${CORE}
    #conversion to .bam
    samtools view --threads ${ADDCORE} -bS "${sra}_paired.sam" > "${sra}_paired.bam"
    #sorting
    samtools sort --threads ${ADDCORE} "${sra}_paired.bam" -o "${sra}_paired_sort.bam"
    #remove temporary unsorted .sam and .bam mapping
    rm "${sra}_paired.sam" "${sra}_paired.bam"
    #orphaned reads
    #mapping to KPN genome
    bowtie2 -q -x "bt2index/KPN" \
        -U "fastq_trimmed_paired/${sra}_orphan_trim.fastq" \
        -S "${sra}_orphan.sam" -p ${ADDCORE}
    #conversion to .bam
    samtools view --threads ${ADDCORE} -bS "${sra}_orphan.sam" > "${sra}_orphan.bam"
    #sorting
    samtools sort --threads ${ADDCORE} "${sra}_orphan.bam" -o "${sra}_orphan_sort.bam"
    #remove temporary unsorted .sam and .bam mapping
    rm "${sra}_orphan.sam" "${sra}_orphan.bam"
    #remove temporary FASTQ
    rm "fastq_trimmed_paired/${sra}_1_trim.fastq" "fastq_trimmed_paired/${sra}_2_trim.fastq" "fastq_trimmed_paired/${sra}_orphan_trim.fastq" "sra_fastq/${sra}_1.fastq" "sra_fastq/${sra}_2.fastq"
    #merge alignments, remove duplicates
    samtools merge --threads ${ADDCORE} "${sra}_merge.bam" "${sra}_paired_sort.bam" "${sra}_orphan_sort.bam"
    samtools sort --threads ${ADDCORE} "${sra}_merge.bam" -o "bam_merge_paired/${sra}_merge_sort.bam"
    rm "${sra}_merge.bam" "${sra}_paired_sort.bam" 
    #conversion to tsv
    samtools index --threads ${ADDCORE} "bam_merge_paired/${sra}_merge_sort.bam" -b "bam_merge_paired/${sra}_merge_sort.bam.bai"
    samtools idxstats --threads ${ADDCORE} "bam_merge_paired/${sra}_merge_sort.bam" > "maps_paired/${sra}_map.tsv"
    rm "${sra}_orphan_sort.bam"
    echo "end process FASTQ: ${sra}"

    #annotate reads
    echo "start annotate reads: ${sra}"
    featureCounts -T ${CORE} -p -t "transcript" -a "../KPN MGH 78578_genome.gtf" -o "reads_annotate/${sra}_annotate.tsv" "bam_merge_paired/${sra}_merge_sort.bam"
    echo "end annotate reads: ${sra}"

    echo "end ${sra}"
done < "sra_list_paired.txt"
echo "END PROCESS PAIRED"

#single reads processing
echo "START PROCESSING SINGLE"
while IFS= read -r sra; do
    echo "start ${sra}"

    #download SRA
    echo "start download SRA: ${sra}"
    fasterq-dump "cache/sra/${sra}.sra" --split-3 -e ${CORE} -O "sra_fastq" 
    rm "cache/sra/${sra}.sra" #clear SRA cache after split to FASTQ
    echo "end download SRA: ${sra}"

    #processFASTQ (single)
    echo "start process FASTQ: ${sra}"
     #trimming
    sickle se -t sanger -q 25 -f "sra_fastq/${sra}.fastq" \
        -o "fastq_trimmed_single/${sra}_trim.fastq"
    #mapping to KPN genome
    bowtie2 -q --fr -x "bt2index/KPN" \
        -U "fastq_trimmed_single/${sra}_trim.fastq" \
        -S "${sra}_single.sam" -p ${CORE}
    #conversion to .bam
    samtools view --threads ${ADDCORE} -bS "${sra}_single.sam" > "${sra}_single.bam"
    #sorting
    samtools sort --threads ${ADDCORE} "${sra}_single.bam" -o "bam_sort_single/${sra}_single_sort.bam"
    #remove temporary unsorted .sam and .bam mapping
    rm "${sra}_single.sam" "${sra}_single.bam"
    #remove temporary fastq
    rm "sra_fastq/${sra}.fastq"
    #remove temporary trimmed FASTQ
    rm -rf "fastq_trimmed_single"
    mkdir -p "fastq_trimmed_single"
    #conversion to tsv
    samtools index --threads ${ADDCORE} "bam_sort_single/${sra}_single_sort.bam" -b "bam_sort_single/${sra}_single_sort.bam.bai"
    samtools idxstats --threads ${ADDCORE} "bam_sort_single/${sra}_single_sort.bam" > "maps_single/${sra}_map.tsv"
    echo "end process FASTQ: ${sra}"

    #annotate reads
    echo "start annotate reads: ${sra}"
    featureCounts -T ${CORE} -t "transcript" -a "../KPN MGH 78578_genome.gtf" -o "reads_annotate/${sra}_annotate.tsv" "bam_sort_single/${sra}_single_sort.bam"
    echo "end annotate reads: ${sra}"

    echo "end ${sra}"
done < "sra_list_single.txt"
echo "END PROCESS SINGLE"
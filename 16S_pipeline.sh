#!/bin/bash

# ========== BƯỚC 1: LỌC TẠP NHIỄM VÀ CHẤT LƯỢNG ========== #
mkdir -p kneaddata_output

for sample in $(ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//'); do
  kneaddata \
    --input1 ${sample}_R1.fastq.gz \
    --input2 ${sample}_R2.fastq.gz \
    --output kneaddata_output/${sample} \
    --reference-db /mnt/d/Do_an_tot_nghiep/NCBIdatabase/bowtie2_index_hg38 \
    --trimmomatic $(dirname $(which trimmomatic)) \
    --threads 4 \
    --remove-intermediate-output \
    --log kneaddata_output/${sample}/kneaddata.log
done

# ========== BƯỚC 2: NỐI CẶP ĐỌC ========== #
mkdir -p vsearch_output

for sample in $(ls kneaddata_output/*/*_kneaddata_paired_1.fastq | sed 's|kneaddata_output/||' | sed 's|/.*_R1_kneaddata_paired_1.fastq||'); do
  echo ">> Merging sample: $sample"

  vsearch \
    --fastq_mergepairs kneaddata_output/${sample}/${sample}_R1_kneaddata_paired_1.fastq \
    --reverse kneaddata_output/${sample}/${sample}_R1_kneaddata_paired_2.fastq \
    --fastqout vsearch_output/${sample}_merged.fastq \
    --fastqout_notmerged_fwd vsearch_output/${sample}_unmerged_R1.fastq \
    --fastqout_notmerged_rev vsearch_output/${sample}_unmerged_R2.fastq \
    --fastq_minovlen 10 \
    --fastq_maxdiffs 5 \
    --relabel ${sample}_
done

# ========== BƯỚC 3: CHUYỂN FASTQ ➜ FASTA & KHỬ TRÙNG LẶP ========== #
mkdir -p fasta_output derep_output

for sample in $(ls vsearch_output/*_merged.fastq | sed 's|vsearch_output/||' | sed 's|_merged.fastq||'); do
  echo ">> Processing sample: $sample"

  cat vsearch_output/${sample}_merged.fastq \
      vsearch_output/${sample}_unmerged_R1.fastq \
      vsearch_output/${sample}_unmerged_R2.fastq \
      > vsearch_output/${sample}_all.fastq

  vsearch \
    --fastq_filter vsearch_output/${sample}_all.fastq \
    --fastaout fasta_output/${sample}.fasta \
    --relabel ${sample}_

  vsearch \
    --derep_fulllength fasta_output/${sample}.fasta \
    --output derep_output/${sample}_derep.fasta \
    --sizeout \
    --minseqlength 100 \
    --relabel ${sample}_
done

# ========== BƯỚC 4: CLUSTER OTU 97% ========== #
mkdir -p otu_rep_seqs otu_tables

ls derep_output/*_derep.fasta | sed 's|derep_output/||' | sed 's|_derep.fasta||' > sample_list.txt

parallel --jobs 8 '
  echo ">> Clustering OTUs for: {1}"
  vsearch \
    --cluster_size derep_output/{1}_derep.fasta \
    --id 0.97 \
    --centroids otu_rep_seqs/{1}_otus.fasta \
    --relabel {1}_OTU_ \
    --threads 1
  vsearch \
    --usearch_global derep_output/{1}_derep.fasta \
    --db otu_rep_seqs/{1}_otus.fasta \
    --id 0.97 \
    --otutabout otu_tables/{1}_otu_table.txt \
    --threads 1 \
    --strand both \
    --sizein \
    --sizeout
' :::: sample_list.txt

# ========== BƯỚC 5: PHÂN LOẠI OTU BẰNG KRAKEN2 ========== #
mkdir -p kraken2_otu_output
cat sample_list.txt > otu_sample_list.txt

cat otu_sample_list.txt | parallel --jobs 8 '
  echo ">> Running Kraken2 for OTU rep seq: {1}"
  kraken2 \
    --db /mnt/d/Do_an_tot_nghiep/NCBIdatabase \
    --output kraken2_otu_output/{1}.out.txt \
    --report kraken2_otu_output/{1}.report.txt \
    --use-mpa-style \
    --threads 1 \
    --memory-mapping \
    otu_rep_seqs/{1}_otus.fasta
  sed -i "s/{1}_//g" otu_rep_seqs/{1}_otus.fasta
  echo ">> Done: {1}"
  echo "----------------------------------------------"
'

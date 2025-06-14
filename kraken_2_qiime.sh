#!/bin/bash

# Read optional taxonomy filter (e.g., "d__Bacteria")
TAX_FILTER=${1:-""}

# Rename space in report file names
sed -i 's/ /_/g' *report.txt
sed -i 's/CLUcontig/OTU/g' *out.txt

# Step 0: Find all *_report.txt files
FILES=( *report.txt )

if [ ${#FILES[@]} -eq 0 ]; then
  echo "❌ No *_report.txt files found in current directory."
  exit 1
fi

echo "🧪 Found ${#FILES[@]} report files."

# Step 1: Collect all unique taxonomy paths
echo "📦 Collecting all unique taxa..."
cut -f1 -d$'\t' "${FILES[@]}" | sort | uniq > all_taxa_raw.txt

# Step 2: Assign Feature IDs (OTU_1, OTU_2, ...) to each unique taxon
if [[ -n "$TAX_FILTER" ]]; then
  echo "🔍 Filtering taxonomy by: '$TAX_FILTER'"
  grep "$TAX_FILTER" all_taxa_raw.txt | awk '{printf "OTU_%d\t%s\n", NR, $0}' > taxa_with_ids.tsv
else
  echo "🔍 No taxonomy filter applied — using all taxa."
  awk '{printf "OTU_%d\t%s\n", NR, $0}' all_taxa_raw.txt > taxa_with_ids.tsv
fi

cut -f1 taxa_with_ids.tsv > feature_ids.txt
cut -f2- taxa_with_ids.tsv > all_taxa.txt

# Step 3: Create taxonomy_qiime.tsv (from k__ to s__, no d__, no leading ;)
echo "📋 Creating taxonomy_qiime.tsv..."
paste feature_ids.txt all_taxa.txt | awk -F'\t' 'BEGIN {
    print "Feature ID\tTaxon\tConfidence"
}
{
    id = $1
    tax = $2
    levels["k__"] = "k__"
    levels["p__"] = "p__"
    levels["c__"] = "c__"
    levels["o__"] = "o__"
    levels["f__"] = "f__"
    levels["g__"] = "g__"
    levels["s__"] = "s__"
    split(tax, ranks, /\|/)
    for (i in ranks) {
        split(ranks[i], kv, "__")
        prefix = kv[1]"__"
        value = (length(kv) > 1 && length(kv[2]) > 0) ? kv[1]"__"kv[2] : prefix
        if (prefix in levels) {
            levels[prefix] = value
        }
    }
    tax_path = levels["k__"]"; "levels["p__"]"; "levels["c__"]"; "levels["o__"]"; "levels["f__"]"; "levels["g__"]"; "levels["s__"]
    print id "\t" tax_path "\t1.0"
}' > taxonomy_qiime.tsv

# Step 4: Prepare abundance matrix without taxonomy
cut -f2- taxa_with_ids.tsv > taxon_key.txt
cp taxon_key.txt taxon_tmp.txt
cut -f1 taxa_with_ids.tsv > feature_tmp.txt

# Step 5: Build abundance counts (only FeatureID + counts)
for file in "${FILES[@]}"; do
    awk -F'\t' '{print $1"\t"$2}' "$file" | sort > "${file}.sorted"
    join -t $'\t' -a1 -e 0 -o auto taxon_tmp.txt "${file}.sorted" > tmp && mv tmp taxon_tmp.txt
done

# Step 6: Create table_qiime.tsv with only abundances
echo "📊 Creating table_qiime.tsv..."
printf "#OTU ID" > table_qiime.tsv
for file in "${FILES[@]}"; do
    sample_name=$(basename "$file" _report.txt)
    fasta_file=$(ls "${sample_name}.fas" 2>/dev/null)
    if [[ -f "$fasta_file" ]]; then
        sample_name=$(basename "$fasta_file" .fas)
    fi
    printf "\t%s" "$sample_name" >> table_qiime.tsv
done
printf "\n" >> table_qiime.tsv

awk 'BEGIN{OFS="\t"} NR==FNR{fid[NR]=$1; next} {
    printf "%s", fid[FNR]
    for(i=2;i<=NF;i++) {
        val = $i
        printf "\t%s", (val=="") ? 0 : val
    }
    printf "\n"
}' feature_tmp.txt taxon_tmp.txt >> table_qiime.tsv

# Step 7.1: generate representative seqs from OTU table
otu_count=$(grep 'OTU_' table_qiime.tsv | wc -l)

for id in $(seq 0.97 0.01 0.99); do
    vsearch --cluster_fast <(cat *.fas) --id $id --centroids tmp.fasta --quiet
    count=$(grep '>' tmp.fasta | wc -l)
    echo "$id identity -> $count clusters"

    if [ "$count" -gt "$otu_count" ]; then
        echo "✅ Found match at $id identity with $count clusters"
        echo "📏 Selecting top $otu_count longest sequences..."

        awk 'BEGIN {RS=">"; ORS=""} NR>1 {
            split($0, lines, "\n");
            header=lines[1];
            seq="";
            for (i=2; i<=length(lines); i++) {
                seq=seq lines[i];
            }
            print length(seq) "\t" seq "\n";
        }' tmp.fasta | sort -nr | head -n "$otu_count" | \
        awk -F'\t' -v OFS="\n" '{printf(">OTU_%d\n%s\n", NR, toupper($2))}' > OTU_Seqs.fas

        echo "✅ Saved to final_representatives.fasta with OTU headers and unwrapped sequences."
        break
    fi
done

# Step 8: Cleanup
#rm -f *.sorted all_taxa.txt all_taxa_raw.txt taxon_tmp.txt taxon_key.txt taxa_with_ids.tsv feature_ids.txt feature_tmp.txt tmp top_seqs.tmp 

# Step 9: Convert to QIIME2 format
biom convert -i table_qiime.tsv -o table_qiime.biom --to-hdf5
qiime tools import --type FeatureTable[Frequency] --input-path table_qiime.biom --output-path table_qiime.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-path taxonomy_qiime.tsv --output-path taxonomy_qiime.qza
#qiime tools import --type 'FeatureData[Sequence]' --input-path OTU_Seqs.fas --output-path OTU_Seqs.qza
# step 8.1: create output for picrust
#magus -i OTU_Seqs.fas -d ./ -o OTU_mg.fas
# remove files
#rm -r decomposition graph log.txt log_debug.txt subalignments tasks
#picrust2_pipeline.py -s OTU_mg.fas -i table_qiime.biom -t epa-ng -m emp_prob --verbose  -o outpicrust2 --stratified --per_sequence_contrib
#cp -r outpicrust2/pathways_out .
#rm -r outpicrust2

# Step 10: Krona
qiime krona collapse-and-plot \
  --i-table table_qiime.qza \
  --i-taxonomy taxonomy_qiime.qza \
  --o-krona-plot AllSamples_krona.qzv
unzip AllSamples_krona.qzv -d krona_files

# Step 11: Bar plot
qiime taxa barplot \
  --i-table table_qiime.qza \
  --i-taxonomy taxonomy_qiime.qza \
  --o-visualization allsamples_Barplot.qzv
unzip allsamples_Barplot.qzv -d barplot_files

# Step 12: Cleanup again
rm taxonomy_qiime.tsv table_qiime.tsv tmp.fasta

# Step 13: Transpose CSVs
for file in barplot_files/*/data/level*.csv; do
  [ -f "$file" ] || continue
  filename=$(basename "$file")
  dirname=$(dirname "$file")
  output="${dirname}/${filename%.csv}_trans.csv"

  echo "🔄 Transposing $file -> $output"

  awk '
  {
    for (i=1; i<=NF; i++) {
      a[NR,i] = $i
    }
  }
  NF > p { p = NF }
  END {
    for (j=1; j<=p; j++) {
      str = a[1,j]
      for (i=2; i<=NR; i++) {
        str = str "," a[i,j]
      }
      print str
    }
  }' FS=',' OFS=',' "$file" > "$output"
  mv "$output" .
done

rm -r barplot_files

echo "✅ All outputs complete. Ready for R or QIIME2 visualization!"
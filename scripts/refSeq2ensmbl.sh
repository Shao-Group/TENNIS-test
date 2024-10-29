i=$1

cvbio UpdateContigNames \
    -i $i  \
    -o $i.ucsc_name.gtf  \
    -m ../../data/ChromosomeMappings/GRCh38_RefSeq2UCSC.txt --comment-chars '#' \
    --columns 0 \
    --skip-missing true


cvbio UpdateContigNames \
    -i $i.ucsc_name.gtf \
    -o $i.ensembl_name.gtf \
    -m ../../data/ChromosomeMappings/GRCh38_UCSC2ensembl.txt \
    --comment-chars '#' \
    --columns 0 \
    --skip-missing true

cat $i.ensembl_name.gtf | awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/ {print}' > $i.MainChrName.ensembl_name.gtf

rm $i.ucsc_name.gtf 
rm $i.ensembl_name.gtf

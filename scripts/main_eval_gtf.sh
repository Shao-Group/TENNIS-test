#!/bin/bash

input_gtf=$1
tennis_gtf=$2
truth_gtf=$3
translate_txt=$4
eval_output_prefix=$5

psi_file=""
chr_translate_file=""
tennis_test_dir=""

# 5 randoms
for i in $(seq 1 5)
do
    $seed="--seed $i"
    if [[ $i -eq 1 ]]; then
        $seed=""
    fi 
    tennis test -f Random1 -x 100 \
        -p 0.0 \
        -o Rand1_run"$i" \
        $seed \
        --xi_gtf_file $tennis_gtf \
        $input_gtf
    tennis test -f RandomX -x 100 \
        -p 0.0 \
        -o RandX_run"$i" \
        $seed \
        --xi_gtf_file $tennis_gtf \
        $input_gtf
done

# PSI1 and PSIX
for j in "1" "X"
do
    tennis test -f PSI"$j" \
        -x 100 -p 0.0 \
        -o PSI"$j"_run"$i" \
        --psi_file $psi_file \
        --xi_gtf_file $tennis_gtf \
        --chr_translate_file $chr_translate_file \
        $input_gtf
done

# gffcompare
for i in $(ls  Rand*_run*.pred.gtf PSI*_run*.pred.gtf)
do
    gffcomapre -r $truth_gtf -o "$eval_output_prefix"_"${i%.pred.gtf}"_eval $i
done


# score TENNIS by PSI
tennis test -f ScorePSI \
    -x 100000000 -p 0.0 \
	-o TENNIS_SAT_psi_score \
	--psi_file $psi_file \
	--predicted_gtf $tennis_gtf \
	--chr_translate_file $chr_translate_file \
	$input_gtf

# TENNIS ROC by PctIn_x_PSI
python $tennis_test_dir/scripts/precision_recall_by_pctIn.py \
    TENNIS_SAT_psi_score.scored.gtf \
    $tennis_eval_tmap  \
    pctIn_PSI_score  \
    x_score.txt

# fig 
# python /Users/xzang/Work/research_code/tennis-test-public/scripts/precision_recall_fig.py \
#   x_score.txt \
#   156.6 24.44 15.17   2.338 \
#   186.8 20    15.0386 1.6199 \
#   232 36.2 \
#   281 30.1 \
#   dm6.lr.pdf

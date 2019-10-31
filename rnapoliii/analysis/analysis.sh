#!/bin/sh

ls imp_tutorial_pol3_xl_cryoem_premodeling/modeling/

ls imp_tutorial_pol3_xl_cryoem_premodeling/modeling/A_output_3/rmfs/

mkdir -p modeling/run1
ln -sf ../../imp_tutorial_pol3_xl_cryoem_premodeling/modeling/A_output_3 \
       modeling/run1/output

python scripts/select_good_scoring_models.py \
    -rd modeling -rp run \
    -sl CrossLinkingMassSpectrometryRestraint_Distance_ \
    -pl ConnectivityRestraint_ABC10alpha \
        CrossLinkingMassSpectrometryRestraint_Data_Score_XL \
        ExcludedVolumeSphere_None \
        GaussianEMRestraint_Total \
        Total_Score \
    -alt 0.9 -aut 1.0 -mlt 0.0 -mut 30.0

head filter/model_ids_scores.txt

python scripts/plot_score.py \
    filter/model_ids_scores.txt \
    GaussianEMRestraint_Total

python scripts/plot_score.py \
    filter/model_ids_scores.txt \
    ConnectivityRestraint_ABC10alpha

python scripts/plot_score.py \
    filter/model_ids_scores.txt \
    ExcludedVolumeSphere_None

python scripts/select_good_scoring_models.py \
    -rd modeling -rp run \
    -sl CrossLinkingMassSpectrometryRestraint_Distance_ \
        GaussianEMRestraint_Total \
        ExcludedVolumeSphere_None \
    -pl ConnectivityRestraint_ABC10alpha \
        CrossLinkingMassSpectrometryRestraint_Data_Score_XL \
        Total_Score \
    -alt 0.9 -50 -50 -aut 1.0 5450 90 \
    -mlt 0.0 0.0 0.0 -mut 30.0 0.0 0.0
    
wc -l filter/model_ids_scores.txt

cat <<END > density_ranges.txt
density_custom_ranges={
    'C82':['C82'],
    'C34':['C34'],
    'C31':['C31'] }
END

python scripts/select_good_scoring_models.py \
    -rd modeling -rp run \
    -sl CrossLinkingMassSpectrometryRestraint_Distance_ \
        GaussianEMRestraint_Total \
        ExcludedVolumeSphere_None \
    -pl ConnectivityRestraint_ABC10alpha \
        CrossLinkingMassSpectrometryRestraint_Data_Score_XL \
        Total_Score \
    -alt 0.9 -50 -50 -aut 1.0 5450 90 \
    -mlt 0.0 0.0 0.0 -mut 30.0 0.0 0.0 -e

python scripts/Master_Sampling_Exhaustiveness_Analysis.py \
       -n rnapoliii -p good_scoring_models/ \
       -d density_ranges.txt -m cpu_omp \
       -c 8 -a -g 0.1 -gp

mkdir -p modeling/run1
ln -sf ../../../../../imp_tutorial_pol3_xl_cryoem_premodeling/modeling/A_output_3 \
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

python scripts/plot_score.py \
    filter/model_ids_scores.txt \
    GaussianEMRestraint_Total

python scripts/select_good_scoring_models.py \
    -rd modeling -rp run \
    -sl CrossLinkingMassSpectrometryRestraint_Distance_ \
        GaussianEMRestraint_Total \
    -pl ConnectivityRestraint_ABC10alpha \
        CrossLinkingMassSpectrometryRestraint_Data_Score_XL \
        ExcludedVolumeSphere_None \
        Total_Score \
    -alt 0.9 -50 -aut 1.0 5440 -mlt 0.0 0.0 -mut 30.0 0.0 -e

python scripts/Master_Sampling_Exhaustiveness_Analysis.py \
       -n rnapoliii -p good_scoring_models/ \
       -d density_ranges.txt -m cpu_omp \
       -c 8 -a -g 0.1 -gp

method=maml
dims=(32 128 256)
for dim in "${dims[@]}";do
    python run.py experiment_subset.feature_class=Promoter \
    exp.name="atacseq_adult_promoter_${method}_dim_${dim}_${dim}" \
    backbone.layer_dim=[${dim},${dim}] 
done

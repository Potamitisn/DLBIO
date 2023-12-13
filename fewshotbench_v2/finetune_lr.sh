method=maml

if [ "$method" == "maml" ]; then
    #lrSamples=(0.02 0.05 0.1)
    lrSamples=(0.07 0.2)
    for lr in "${lrSamples[@]}";do
        python run.py experiment_subset.feature_class=Promoter exp.name="atacseq_adult_promoter_${method}_${lr}lr" method.maml_inner_lr=${lr} 
    done
else 
    lrSamples=(0.002 0.005 0.01)
    for lr in "${lrSamples[@]}"; do
        python run.py experiment_subset.feature_class=Promoter exp.name="atacseq_adult_promoter_${method}_${lr}lr" lr=${lr} 
    done
fi
nWaySamples=(1 10 20 35 50)
method=maml

for nWay in "${nWaySamples[@]}"
do
    python run.py experiment_subset.feature_class=Promoter exp.name="atacseq_adult_promoter_${method}_${nWay}way" method=${method} n_way=${nWay} 
done
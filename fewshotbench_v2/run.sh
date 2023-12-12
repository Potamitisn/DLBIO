method=baseline
param=n_shot # [n_way, n_shot]
nWaySamples=(1 35 50)
nShotSamples=(50 100)


if [ "$param" == "n_way" ]; then
    echo "Studying ${method} for different n_way"
    for nWay in "${nWaySamples[@]}"
    do
        python run.py experiment_subset.feature_class=Promoter exp.name="atacseq_adult_promoter_${method}_${nWay}way" method=${method} n_way=${nWay} 
    done
else 
    echo "Studying ${method} for different n_shot"
    for nShot in "${nShotSamples[@]}"
    do
        python run.py experiment_subset.feature_class=Promoter exp.name="atacseq_adult_promoter_${method}_${nShot}shot" method=${method} n_shot=${nShot} 
    done
fi
# DLBIO

This repository includes our implementations for the final project of the course Deep learning in iomedicine , EPFL.

Group 8 : [Hanrong Hu](https://github.com/hanrong498), [Nearchos Potamitis](https://github.com/Potamitisn), [Karim Abi Said](https://github.com/KarimAbiSaid).

## The dataset

In this report, we evaluate several few-shot learning models
on a newly compiled ATAC-seq dataset of human cells
([Zhang et al., 2021](https://www.sciencedirect.com/science/article/pii/S0092867421012794)), encompassing fetal and adult samples.
The dataset includes 30 different tissues from adult samples
with 111 cell types, and 16 tissues from fetal samples with
the same number of cell types. Overall, over 1.3 million
single nuclei were sequenced, identifying approximately
1.2 million candidate CREs. The candidate CREs were
classified in to promoter (-200 to 200 bp to transciption
starting site), promoter proximal and distal.

### Acquiring the dataset
There are 2 options to acquire the necessary data:

### EDA
🔴 Say a few words about our eda notebooks after we have their final versions.🔴 

#### 1. Manually
Download the following files (by clicking on the links) and move them to `DLBIO/fewshotbench_v2/data/atacseq/`.
- [Cell_metadata.tsv.gz](https://storage.googleapis.com/atacseq-np/data/atacseq/Cell_metadata.tsv.gz)
- [cCRE_hg38.tsv.gz](https://storage.googleapis.com/atacseq-np/data/atacseq/cCRE_hg38.tsv.gz)
- [matrix.h5ad](https://storage.googleapis.com/atacseq-np/data/atacseq/matrix.h5ad)

#### 2. Google Cloud
If you're already using a google cloud VM instance then the preferrable way would be to directly download the data from our Cloud Storage Bucket to your instance. To do that SSH into your instance and run the following commands in the terminal:
```
gcloud storage cp gs://atacseq-np/data/atacseq/Cell_metadata.tsv.gz DLBIO/fewshotbench_v2/data/atacseq/
gcloud storage cp gs://atacseq-np/data/atacseq/cCRE_hg38.tsv.gz DLBIO/fewshotbench_v2/data/atacseq/
gcloud storage cp gs://atacseq-np/data/atacseq/matrix.h5ad DLBIO/fewshotbench_v2/data/atacseq/
```

## Running the benchmarks
Our work has been fully integrated with the given "Few Shot Benchmark for Biomedical Datasets" (`fewshotbench_v2`). If you're not familiar with that take a look at its [readme file](https://github.com/Potamitisn/DLBIO/blob/main/fewshotbench_v2/README.md) and make sure to create and activate the conda environment mentioned. Otherwise, the `dataset` of the `main.yaml` file has already been set to `atacseq` and as soon as the data is downloaded the code is ready to be tested! For example you can test the MAML method simply by running the following command in the terminal:
```
python fewshotbench_v2/run.py exp.name=test_grou_8 method=maml dataset=atacseq
```

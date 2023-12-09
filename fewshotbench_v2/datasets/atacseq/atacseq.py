from abc import ABC

import numpy as np
from torch.utils.data import DataLoader

from datasets.dataset import *
from datasets.atacseq.utils import AtacData



class AtSDataset(FewShotDataset, ABC):
    _dataset_name = 'atacseq'
    #_dataset_url = "http://catlas.org/catlas_downloads/humantissues/"
    _dataset_url = "https://drive.google.com/uc?export=download&id=1MCCpTq1Xi6uQ-oHgVhOsPAZi9y5zdeZH"

    def load_atac_seq(self, mode='train', min_samples=20):
        dclass = AtacData()
        adata = dclass.adata
        
        if dclass.life_stage == "Adult":
            test_tissues = ["pancreas", "adrenal_gland", "thyroid", "islet", "ovary"]
            val_tissues = ["heart_atrial_appendage","heart_la","heart_lv","heart_ra","heart_rv"]
        elif dclass.life_stage == "Fetal":
            test_tissues = ["cerebrum", "cerebellum", "standard"]
            val_tissues = ["eye", "intestine", "thymus"]
        
        train_tissues = []
        for tissue_type in adata.obs["tissue"].unique():
            if tissue_type not in val_tissues+test_tissues :
                train_tissues.append(tissue_type)

        split = {'train': train_tissues,
                 'val': val_tissues,
                 'test': test_tissues}
        
        ## get an indicator of what data/mode is being processed
        print("===== Current data mode is:", mode, "=======")
        tissues = split[mode]
        # subset data based on target tissues
        print("Subset data based on type of tissues")
        adata = adata[adata.obs['tissue'].isin(tissues)]

        print("Subset data based on number of samples")
        filtered_index = adata.obs.groupby(["label"]) \
            .filter(lambda group: len(group) >= min_samples) \
            .reset_index()['barcodes']
        adata = adata[filtered_index]

        # convert gene to torch tensor x
        print("Convert gene to torch tensor")
        samples = adata.to_df().to_numpy(dtype=np.float32)
        # convert label to torch tensor y
        print("Convert label to torch tensor")
        targets = adata.obs['label'].to_numpy(dtype=np.int32)
        
        print("====== Loading", mode, "data done :)")
        return samples, targets


class AtSSimpleDataset(AtSDataset):
    def __init__(self, batch_size, root='./data/', mode='train', min_samples=20):
        self.initialize_data_dir(root, download_flag=False)
        self.samples, self.targets = self.load_atac_seq(mode, min_samples)
        self.batch_size = batch_size
        super().__init__()

    def __getitem__(self, i):
        return self.samples[i], self.targets[i]

    def __len__(self):
        return self.samples.shape[0]

    @property
    def dim(self):
        return self.samples.shape[1]

    def get_data_loader(self) -> DataLoader:
        data_loader_params = dict(batch_size=self.batch_size, shuffle=True, num_workers=4, pin_memory=True)
        data_loader = torch.utils.data.DataLoader(self, **data_loader_params)

        return data_loader


class AtSSetDataset(AtSDataset):

    def __init__(self, n_way, n_support, n_query, n_episode=100, root='./data', mode='train'):
        self.initialize_data_dir(root, download_flag=False)

        self.n_way = n_way
        self.n_episode = n_episode
        min_samples = n_support + n_query

        samples_all, targets_all = self.load_atac_seq(mode, min_samples)
        self.categories = np.unique(targets_all)  # Unique cell labels
        self.x_dim = samples_all.shape[1]

        self.sub_dataloader = []

        sub_data_loader_params = dict(batch_size=min_samples,
                                      shuffle=True,
                                      num_workers=0,  # use main thread only or may receive multiple batches
                                      pin_memory=False)
        for cl in self.categories:
            samples = samples_all[targets_all == cl, ...]
            sub_dataset = FewShotSubDataset(samples, cl)
            self.sub_dataloader.append(torch.utils.data.DataLoader(sub_dataset, **sub_data_loader_params))

        super().__init__()

    def __getitem__(self, i):
        return next(iter(self.sub_dataloader[i]))

    def __len__(self):
        return len(self.categories)

    @property
    def dim(self):
        return self.x_dim

    def get_data_loader(self) -> DataLoader:
        sampler = EpisodicBatchSampler(len(self), self.n_way, self.n_episode)
        data_loader_params = dict(batch_sampler=sampler, num_workers=4, pin_memory=True)
        data_loader = torch.utils.data.DataLoader(self, **data_loader_params)
        return data_loader
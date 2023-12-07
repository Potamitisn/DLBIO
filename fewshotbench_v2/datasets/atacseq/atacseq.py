from abc import ABC

import numpy as np
from torch.utils.data import DataLoader

from datasets.dataset import *
from datasets.atacseq.utils import AtacData



class ASDataset(FewShotDataset, ABC):
    _dataset_name = 'atacseq'

    def load_atac_seq(self, mode='train', min_samples=20):
        adata = AtacData().adata
        train_tissues = []
        test_tissues = ["pancreas", "adrenal_gland", "thyroid", "islet", "ovary"]
        val_tissues = ["heart_atrial_appendage","heart_la","heart_lv","heart_ra","heart_rv"]

        for tissue_type in adata.obs["tissue"].unique():
            if tissue_type not in val_tissues+test_tissues :
                train_tissues.append(tissue_type)

        split = {'train': train_tissues,
                 'val': val_tissues,
                 'test': test_tissues}
        
        tissues = split[mode]
        # subset data based on target tissues
        adata = adata[adata.obs['tissue'].isin(tissues)]

        filtered_index = adata.obs.groupby(["label"]) \
            .filter(lambda group: len(group) >= min_samples) \
            .reset_index()['index']
        adata = adata[filtered_index]

        # convert gene to torch tensor x
        samples = adata.to_df().to_numpy(dtype=np.float32)
        # convert label to torch tensor y
        targets = adata.obs['label'].cat.codes.to_numpy(dtype=np.int32)
        # go2gene = get_go2gene(adata=adata, GO_min_genes=32, GO_max_genes=None, GO_min_level=6, GO_max_level=1)
        # go_mask = create_go_mask(adata, go2gene)
        return samples, targets


class TMSimpleDataset(ASDataset):
    def __init__(self, batch_size, root='./data/', mode='train', min_samples=20):
        self.initialize_data_dir(root, download_flag=True)
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


class TMSetDataset(ASDataset):

    def __init__(self, n_way, n_support, n_query, n_episode=100, root='./data', mode='train'):
        self.initialize_data_dir(root, download_flag=True)

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
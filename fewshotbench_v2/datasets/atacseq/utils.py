# util function for human atac-seq data pre-processing 
import os, requests, sys
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')


class AtacData():
    def __init__(self, experiment_subset, src_file = "data/atacseq/matrix.h5ad", ) -> None:

        """
        Args:
            src_file: h5ad file 
            life_stage: choose between 'Adult' or 'Fetal'. default: Adult.
            pre_processing: True or False; if filter unimportant cis-regulatory elements and low-quality cells. default: True
            feature_class: Select between None, "Promoter", "Promoter Proximal" and "Distal". default: None.
                Promoter (-200 to +200 of TSS), Promoter Proximal (less) or Distal
            subset_fraction: float (0 - 1), fraction to subset the data. default: None
        """
        life_stage = experiment_subset.life_stage
        feature_class = experiment_subset.feature_class
        subset_fraction = experiment_subset.subset_fraction
        pre_processing = experiment_subset.pre_processing

        ##### Check input args #####
        # Check for valid life_stage
        if life_stage not in ["Adult", "Fetal"]:
            raise ValueError("Life stage must be 'Adult' or 'Fetal'")
        # Check for valid feature_class
        valid_feature_classes = [None, "Promoter", "Promoter Proximal", "Distal"]
        if feature_class not in valid_feature_classes:
            raise ValueError("Feature class must be None, 'Promoter', 'Promoter Proximal', or 'Distal'")
        # Check for valid subset_fraction
        if subset_fraction is not None and not (0 < subset_fraction <= 1):
            raise ValueError("Subset fraction must be a value between 0 and 1, or None")

        ##### Data path #####
        # Initialize the base path
        pp_out_path = "data/atacseq/processed_{}".format(life_stage.lower())

        # Conditionally add feature_class and subset_fraction
        if feature_class is not None:
            pp_out_path += "_{}".format(feature_class.replace(" ", "_").lower())
        if subset_fraction is not None:
            pp_out_path += "_{:.0f}fraction".format(subset_fraction*100)
        # Append the file extension
        pp_out_path += "_matrix.h5ad"
        print("Processed data will be saved to {}".format(pp_out_path))

        #assert os.path.exists(pp_out_path)

        ##### Check if processed data exits, or process it #####
        if os.path.exists(pp_out_path) and pre_processing:
            print("Loading the pre-processed data from memory")
            self.adata = sc.read_h5ad(pp_out_path)
        else:
            print("Preparing data...")
            self.set_up_file(src_file)
            print("Loading the data from memory")
            self.adata = sc.read_h5ad(src_file, backed="r") # data do not load to memory 

            print("Adding metadata information")
            self.adata = self.add_cell_metadata(self.adata)
            self.adata = self.add_CRE_metadata(self.adata)

            # Select Adult or Fetal cells 
            print("Selecting {} cells".format(life_stage))
            if life_stage is not None:
                self.adata = self.adata[self.adata.obs["Life stage"] == life_stage]
                if life_stage == "Adult":
                    self.adata = self.adata.to_memory()[:,self.adata.var["Present in adult tissues"] == "yes"]
                elif life_stage == "Fetal":
                    self.adata = self.adata.to_memory()[:,self.adata.var["Present in fetal tissues"] == "yes"]

            # Select CRE (feature) class
            if feature_class is not None:
                print("Selecting feature class {}".format(feature_class))
                self.adata = self.adata.to_memory()[:,self.adata.var["Class"] == feature_class]
        
            # Pre-processing 
            if pre_processing:
                print("Pre-processing the data")
                self.adata = self.pre_process(self.adata)

            # Save it
            self.adata.write_h5ad(pp_out_path) 

        self.add_labels()

        if subset_fraction:
            print("Subsetting the data to {:.1f}%".format(subset_fraction*100))
            subset_index = self.adata.obs.groupby("cell type", group_keys=False).apply(lambda x : self.sample_fraction(x, subset_fraction)).index
            self.add_labels()
            self.adata = self.adata[subset_index]
            # Save it 
            self.adata.write_h5ad(pp_out_path) 
        print("Dataclass is ready!")
        
    def add_cell_metadata(self, adata, src_file = "data/atacseq/Cell_metadata.tsv.gz"):
        """
        Args:
            src_file: Cell_metadata in tsv.gz format
        """
        cell_metadata = pd.read_csv(src_file, sep='\t', compression='gzip')
        # Rename cellID to barcodes, to be consistant with adata.obs
        cell_metadata.rename(columns={"cellID": "barcodes"}, inplace=True)

        #### Get tissue names: need some mannual work
        cell_metadata['tissue_full'] = cell_metadata['tissue']
        cell_metadata['tissue'] = cell_metadata['tissue'].str.rsplit('_', n=1).str[0]
        # Remove '_sample'
        cell_metadata['tissue'] = cell_metadata['tissue'].str.replace('_sample', '', regex=False)
        # Remove '_CARE' followed by any digit
        cell_metadata['tissue'] = cell_metadata['tissue'].str.replace('_CARE.*', '', regex=True)
        # Remove 'Map'
        cell_metadata['tissue'] = cell_metadata['tissue'].str.replace('Map', '', regex=False)
        # Change to lower 
        cell_metadata['tissue'] = cell_metadata['tissue'].str.lower()

        #### Add cell metadata
        adata.obs = adata.obs.merge(cell_metadata, on='barcodes', how='left')
        # Adata prefers the index to be strings, not [1,2,3...] so we use barcodes as the index
        adata.obs.index = adata.obs.barcodes.astype(str)
        adata.obs = adata.obs.drop(columns=['barcodes'])
        # Change the dtype of columns 
        for col in adata.obs.columns:
            if adata.obs[col].dtype != 'float64':
                adata.obs[col] = adata.obs[col].astype('category')

        return adata
    
    def add_CRE_metadata(self, adata, src_file = "data/atacseq/cCRE_hg38.tsv.gz"):
        """
        Args:
            src_file: cCRE_hg38 file in tsv.gz format
        """
        cCRE_hg38 = pd.read_csv(src_file, sep='\t', compression='gzip')
        # Add a 'Feature_ID' column so that we can add it to the adata.var
        cCRE_hg38['Feature_ID'] = cCRE_hg38['#Chromosome'].astype(str) + ':' + cCRE_hg38['hg38_Start'].astype(str) + '-' + cCRE_hg38['hg38_End'].astype(str)
        
        #### Add CRE information 
        adata.var = adata.var.merge(cCRE_hg38, on='Feature_ID', how='left')
        # Adata prefers the index to be strings
        adata.var.index = adata.var.Feature_ID.astype(str)
        adata.var = adata.var.drop(columns = "Feature_ID")
        
        return adata
    
    def pre_process(self, bdata, out_path=None):
        # Step1: filtering based on n_count and n_cCRE
        sc.pp.filter_cells(bdata, min_counts=500)
        sc.pp.filter_cells(bdata, min_genes=500)
        bdata.obs.rename(columns = {"n_genes" : "n_cCRE"}, inplace = True)

        # Step2: Normalize to library size 
        sc.pp.normalize_total(bdata, target_sum = 5e3)

        # Step3: Select hight variable CREs
        sc.pp.log1p(bdata)
        sc.pp.highly_variable_genes(bdata)
        bdata.raw = bdata
        bdata = bdata[:,bdata.var.highly_variable]
        if out_path:
            bdata.write_h5ad(out_path)

        return bdata
    
    def add_labels(self):
        """
        Add labels to the adata.obs
        """
        # Add labels to the adata.obs
        self.adata.obs["label"] = self.adata.obs['cell type'].astype('category').cat.codes

    def set_up_file(self, src_file):
        server = "http://catlas.org/catlas_downloads/humantissues/"
        link = server + "/".join(src_file.split("/")[1:])
        if not os.path.exists(src_file):
            os.makedirs("/".join(src_file.split("/")[:-1]), exist_ok=True)
            with open(src_file, "wb") as f:
                print("Downloading %s" % src_file)
                response = requests.get(link, stream=True)
                total_length = response.headers.get('content-length')

                if total_length is None: # no content length header
                    f.write(response.content)
                else:
                    dl = 0
                    total_length = int(total_length)
                    for data in response.iter_content(chunk_size=4096):
                        dl += len(data)
                        f.write(data)
                        done = int(50 * dl / total_length)
                        sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)) )    
                        sys.stdout.flush()
        else:
            print("File exists: %s (no need to download)" % src_file)

    def sample_fraction(self, group, fraction):
        return group.sample(frac=fraction)



        

        

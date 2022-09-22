import pandas as pd
import mygene
import subprocess
import itertools
from utils import Utils
import glob

class GeneExp():
    def __init__(self,set_name = "all"):
        self.set_name = set_name
        self.set_raw_data = None

        self.id_type_to_sets = {}
        self.local_ens_to_syb_map = {}

    #region Public

    def return_gene_enrch_sets(self,p_val_trh = 0.05,id_type = "symbol"):
        if id_type in self.id_type_to_sets.keys() : return self.id_type_to_sets[id_type]

        enrch_genes, non_enrch_genes = self._return_gene_enrch_sets(p_val_trh)

        if id_type == "ensembl.gene" :
            self.id_type_to_sets[id_type] = (enrch_genes, non_enrch_genes)
        if id_type == "symbol" :
           (_enrch_genes, _non_enrch_genes),self.local_ens_to_syb_map = Utils.convert_ensembl_to_symbols(
               [enrch_genes,non_enrch_genes],self.local_ens_to_syb_map)
           self.id_type_to_sets[id_type] = (_enrch_genes, _non_enrch_genes)

        return self.id_type_to_sets[id_type]

    #endregion

    #region Private

    def _return_gene_enrch_sets(self,p_val_trh = 0.05):
        full_table = self._loda_genes_from_dataset()
        if type(full_table) is list :
            enrch_genes,non_enrch_genes = [],[]
            for table in full_table :
                enrch_genes.append(table[table["pval"] < p_val_trh]["id"].to_list())
                non_enrch_genes.append(table[table["pval"] >= p_val_trh]["id"].to_list())

        else :
            enrch_genes = full_table[full_table["pval"] < p_val_trh]["id"].to_list()
            non_enrch_genes = full_table[full_table["pval"] >= p_val_trh]["id"].to_list()

        return enrch_genes, non_enrch_genes

    def _loda_genes_from_dataset(self):
        if self.set_raw_data is not None :
            return self.set_raw_data

        if self.set_name != "all" :
            full_table = pd.read_table(f"../../data/From_DOMINO/GE/{self.set_name}.tsv", header=0)
        else :
            all_tsvs_files = glob.glob(f"../../data/From_DOMINO/GE/*.tsv")
            full_table = []
            for tsf in all_tsvs_files :
                table = pd.read_table(tsf, header=0)
                full_table.append(table)

        self.set_raw_data = full_table
        return full_table

    #endregione
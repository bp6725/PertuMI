import pandas as pd
import gseapy as gp




class Pathways():
    def __init__(self,gmt_path = "/home/benny/Repos/knowledge_base_clustering/data/Pathways/c2.cp.kegg.v7.5.1.symbols.gmt"):
        self.gmt_file = gp.parser.gsea_gmt_parser(gmt_path)

    #region Public

    @staticmethod
    def get_pathways_lib_names(condition_selector=None):
        """Returns the names of the KEGG pathways associated to the selected condition.
        Parameters
        ----------
        condition_selector : ConditionSelector
            Specifies for which condition the associated pathways should be loaded.
        Returns
        -------
        pathways : list of str
            Names of phenotype-related KEGG pathways.
        """
        if condition_selector == "ALS":
            return ['hsa05014']
        elif condition_selector == "LC":
            return ['hsa05223']
        elif condition_selector == "UC":
            return ['hsa04060', 'hsa04630', 'hsa05321']
        elif condition_selector == "HD":
            return ['hsa05016']
        elif condition_selector == "CD":
            return ['hsa04621', 'hsa04060', 'hsa04630', 'hsa05321', 'hsa04140']

        return pd.read_csv("../../data/all_human_path_kegg", header=None).iloc[0].to_list()

    def get_gmt_file(self,path = None, gene_id_types = "symbols"):
        if gene_id_types == "symbols" :
            gmt =  self.gmt_file if path is None else gp.parser.gsea_gmt_parser(path)
        if gene_id_types == "ensembl.gene" :
            raise NotImplemented()

        return gmt

    #endregion
import numpy as np
import random
from DAL.pathways import Pathways
import itertools

class BuilderSyntheticData():
    def __init__(self):
        pass

    # region Public

    @staticmethod
    def generate_enriched_ge_for_pathways(ppi_gene_set, all_pathways_dict,pathways=3, random_samples_perc = 0.15):
        '''
        if pathways is int pick random pathways, if list pick those pathways.
        :param pathways:
        :return:
        '''
        if type(pathways) is int:
            picked_pathways, picked_pathways_names = BuilderSyntheticData._pick_n_random_pathways(pathways,all_pathways_dict)
            n_genes_per_pathway = list(map(lambda x: random.randint(int(len(x) / 2), int(len(x))), picked_pathways))
            picked_genes_per_pathway = [random.sample(path, n_genes) for path, n_genes in
                                        zip(picked_pathways, n_genes_per_pathway)]

            random_samples = random.sample(ppi_gene_set,int(len(ppi_gene_set)*random_samples_perc))
            picked_pathways_names.append("random samples")
            picked_genes_per_pathway.append(random_samples)
            return picked_pathways_names,picked_genes_per_pathway



    # endregion

    @staticmethod
    def _pick_n_random_pathways(n_paths,all_pathways_dict):
        picked_paths_items = random.sample(all_pathways_dict.items(),n_paths)
        return list(map(lambda x:x[1],picked_paths_items)), list(map(lambda x:x[0],picked_paths_items))
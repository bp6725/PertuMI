from unittest import TestCase
import numpy as np
import random

from DAL.ppis import PPIS
from DAL.gene_exp import GeneExp
from DAL.pathways import Pathways
from models.community_detection import CommunityDetection
from models.pcst import PCST
from models.gesa_infr import GSEaInfr


class Testgesa_infr(TestCase):
    def test_compute_neg_log_gsea_p_value(self):
        def __pick_genes_from_random_pathway(all_pathways_dict, perc_to_pick):
            random_phws = np.random.choice(list(all_pathways_dict.keys()))
            genes = all_pathways_dict[random_phws]

            n_to_pick = int(perc_to_pick * len(genes))
            return random.sample(genes,n_to_pick), random_phws

        # genes to test
        phws = Pathways()
        all_pathways_dict = phws.get_gmt_file()

        genes_enr, phws_enr = __pick_genes_from_random_pathway(all_pathways_dict,0.8)
        genes_noise, phws_noise = __pick_genes_from_random_pathway(all_pathways_dict, 0.3)
        genes_noise1, phws_noise1 = __pick_genes_from_random_pathway(all_pathways_dict, 0.3)

        all_genes = list(set(genes_enr+genes_noise+genes_noise1))

        # gmt file
        gmt_path = phws.gmt_file

        #GSEA
        gi = GSEaInfr()

        single_pathway = {phws_enr:all_pathways_dict[phws_enr]}
        res = gi.compute_pval_local_single_pathway(all_genes,single_pathway, background = 1000)
        print(res)

        self.fail()

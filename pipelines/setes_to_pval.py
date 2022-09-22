from tqdm import tqdm

from models.pcst import PCST
from DAL.ppis import PPIS
from utils import Utils
from DAL.gene_exp import GeneExp
from DAL.pathways import Pathways
from models.community_detection import CommunityDetection
from models.gesa_infr import GSEaInfr
from pipelines.pipelines_codebase import PipelinesCodebase

class SetsToPval():
    '''
    We get genes set, currently relevant pathways and return pval per pathways as before
    '''
    def __init__(self):
        pass

    @staticmethod
    def run_pipeline(enriched_genes_set,non_enriched_genes_set,
                     relevant_pathways_dict, ppi_nx, all_clusters_sets, level_mapping ,background,n_permut = 5,
                     pct = PCST(), gi = GSEaInfr(), cd = CommunityDetection(),with_original = True):
        final_results = {}

        if with_original :
            ge_ami = PipelinesCodebase.run_ami_detection(enriched_genes_set, non_enriched_genes_set, ppi_nx,
                                                         relevant_pathways_dict, all_clusters_sets, level_mapping,background,
                                                         cd, pct, gi)

            final_results["original"] = ge_ami

        for permu_index in range(n_permut) :
            _peru_enrch_genes, _peru_non_enrch_genes = PipelinesCodebase.build_permuted_sets(enriched_genes_set,
                                                                                             non_enriched_genes_set)

            _peru_ge_ami = PipelinesCodebase.run_ami_detection(_peru_enrch_genes, _peru_non_enrch_genes, ppi_nx,
                                                               relevant_pathways_dict, all_clusters_sets,
                                                               level_mapping,background, cd, pct, gi)

            final_results[permu_index] = _peru_ge_ami


        return final_results
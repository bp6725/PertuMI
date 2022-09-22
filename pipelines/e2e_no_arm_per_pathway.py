from tqdm import tqdm

from models.pcst import PCST
from DAL.ppis import PPIS
from DAL.gene_exp import GeneExp
from DAL.pathways import Pathways
from models.community_detection import CommunityDetection
from models.gesa_infr import GSEaInfr
from pipelines_codebase import PipelinesCodebase

class E2eNoArmPerPpath():
    def __init__(self):
        pass

    @staticmethod
    def run_pipeline(ppi_name = "huri" ,ge_name = "cbx",n_permut = 500,min_cluster_size=15, ge_enrich_trh = 0.05):

        # Load datasets
        ## ppi network
        ppi_nx = PPIS.read_ppi_network(network=ppi_name,min_cluster_size=min_cluster_size,node_id_type="symbols")
        ppi_genes_set = set(ppi_nx.nodes)

        ## enriched genes
        ge = GeneExp(ge_name)
        enrch_genes,non_enrch_genes = ge.return_gene_enrch_sets(p_val_trh = ge_enrich_trh)

        phws = Pathways()
        _all_pathways_dict = phws.get_gmt_file()
        # We filter the pathways by removing genes which are not in the ppi
        all_pathways_dict = {path_name:list(filter(lambda x:x in ppi_genes_set,path_genes))
                             for path_name,path_genes in _all_pathways_dict.items()}

        # run cluster finding
        ## initial community detection
        cd = CommunityDetection()
        all_clusters_sets, level_mapping  = cd.get_communities_all_hierarchies(ppi_nx,True)

        # use pcst
        pct = PCST()

        # use GSEA on pcst results
        gi = GSEaInfr()

        print("#start")

        final_results = {}
        ge_ami = PipelinesCodebase.run_ami_detection(enrch_genes, non_enrch_genes, ppi_nx, all_pathways_dict,all_clusters_sets,
                                  level_mapping,cd,pct,gi)

        final_results["original"] = ge_ami

        with tqdm(n_permut) as p :
            for permu_index in range(n_permut) :
                _peru_enrch_genes, _peru_non_enrch_genes = PipelinesCodebase.build_permuted_sets(enrch_genes,non_enrch_genes)

                _peru_ge_ami = PipelinesCodebase.run_ami_detection(_peru_enrch_genes, _peru_non_enrch_genes, ppi_nx,
                                                                   all_pathways_dict,all_clusters_sets,level_mapping,cd,pct,gi)

                final_results[permu_index] = _peru_ge_ami
                p.update(1)

        return final_results


if __name__ == '__main__':
    E2eNoArmPerPpath.run_pipeline(ge_name = "ift" , ge_enrich_trh = 1.1 )
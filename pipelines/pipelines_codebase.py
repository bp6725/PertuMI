from DAL.ppis import PPIS
from DAL.gene_exp import GeneExp
from DAL.pathways import Pathways
from models.community_detection import CommunityDetection
from models.pcst import PCST
from models.gesa_infr import GSEaInfr
from tqdm import tqdm
import random
import itertools
from DAL.builder_synthetic_data import BuilderSyntheticData
from utils import Utils


class PipelinesCodebase():
    def __init__(self):
        pass

    @staticmethod
    def run_ami_detection(enrch_genes, non_enrch_genes, ppi_nx, all_pathways_dict, all_clusters_sets,
                          level_mapping, mutual_genes, cd, pct, gi):

        all_res = {}
        if type(enrch_genes[0]) is list:
            for ix, (_enrch_genes, _non_enrch_genes) in enumerate(zip(enrch_genes, non_enrch_genes)):
                ## background genes for the GSEA
                background = mutual_genes

                with tqdm(len(all_pathways_dict)) as p:
                    pathway_to_pval = {}
                    for pat_name, pat_genes in all_pathways_dict.items():
                        pval = PipelinesCodebase.run_over_single_pathway({pat_name: pat_genes}, _enrch_genes,
                                                                         all_clusters_sets,
                                                                         level_mapping, ppi_nx, background,
                                                                         cd=cd, pct=pct, gi=gi)
                        pathway_to_pval[pat_name] = pval
                        p.update(1)

                all_res[ix] = pathway_to_pval
        else:
            ## background genes for the GSEA
            background = set(ppi_nx.nodes).intersection(set(enrch_genes).union(set(non_enrch_genes)))

            # with tqdm(len(all_pathways_dict)) as p:
            pathway_to_pval = {}
            for pat_name, pat_genes in all_pathways_dict.items():
                pval = PipelinesCodebase.run_over_single_pathway({pat_name: pat_genes}, enrch_genes,
                                                                 all_clusters_sets,
                                                                 level_mapping, ppi_nx, background,
                                                                 cd=cd, pct=pct, gi=gi)
                if pval is None:
                    # pval = 1
                    # p.update(1)
                    continue

                pathway_to_pval[pat_name] = pval
                # p.update(1)

            all_res[0] = pathway_to_pval
        return all_res

    @staticmethod
    def run_over_single_pathway(pathway_dict, enrch_genes, all_clusters_sets,
                                level_mapping, ppi_nx, background,
                                cd=CommunityDetection(), pct=PCST(), gi=GSEaInfr()):
        pathway_name = list(pathway_dict.keys())[0]
        pathway_genes = pathway_dict[pathway_name]
        enrch_genes_in_pathway = set(pathway_genes).intersection(set(enrch_genes))

        if len(enrch_genes_in_pathway) < 10:
            return None

        sub_gs = cd.get_genes_cover_for_lower_level(enrch_genes_in_pathway, all_clusters_sets, level_mapping)[0]

        # print(f"{pathway_name}:{len(pathway_genes)} ; enrch_genes_in_pathway : {len(enrch_genes_in_pathway)} ")
        # if sub_gs is None :
        #     print("None")
        # else :
        #     print(
        #     f"sub_gs : {len(sub_gs)} ; interseaction_sub_gs : {len(set(sub_gs).intersection(set(enrch_genes_in_pathway)))} ")

        if sub_gs is None:
            return None

        sub_gs_nx = ppi_nx.subgraph(sub_gs)
        active_model, active_edge = pct.find_prize_tree(sub_gs_nx, enrch_genes)

        if len(active_model) < 10:
            return None

        # Important to notice that we use the HG only on "mutual genes" (= background)
        single_pathway = {pathway_name: pathway_dict[pathway_name]}
        return gi.compute_pval_local_single_pathway(active_model, single_pathway, background=background)

    @staticmethod
    def build_permuted_sets(enrch_genes, non_enrch_genes):
        all_genes = set(set(enrch_genes).union(set(non_enrch_genes)))
        permuted_enrch_genes = random.sample(all_genes, len(enrch_genes))
        permuted_non_enrch_genes = set(all_genes).difference(set(permuted_enrch_genes))

        return permuted_enrch_genes, permuted_non_enrch_genes

    @staticmethod
    def build_synthetic_experiment_sets(n_enriched_paths,perc_random_genes, ppi_name="huri", min_cluster_size=15):
        ## ppi network
        ppi_nx = PPIS.read_ppi_network(network=ppi_name, min_cluster_size=min_cluster_size, node_id_type="symbols")
        ppi_genes_set = set(ppi_nx.nodes)

        phws = Pathways()
        _all_pathways_dict = phws.get_gmt_file()

        # We filter the pathways by removing genes which are not in the ppi
        all_pathways_dict = {path_name: list(filter(lambda x: x in ppi_genes_set, path_genes))
                             for path_name, path_genes in _all_pathways_dict.items()}

        picked_pathways, picked_genes_per_pathway = BuilderSyntheticData.generate_enriched_ge_for_pathways(
            ppi_genes_set, all_pathways_dict, n_enriched_paths, perc_random_genes)

        enriched_genes = list(itertools.chain(*(picked_genes_per_pathway)))
        non_enriched_genes = list(set(ppi_genes_set).difference(set(enriched_genes)))

        print(picked_pathways)

        ppi_nx, (enriched_genes_set, non_enriched_genes_set), relevant_pathways_dict, background = \
            Utils.build_matched_sets(ppi_nx, [enriched_genes, non_enriched_genes], all_pathways_dict)

        cd = CommunityDetection()
        all_clusters_sets, level_mapping = cd.get_communities_all_hierarchies(ppi_nx, True)

        return enriched_genes_set, non_enriched_genes_set, relevant_pathways_dict, ppi_nx,\
               all_clusters_sets, level_mapping, background
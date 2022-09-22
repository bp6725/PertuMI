import numpy as np
import networkx as nx
import pandas as pd
from DAL.ppis import PPIS
from DAL.gene_exp import GeneExp
from DAL.pathways import Pathways
from models.community_detection import CommunityDetection
from models.pcst import PCST
from models.gesa_infr import GSEaInfr

# region Config


# endregion


# Load datasets
## ppi network
ppi_nx = PPIS.read_ppi_network(min_cluster_size=30,node_id_type="symbols")

## enriched genes
ge = GeneExp()
enrch_genes,non_enrch_genes = ge.return_gene_enrch_sets()

## background genes for the GSEA
background = set(ppi_nx.nodes).intersection(set(enrch_genes).union(set(non_enrch_genes)))

phws = Pathways()
all_pathways_dict = phws.get_gmt_file()

# for each pathway p find p.inters(ge)
enrch_genes_in_pathways = list(map(lambda x:set(x).intersection(set(enrch_genes)),all_pathways_dict.values()))

# run cluster finding
## initial community detection
cd = CommunityDetection()
all_clusters_sets, level_mapping  = cd.get_communities_all_hierarchies(ppi_nx,True)

##
all_sub_gs = cd.get_genes_cover_for_lower_level(enrch_genes_in_pathways,all_clusters_sets, level_mapping)[0]
all_sub_gs_nx = list(map(lambda x:ppi_nx.subgraph(x), all_sub_gs))

# use pcst
pct = PCST()
active_models_and_edges = pct.find_prize_tree(all_sub_gs_nx,enrch_genes)
active_models = list(map(lambda x:x[0], active_models_and_edges))

# use GSEA on pcst results
gi = GSEaInfr()
single_pathway = {phws_enr:all_clusters_sets[phws_enr]}
res = gi.compute_pval_local_single_pathway(all_genes,single_pathway, background = 1000)


print("pess")



print("end")

import pandas as pd
import networkx as nx
from utils import Utils

class PPIS():
    def __init__(self):
        pass

    #region Public
    @staticmethod
    def read_ppi_network(network="huri",node_id_type = "ensembl.gene",min_cluster_size = 10):
        raw_ppi = PPIS._read_raw_ppi(network)

        if node_id_type == "symbols" :
            raw_ppi = Utils.convert_ensembl_to_symbols_in_df(raw_ppi)

        G = nx.from_pandas_edgelist(raw_ppi)

        for component in list(nx.connected_components(G)):
            if len(component) < min_cluster_size:
                for node in component:
                    G.remove_node(node)

        return G

    #endregion


    #region Private

    @staticmethod
    def _read_raw_ppi(network="huri"):
        if network == "huri":
            pp_network = pd.read_table("../../data/From_DOMINO/PPIs/HuRI.tsv", header=None)
        if network == "dip":
            pp_network = pd.read_table("../../data/From_DOMINO/PPis/dip.tsv", header=None)
        if network == "string":
            pp_network = pd.read_table("../../data/From_DOMINO/PPis/string.sif", header=None)

        return pp_network.rename(columns={0: "source", 1: "target"})

    #endregion
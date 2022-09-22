import itertools
import mygene
import sys
import os
import pickle

class Utils():
    CACHE_PATH = "/home/benny/Repos/knowledge_base_clustering/pertuMi/cache/local_ens_to_syb_map.pkl"

    @staticmethod
    def build_matched_sets(ppi_network, genes_lists, pathways_dict):
        '''
        Making sure that PPI, Pathways and GE can play together.
        The idea : 1) PPI and GE should have the same genes(aka Background). 2) all genes in the pathways should be in
         Background but not the other way around.
        :param ppi_network:
        :param genes_lists:
        :param pathways_dict:
        :return: ppi_network, genes_lists, pathways_dict
        '''
        all_genes_from_lists = set(itertools.chain(*genes_lists))
        all_genes_from_ppi = set(ppi_network.nodes)

        backgroung_genes = all_genes_from_lists.intersection(all_genes_from_ppi)
        filtered_genes_list = [list(filter(lambda x: x in backgroung_genes, gene_list)) for gene_list in genes_lists]
        filtered_nx = ppi_network.subgraph(backgroung_genes)

        filtered_pathways_dict = {path_name: list(filter(lambda x: x in backgroung_genes, path_genes))
                                  for path_name, path_genes in pathways_dict.items()}
        return filtered_nx, tuple(filtered_genes_list), filtered_pathways_dict, backgroung_genes

    @staticmethod
    def convert_ensembl_to_symbols(list_of_genes_list, local_ens_to_syb_map = {}) :
        known_ens_to_syb_map = Utils._load_known_ensembl_to_symbols_map()
        current_ens_to_syb_map = {**local_ens_to_syb_map,**known_ens_to_syb_map}

        if type(list_of_genes_list[0]) is list:
            if type(list_of_genes_list[0][0]) is list:
                new_lists = []
                for l in list_of_genes_list :
                    _new_list,_ = Utils.convert_ensembl_to_symbols(l,current_ens_to_syb_map)
                    new_lists.append(_new_list)
                return new_lists, current_ens_to_syb_map

        all_genes = list(itertools.chain(*list_of_genes_list))
        unknown_genes = [g for g in all_genes if (not (g in current_ens_to_syb_map.keys()))]

        if len(unknown_genes) > 0 :
            mg = mygene.MyGeneInfo()
            out = mg.querymany(unknown_genes, scopes='ensembl.gene', fields='symbol', species='human', verbose=False)
            for line in out:
                try:
                    current_ens_to_syb_map[line["query"]] = line["symbol"]
                except KeyError:
                    current_ens_to_syb_map[line["query"]] = None
        Utils._dump_current_ensembl_to_symbols_map(current_ens_to_syb_map)
        new_lists = [[current_ens_to_syb_map[g] for g in list if g in current_ens_to_syb_map.keys()]
                     for list in list_of_genes_list]
        return new_lists, current_ens_to_syb_map

    @staticmethod
    def convert_ensembl_to_symbols_in_df(df) :
        all_genes_in_df = set(list(itertools.chain(*[ele[1].to_list() for ele in df.items()])))
        _,gene_to_sym_map = Utils.convert_ensembl_to_symbols([all_genes_in_df])

        new_df = df.applymap(lambda x:gene_to_sym_map[x] if x in gene_to_sym_map.keys() else None )
        return new_df.dropna()

    @staticmethod
    def _load_known_ensembl_to_symbols_map():
        local_ens_to_syb_map = {}
        if os.path.exists(Utils.CACHE_PATH):
            with open(Utils.CACHE_PATH, 'rb') as f:
                local_ens_to_syb_map = pickle.load(f)

        return local_ens_to_syb_map

    @staticmethod
    def _dump_current_ensembl_to_symbols_map(current_mapping):
        with open(Utils.CACHE_PATH, 'wb') as f:
            pickle.dump(current_mapping,f)



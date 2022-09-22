import numpy as np
import networkx as nx
import pcst_fast



class PCST():
    def __init__(self):
        pass

    #region Public

    def find_prize_tree(self,sub_g_nx, rel_genes, prize_val=3, cost_val=0.1):
        '''
        Find the pcst subtree. If 'sub_g_nx' or 'rel_genes' is list, run this method on every element
        :param sub_g_nx:
        :param rel_genes:
        :param prize_val:
        :param cost_val:
        :return:
        '''
        if type(sub_g_nx) is list :
            return list(map(lambda x:self.find_prize_tree(x[0],x[1],prize_val,cost_val),zip(sub_g_nx,rel_genes)))

        # the package required indexs insted of gene names
        node_to_number = {node: i for i, node in enumerate(sub_g_nx.nodes)}
        edges = np.array([(node_to_number[edge[0]], node_to_number[edge[1]]) for edge in sub_g_nx.edges])

        # set prizes and costs
        prizes = np.array([prize_val * int(node in rel_genes) for node in node_to_number.keys()])
        costs = np.array([cost_val for node in edges])

        # run PCST
        _vertices_res, _edges_res = pcst_fast.pcst_fast(edges, prizes,
                                                        costs, -1, 1, "strong", 0)

        # move back to gene names
        number_to_nodes = {a: b for b, a in node_to_number.items()}
        vertices_res = list(map(lambda x: number_to_nodes[x], _vertices_res))
        edges_res = [(number_to_nodes[edges[e][0]], number_to_nodes[edges[e][1]]) for e in _edges_res]

        return vertices_res, edges_res

    #endregion

    #region Private



    #endregion
from community import community_louvain
from collections import defaultdict
import numpy as np
import itertools


class CommunityDetection():
    def __init__(self):
        pass

    # region Public

    def get_communities_all_hierarchies(self, G, return_level_mapping=False):
        # build all communites by hier
        dendrogram = community_louvain.generate_dendrogram(G, resolution = 10)

        # build all clusters sets
        return self._return_all_cluster_sets(dendrogram, return_level_mapping)

    @staticmethod
    def get_smallest_cluster_over_score(gene_set, all_clusters_sets, score_trh=0.7):
        '''
        Find the smallest cluster with intersection score(IS) higher than score_trh.
        If there is no cluster with IS bigger then score_trh, lower the score trh recrosuvly.
        '''
        clusters_intersection_results = CommunityDetection._return_intersection_per_cluster(gene_set, all_clusters_sets,
                                                                                            intersection_trh=0.01)

        return CommunityDetection._get_smallest_cluster_over_score(clusters_intersection_results, score_trh)

    @staticmethod
    def get_genes_cover_for_lower_level(set_to_cover, all_clusters_sets,
                                        level_mapping, intersection_with_gs_trh=0.01,
                                        gs_cover_trh=0.7):
        '''
        If set_to_cover is a list run this function over each item in the list and return list of results.
        Find all the clusters with intersection score of more than intersection_with_gs_trh with set_to_cover.
        First we start with all clusters of level 0, if the final intersection score of all clusters with set_to_cover
        is <  gs_cover_trh we add clusters from higher level.
        :param set_to_cover:
        :param all_clusters_sets:
        :param level_mapping:
        :param intersection_with_gs_trh:
        :param gs_cover_trh:
        :return:
        '''
        if type(set_to_cover) is list :
            return list(map(lambda sin_set_to_cover:
                            CommunityDetection.get_genes_cover_for_lower_level(sin_set_to_cover,all_clusters_sets,
                                        level_mapping, intersection_with_gs_trh,gs_cover_trh), set_to_cover))

        if len(set_to_cover) == 0:
            return None, None

        intersection_score = lambda x: len(x[0].intersection(x[1])) / min(len(x[0]), len(x[1])) if min(len(x[0]), len(
            x[1])) > 0 else 0

        all_inter_clusters = []
        all_inter_clusters_genes = []
        # we try the smallest resolution last, because it will fit everything
        for li in list(range(1,len(level_mapping.keys()))) + [0]:
            all_clusters_in_level = [(c, all_clusters_sets[c]) for c in level_mapping[li]]

            for c_name, c_genes in all_clusters_in_level:
                if intersection_score((set_to_cover, c_genes)) > intersection_with_gs_trh:
                    all_inter_clusters_genes.append(c_genes)
                    all_inter_clusters.append(c_name)

            # cover set so far
            curr_cover_set = set(itertools.chain(*all_inter_clusters_genes))
            curr_cover_set_score = intersection_score((set_to_cover, curr_cover_set))

            if curr_cover_set_score > gs_cover_trh: return (curr_cover_set, curr_cover_set_score, all_inter_clusters)

        return None,None

    # endregion

    # region Private

    def _return_all_cluster_sets(self, dendrogram, return_level_mapping=False):
        all_clusters = defaultdict(list)
        level_mapping = defaultdict(list)
        for li in range(len(dendrogram)):
            level = community_louvain.partition_at_level(dendrogram, li)
            for g, clu in level.items():
                all_clusters[(li, clu)].append(g)
                level_mapping[li].append((li, clu))

        all_clusters_sets = {k: set(v) for k, v in all_clusters.items()}
        level_mapping = {k: set(v) for k, v in level_mapping.items()}
        return all_clusters_sets if (not return_level_mapping) else (all_clusters_sets, level_mapping)

    @staticmethod
    def _return_intersection_per_cluster(gene_set, all_clusters_sets, intersection_trh=0.1):
        intersection_score = lambda x: len(x[0].intersection(x[1])) / len(x[0])

        # keep all clusters with intersection greater than intersection_trh
        inter_clu_size = list(map(lambda x: (x[0], x[0][0], intersection_score((gene_set, x[1])), len(x[1])),
                                  filter(lambda x: intersection_score((gene_set, x[1])) > intersection_trh,
                                         all_clusters_sets.items())))

        # each row is : [cluster name,cluster level,intersection score,cluster_size]
        return np.array(inter_clu_size, dtype=object)

    @staticmethod
    def _get_smallest_cluster_over_score(clusters_intersection_results, score_trh):
        '''
        This is the recursion function.
        Find the smallest cluster with intersection score(IS) higher than score_trh.
        If there is no cluster with IS bigger than score_trh, lower the score trh recursively.
        '''
        smallest_size = 20000
        smallest_cluster_name = None
        smallest_cluster_trh = 1
        for c_name, level, i_score, size in clusters_intersection_results:
            if i_score < score_trh: continue
            if size < smallest_size:
                smallest_size = size
                smallest_cluster_name = c_name
                smallest_cluster_trh = i_score

        if smallest_cluster_name is not None:
            return smallest_cluster_name, smallest_cluster_trh, smallest_size
        else:
            return CommunityDetection._get_smallest_cluster_over_score(clusters_intersection_results, score_trh - 0.1)

    # endregion
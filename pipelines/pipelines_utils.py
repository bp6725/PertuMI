from collections import Counter


class PipelinesUtils():

    @staticmethod
    def return_best_am_for_each_permut(final_results):
        best_of_each_permut = list(
            map(lambda x: (x[0][min(x[0], key=x[0].get)], min(x[0], key=x[0].get)), final_results.values()))
        fdr_pval = sum([(pval<best_of_each_permut[0][0]) for pval,pat in best_of_each_permut[1:]])/len(final_results)

        path_counter = Counter([pat for pval, pat in best_of_each_permut[1:]])
        original_count = Counter([pat for pval, pat in best_of_each_permut[1:]])[best_of_each_permut[0][1]]

        return best_of_each_permut, fdr_pval, path_counter, original_count
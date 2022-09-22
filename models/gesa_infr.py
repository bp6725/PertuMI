import sklearn.feature_selection as skf
import numpy as np
import gseapy
from gseapy.enrichr import Enrichr
import mygene
import subprocess


class GSEaInfr():
    def __init__(self):
        pass

    #region Public

    @staticmethod
    def compute_neg_log_gsea_p_value(pathways, result_genes,cutoff=0.05):
        """Computes gene set enrichment score for result genes w.r.t. phenotype-related pathways.
        Parameters
        ----------
        pathways : list of str
            Names of phenotype-related pathways.
        result_genes : list of str
            Set of genes computed by network enrichment algorithm.
        Returns
        -------
        mean neg_log_gsea_p_value : float
            Negative log-transformed p-value of gene set enrichment analysis.
        """
        if len(result_genes) == 0:
            return 0.0
        try:
            res = gseapy.enrichr(gene_list=result_genes, description='pathway', gene_sets=pathways, cutoff=cutoff,
                                 outdir='../temp/enrichment', no_plot=True)
            full_results = res.results
            terms = list(full_results.Term)
            terms = [x.split(' ')[-1] for x in terms]
            p_values = []
            for i in range(len(terms)):
                p_values.append((terms[i], -np.log10(full_results['Adjusted P-value'][i])))
            subprocess.call('rm -rf ../temp/enrichment/', shell=True)
            if len(p_values) > 0:
                return p_values
            else:
                return []
        except:
            return -1

    @staticmethod
    def compute_pval_local_single_pathway(result_genes,pathway,background,cutoff=0.05):
        # Important to notice that we use the HG only on "mutual genes" (= background)
        result_genes = list(filter(lambda x : x in background,result_genes))
        pathway = {k:list(filter(lambda x : x in background,v)) for k,v in pathway.items()}

        enr = Enrichr(result_genes, pathway, organism='human',outdir='Enrichr', background=background, cutoff = cutoff,
                      format='pdf', figsize=(8, 6), top_term=10, no_plot=False, verbose=False)
        enr.set_organism()
        enr.parse_genelists()
        if enr.enrich(pathway) is None :
            return 1
        return enr.enrich(pathway)["P-value"].astype(float).iloc[0]

    #endregion


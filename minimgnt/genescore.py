#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.stats
import six

from .utils import *
from .stepwisefit import *

def _calc_uncorr_gene_score(gene, input_gene, input_snp, pruned_snps, hotspots):
    # find local snps given a gene
    cond_snps_near_gene = logical_and(np.equal(input_snp[:, 0], input_gene[gene, 0]),
                                      np.greater_equal(input_snp[:, 1], input_gene[gene, 1]),
                                      np.less_equal(input_snp[:, 1], input_gene[gene, 2]))
    # if no snps found
    if not np.any(cond_snps_near_gene):
        return (np.nan, 0, 1, 0, 0)

    n_snps_zscore_finite = np.sum(np.isfinite(input_snp[cond_snps_near_gene][:, 3]))
    # if no snps with finite zcore
    if n_snps_zscore_finite == 0:
        return (np.nan, 0, 1, 0, 0)

    n_snps_per_gene = n_snps_zscore_finite

    # use p-value to find most significant SNP
    idx_min_pval = np.nanargmin(input_snp[cond_snps_near_gene][:, 3])

    uncorr_score = input_snp[cond_snps_near_gene][idx_min_pval, 2]

    # count number of independent SNPs per gene
    n_indep_snps_per_gene = np.sum(logical_and(np.equal(pruned_snps[:, 0], input_gene[gene, 0]),
                                               np.greater_equal(pruned_snps[:, 1], input_gene[gene, 1]),
                                               np.less_equal(pruned_snps[:, 1], input_gene[gene, 2])))

    # count number of hotspots per gene
    n_hotspots_per_gene = np.sum(np.logical_and(np.equal(hotspots[:, 0], input_gene[gene, 0]),
                                                np.greater(np.fmin(hotspots[:, 2], input_gene[gene, 2])
                                                           - np.fmax(hotspots[:, 1], input_gene[gene, 1]), 0)))
    return (uncorr_score, n_snps_per_gene, 0, n_indep_snps_per_gene, n_hotspots_per_gene)


def calc_uncorr_gene_score(input_gene, input_snp, pruned_snps, hotspots, n_cpus):
    f = lambda gene: _calc_uncorr_gene_score(gene, input_gene, input_snp, pruned_snps, hotspots)

    if n_cpus > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=n_cpus) as e:
            ret = e.map(f, six.moves.range(len(input_gene)))
    else:
        ret = (f(gene) for gene in six.moves.range(len(input_gene)))

    uncorr_score, n_snps_per_gene, n_genes_score_nan, n_indep_snps_per_gene, n_hotspots_per_gene = [np.array(x) for x in zip(*ret)]

    return (np.fabs(uncorr_score), n_snps_per_gene, np.sum(n_genes_score_nan), n_indep_snps_per_gene, n_hotspots_per_gene)


def calc_corr_gene_score(confounders, uncorr_score, cutoff=0.05):
    n_genes, n_confounders = confounders.shape
    cond_score_finite = np.isfinite(uncorr_score)
    beta, se, pval, inmodel, stats, nextstep, history = stepwisefit(confounders[cond_score_finite], uncorr_score[cond_score_finite])
    residuals = uncorr_score - stats.intercept
    for j in six.moves.range(n_confounders):
        if pval[j] <= cutoff:
            residuals = residuals - beta[j] * confounders[:, j]
    # residuals must be finite for np.mean/std
    #                   not nan for cdf (avoid warnings)
    cond_residuals_finite = np.isfinite(residuals)
    cond_residuals_nnan = np.logical_not(np.isnan(residuals))
    corr_score = np.tile(np.nan, len(residuals))
    corr_score[cond_residuals_nnan] = 1 - sp.stats.norm.cdf(residuals[cond_residuals_nnan],
                                                            np.mean(residuals[cond_residuals_finite]),
                                                            np.std(residuals[cond_residuals_finite]))
    return (corr_score, residuals)


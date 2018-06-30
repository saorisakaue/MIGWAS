#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse
import os
import time

from minimgnt.genescore import *
from minimgnt.utils import *


__version__ = '0.9'

BASEDIR = os.path.dirname(__file__)
DATADIR = os.path.join(BASEDIR, 'minimgnt', 'data')
REFSEQ_GENE_FILEPATH = os.path.join(DATADIR, 'AllHumanGeneChrPosStrandNames_RefSeq_hg19_072111.txt')
MIRBASE_MIR_FILEPATH = os.path.join(DATADIR, 'AllHumanMiRNAChrPosStrandNames_miRBase_hg19_062413.txt')
PRUNED_SNPS_FILEPATH = os.path.join(DATADIR, 'CEU_HapMap_pruned_SNPs_ChrNumPos_hg19_072111.txt')
HOTSPOT_BOUNDARIES_FILEPATH = os.path.join(DATADIR, 'hotspot_boundaries_b37_hg19_072111.txt')
OUTMIRDIR = os.path.join(BASEDIR, 'miRNA_P')
OUTGENEDIR = os.path.join(BASEDIR, 'Gene_P')

def minimgnt(args):
    logger = Logger(args.out + '.log')
    logger.log('****************************************************************\n')
    logger.log('* minimgnt (miniMAGENTA) ver. {V} written by Masahiro Kanai.\n', V=__version__)
    logger.log('* https://github.com/mkanai/minimgnt\n')
    logger.log('*\n')
    logger.log('* MAGENTA was written by Ayellet Segre, Altshuler and Daly Labs.\n')
    logger.log('* https://www.broadinstitute.org/mpg/magenta/\n')
    logger.log('****************************************************************\n\n')

    logger.log('Call: {S}\n', S=' '.join(sys.argv))
    logger.log('Analysis started at {T}\n\n', T=time.ctime())

    # load gene/mir files
    # column 0: chrom
    #        1: txStart (bp)
    #        2: txEnd (bp)
    #        3: strand (1: forward, 0: reverse)
    #        4: name
    refseq_gene = pd.read_csv(REFSEQ_GENE_FILEPATH, header=None, delim_whitespace=True)
    n_genes = len(refseq_gene.index)
    logger.log('{N} genes loaded from {F}\n', N=n_genes, F=REFSEQ_GENE_FILEPATH)

    mirbase_mir = pd.read_csv(MIRBASE_MIR_FILEPATH, header=None, delim_whitespace=True)
    n_mirs = len(mirbase_mir.index)
    logger.log('{N} miRNAs loaded from {F}\n\n', N=n_mirs, F=MIRBASE_MIR_FILEPATH)

    if args.remove_HLA:
        refseq_gene = remove_HLA_region(refseq_gene, args.HLA_start, args.HLA_end)
        mirbase_mir = remove_HLA_region(mirbase_mir, args.HLA_start, args.HLA_end)
        n_genes2, n_mirs2 = len(refseq_gene.index), len(mirbase_mir.index)
        logger.log('--remove-HLA: {G} genes and {M} miRNAs were removed.\n',
                   G=n_genes - n_genes2, M=n_mirs - n_mirs2)
        n_genes, n_mirs = n_genes2, n_mirs2
        logger.log('              {G} genes and {M} miRNAs will be used in the following analysis.\n\n',
                   G=n_genes, M=n_mirs)

    # rbind gene and mir
    all_gene_mir = np.vstack((refseq_gene.iloc[:, :4].values, mirbase_mir.iloc[:, :4].values))
    all_gene_mir = extend_gene_site(all_gene_mir, args.boundr_upstr, args.boundr_downstr)

    # load other files
    pruned_snps = pd.read_csv(PRUNED_SNPS_FILEPATH, header=None, delim_whitespace=True).values
    logger.log('Pruned SNPs loaded from {F}\n', F=PRUNED_SNPS_FILEPATH)
    hotspot_boundaries = pd.read_csv(HOTSPOT_BOUNDARIES_FILEPATH, header=None, delim_whitespace=True).values
    logger.log('Hotspot boundaries loaded from {F}\n\n', F=HOTSPOT_BOUNDARIES_FILEPATH)

    # load GWAS SNP score file
    # column X0: rsID (omitted)
    #         0: chr
    #         1: position
    #         2: z-score (optional)
    #         3: p-value
    input_snp = pd.read_csv(args.score_filename, header=None, delim_whitespace=True)
    if args.rsid:
        input_snp = input_snp.iloc[:, 1:].values
    else:
        input_snp = input_snp.values

    n_snps, n_c = input_snp.shape
    logger.log('{N} SNPs loaded from {F}\n\n', N=n_snps, F=args.score_filename)

    if n_c < 3:
        raise ValueError('Cannot find p-value in {F}'.format(F=args.score_filename))
    elif n_c == 3:
        # if pval only
        if np.nanmax(input_snp[:, 2]) <= 1:
            # remove NaN
            input_snp = input_snp[np.isfinite(input_snp[:, 2])]
            input_snp = convert_twotailed_pval_to_zscore(input_snp)
            n_snps2 = len(input_snp)
            logger.log('Converted {N} SNP p-values to z-score ({M} SNPs were removed due to NA)\n\n',
                       N=n_snps2, M=n_snps - n_snps2)
            n_snps = n_snps2
        else:
            raise ValueError('p-values must be 0 <= p <= 1')

    # calc uncorrected score from local SNP association z-scores Best SNP per gene scoring metric
    uncorr_score, n_snps_per_gene, n_genes_score_nan, n_indep_snps_per_gene, n_hotspots_per_gene = calc_uncorr_gene_score(
        all_gene_mir, input_snp, pruned_snps, hotspot_boundaries, args.cpus)
    logger.log('Calculated uncorrected "gene association score" from best SNP association z-scores per gene\n')
    logger.log('{N} out of {T} genes/mirs were not assigned the score.\n\n',
               N=n_genes_score_nan, T=n_genes + n_mirs)

    # calc corrected score according to confounders using step-wise multivariate linear regression analysis
    gene_size_plus_interval_kb = (all_gene_mir[:, 2] - all_gene_mir[:, 1]) / 1000.0

    confounders = np.vstack((gene_size_plus_interval_kb,
                             n_snps_per_gene / gene_size_plus_interval_kb,
                             n_indep_snps_per_gene / gene_size_plus_interval_kb,
                             n_hotspots_per_gene / gene_size_plus_interval_kb)).T

    corr_score, residuals = calc_corr_gene_score(confounders, uncorr_score)
    logger.log('Calculated corrected "gene association score" using step-wise multivariate linear regression analysis\n\n')

    # output results
    n_genes = len(refseq_gene)
    gene_result = pd.DataFrame({'GENE': refseq_gene.iloc[:, 4], 'P': corr_score[:n_genes]})
    mir_result = pd.DataFrame({'MIR': mirbase_mir.iloc[:, 4], 'P': corr_score[n_genes:]})

    if args.remove_NA:
        gene_result = gene_result.dropna()
        mir_result = mir_result.dropna()

    gene_result.to_csv(os.path.join(OUTGENEDIR, args.out + '.gene.pval.txt'), sep='\t', na_rep='NA', header=False, index=False, float_format='%.6g')
    mir_result.to_csv(os.path.join(OUTMIRDIR, args.out + '.mir.pval.txt'), sep='\t', na_rep='NA', header=False, index=False, float_format='%.6g')

    logger.log('Results were written to {F}.{{gene, mir}}.pval.txt\n\n', F=args.out)
    logger.log('Analysis finished at {T}.\n', T=time.ctime())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate the corrected gene association score from a GWAS result, according to MAGENTA\'s method.')
    parser.add_argument('score_filename', type=str,
                        help='GWAS SNP score filename containing one SNP information per row. '
                             '[columns] (1) rsID, (2) chromosome, (3) bp, (4) z-score (optional), (5) p-value.')
    parser.add_argument('--out', default='your_phenotype', type=str, help='output filename prefix. This should preferably be a phenotype name.')
    parser.add_argument('-j', '--cpus', default=1, type=int, help='a number of cpus used for computation')
    parser.add_argument('--not-remove-HLA', dest='remove_HLA', action='store_false',
                        help='do not remove genes in HLA region from a result. (default: False)')
    parser.add_argument('--HLA-start', default=25000000, type=int,
                        help='start position of HLA region in chr6. (default: 25,000,000)')
    parser.add_argument('--HLA-end', default=35000000, type=int,
                        help='end position of HLA region in chr6. (default: 35,000,000)')
    parser.add_argument('--boundary-upstream', dest='boundr_upstr', default=110000, type=int,
                        help='added distance (bp) upstream to gene\'s start position. (default: 110,000)')
    parser.add_argument('--boundary-downstream', dest='boundr_downstr', default=40000, type=int,
                        help='added distance (bp) downstream to gene\'s end position. (default: 40,000)')
    parser.add_argument('--remove-NA', action='store_true',
                        help='remove genes with NA score from the output (default: False)')
    parser.add_argument('--no-rsid', dest='rsid', action='store_false',
                        help='use this flag when a score file doesn\'t contain a rsID column. (default: False)')

    args = parser.parse_args()

    if args.cpus > 1:
        try:
            import concurrent.futures
        except ImportError:
            raise ImportError('For multithreading, please install "futures" package via `pip install futures`')
    if args.HLA_start < 0:
        raise ValueError('--HLA-start must be a positive intenger.')
    if args.HLA_end < 0:
        raise ValueError('--HLA-end must be a positive intenger.')
    if args.boundr_upstr < 0:
        raise ValueError('--boundary-upstream must be a positive intenger.')
    if args.boundr_downstr < 0:
        raise ValueError('--boundary-downstream must be a positive intenger.')

    minimgnt(args)


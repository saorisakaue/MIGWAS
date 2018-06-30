#!/usr/bin/env python
# -*- coding: utf-8 -*-

import functools
import numpy as np
import scipy as sp
import scipy.stats
import sys


class Logger(object):

    # bufsize=1 to write in a moment
    def __init__(self, filename, bufsize=1):
        self.file = open(filename, 'w', bufsize)

    def log(self, msg, *args, **kwargs):
        msg = msg.format(*args, **kwargs)
        self.file.write(msg)
        sys.stdout.write(msg)


def logical_and(*args):
    return functools.reduce(np.logical_and, args)


def logical_or(*args):
    return functools.reduce(np.logical_or, args)


def remove_HLA_region(input_gene, HLA_start, HLA_end):
    cond_not_HLA = np.logical_or(np.not_equal(input_gene.iloc[:, 0], 6),
                                 np.less_equal(np.fmin(input_gene.iloc[:, 2], HLA_end)
                                               - np.fmax(input_gene.iloc[:, 1], HLA_start), 0))
    return input_gene[cond_not_HLA]


def extend_gene_site(input_gene, boundr_upstr, boundr_downstr):
    # whether positive strand or not
    strand = np.equal(input_gene[:, 3], 1)
    input_gene[:, 1] -= np.where(strand, boundr_upstr, boundr_downstr)
    input_gene[:, 2] += np.where(strand, boundr_downstr, boundr_upstr)

    return input_gene


def convert_twotailed_pval_to_zscore(input_snp):
    pval = input_snp[:, 2]
    zscores = np.sqrt(sp.stats.chi2.ppf(1 - pval, 1))
    return np.concatenate((input_snp[:, :2], zscores[:, None], input_snp[:, 2:None]), axis=1)

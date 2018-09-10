#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import gzip
import pandas as pd
import numpy as np
from itertools import repeat
from scipy.stats import norm
import math
from multiprocessing import Pool
import argparse
import warnings; warnings.filterwarnings('ignore')

# consts
BASEDIR = os.path.dirname(__file__)
DATADIR = os.path.normpath(os.path.join(BASEDIR, 'data'))
MIRPDIR = os.path.normpath(os.path.join(BASEDIR, 'miRNA_P'))
GENEPDIR = os.path.normpath(os.path.join(BASEDIR, 'gene_P'))

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--phenotype', '-p', default=None, type=str,
    help='Name of the phenotype of interest (file name prefix from minimgnt output).',
    required=True)
parser.add_argument('--out', '-o', default="your_migwas", type=str,
    help='Output file prefix.')
parser.add_argument('--cpus', '-j', default=1, type=int,
    help='Number of CPUs to be used.')
parser.add_argument('--tsi', '-t', default=0.7, type=float,
    help='Tissue specificity index threshold for partitioning miRNA enrichment signal.',
    required=False)
parser.add_argument('--iterations', '-i', default=20000, type=int,
    help='Number of iterations to simulate null distributions.')
parser.add_argument('--output-candidate', '-c', default=False, action='store_true',
    help='If you want to output a list of candidate miRNAs and genes associated with the trait, set this flag.',
    required=False)
args = parser.parse_args()
iteration = args.iterations
num_threads = int(args.cpus)
MIRFILE = os.path.join(MIRPDIR,args.phenotype+".mir.pval.txt")
GENEFILE = os.path.join(GENEPDIR,args.phenotype+".gene.pval.txt")
# define arguments and values
EXPRESSION_DATA = os.path.join(DATADIR,"qn_exp.csv")
TISSUELIST = os.path.join(DATADIR,"human.srna.cellontology_nohead_rev.tsv")
IDFILE = os.path.join(DATADIR,"ID_MIMAT_MI.txt")
Exp = pd.read_csv(EXPRESSION_DATA, index_col=0)
Top = 90
mi_alpha = 0.01
ge_alpha = 0.01
thresflag = 5
THRES_list = {}
Algo_list = ["DIANAv5", "miRDB", "PITA_0_0", "TSH_Ncons"]
THRES_list["DIANAv5"] = [0.89, 0.916, 0.937, 0.955, 0.968, 0.979, 0.987, 0.993, 0.997]
THRES_list["miRDB"] = [65.359, 70.0753, 74.3884, 78.3834284889, 82.1096, 85.6510723331, 88.9973910801, 92.2397717152, 95.2520361716]
THRES_list["PITA_0_0"] = [66.12, 67.13, 68.22, 69.42, 70.84, 72.72, 75.83]
THRES_list["TSH_Ncons"] = [0.324, 0.358, 0.391, 0.424, 0.458, 0.496, 0.538, 0.593, 0.675]
percentile_list = (99,99.25010579,99.43765867,99.5783035,99.68377223,99.76286263,99.82217206,99.86664786,99.9)
# input
P_RNA = {}
with open(MIRFILE, 'rt') as f:
	for line in f:
		line = line.rstrip().split("\t")
		if line[1] != "NA":
			P_RNA[line[0]] = float(line[1]) # pre-hsa to P value
P_Gene = {}
with open(GENEFILE, 'rt') as f:
	for line in f:
		line = line.rstrip().split("\t")
		if line[1] != "NA":
			P_Gene[line[0]] = float(line[1]) # gene symbol to P value
pre_hsa = {}
MIMAT = {}
with open(IDFILE, 'rt') as f:
	for line in f:
		line = line.rstrip().split("\t")
		if line[4] != "-":
			pre_hsa[line[0]] = line[4] # MIMAT to pre-hsa
			MIMAT[line[2]] = line[0] # mature-hsa to  MIMAT
# get mean value for each tissue and tissue-sample dict
Matrix_mean_exp = pd.DataFrame({'Dummy':[i for i in repeat(0, len(Exp.index))]},index=Exp.index)
with open(TISSUELIST, 'rt') as f:
	for line in f:
		line = line.rstrip().split("\t")
		tissue_name = line[0]
		sample_list = []
		sample_list = list(line[1].split(","))
		Sample_exp = Exp.loc[:,sample_list]
		Sample_mean_exp = pd.DataFrame({tissue_name:Sample_exp.mean(axis=1)},index=Sample_exp.index)
		Matrix_mean_exp = pd.concat([Matrix_mean_exp, Sample_mean_exp],axis=1)
Matrix_mean_exp = Matrix_mean_exp.drop("Dummy",axis=1)
# get TSI and TSI significant miRNA list
max_row = Matrix_mean_exp.apply(max,axis=1)
max_row = max_row.values
Exp_diff = 1-Matrix_mean_exp.div(max_row,axis=0)
TSI = pd.Series(Exp_diff.apply(sum,axis=1)/(len(Matrix_mean_exp.columns)-1),index=Exp.index)
#get tissue specific miRNA list and exp matrix
Sig_TSI = TSI[TSI>args.tsi]
Sig_TSI_miRNA = list(Sig_TSI.index)
Analyze_miRNA = []
for mirna in Sig_TSI_miRNA:
	if mirna in MIMAT:
		if MIMAT[mirna] in pre_hsa:
			if pre_hsa[MIMAT[mirna]] in P_RNA:
				Analyze_miRNA.append(mirna)
# make a tissue specific high exp miRNA list (by tissue, only for TSI high miRNAs)
Specific_exp = Matrix_mean_exp.loc[Analyze_miRNA,:]
Ranked_exp = Specific_exp.rank()
Tissue_sig_MIMAT = {}
for tissue in Ranked_exp.columns:
	Rank = pd.Series(Ranked_exp.loc[:,tissue],index=Ranked_exp.index)
	Sig_exp = Rank[Rank>=len(Ranked_exp.index)*(Top/100)]
	for mature_miR in Sig_exp.index:
		Tissue_sig_MIMAT.setdefault(tissue, []).append(MIMAT[mature_miR])
# define functions
def np_random_shuffle(some_dic):
	keys = list(some_dic.keys())
	np.random.shuffle(keys)
	return dict(zip(keys, some_dic.values()))

def migwas_func(database):
	TARGETFILE = os.path.join(DATADIR,database+"_miRNAEnrichment.txt.gz")
	Target = pd.read_csv(TARGETFILE,delimiter='\t',compression='gzip',index_col=0)
	Target_thres = THRES_list[database]
	Target_mir = []
	Target_sym = []
	for  mimat in Target.columns:
		if mimat in pre_hsa:
			if pre_hsa[mimat] in P_RNA:
				Target_mir.append(mimat)
	for symbol in Target.index:
		if symbol in P_Gene:
			Target_sym.append(symbol)
	Target_ana = Target.loc[Target_sym,Target_mir]
	Tissue_sig_MIMAT["all"] = Target_mir
	Count_tissue = {}
	Count_i_tissue = {}
	P_count_tissue = {}
	SUM_iter_count_tissue = {}
	Pval_tissue = {}
	Fold_tissue = {}
	Mean_tissue = {}
	Top_assoc_pair = []
	Tissue_sig_MIMAT_algo = {}
	for tissue in Tissue_sig_MIMAT.keys():
		Count_tissue[tissue] = np.zeros((len(Target_thres)))
		Count_i_tissue[tissue] = []
		P_count_tissue[tissue] = [i for i in repeat(0,len(Target_thres))]
		SUM_iter_count_tissue[tissue] = [i for i in repeat(0,len(Target_thres))]
		Pval_tissue[tissue] = [i for i in repeat(0,len(Target_thres))]
		Fold_tissue[tissue] = [i for i in repeat(0,len(Target_thres))]
		Mean_tissue[tissue] = [i for i in repeat(0,len(Target_thres))]
		for MIMAT in Tissue_sig_MIMAT[tissue]:
			if MIMAT in Target_mir:
				Tissue_sig_MIMAT_algo.setdefault(tissue, []).append(MIMAT)
	Dis_sig_mir = [x for x in Target_ana.columns if P_RNA[pre_hsa[x]] < mi_alpha]
	Dis_sig_gene = [y for y in Target_ana.index if P_Gene[y] < ge_alpha]
	Dis_sig_Target = Target_ana.loc[Dis_sig_gene,Dis_sig_mir]
	for t in range(len(Target_thres)):
		Top_assoc_pair.append([])
		tmp = Dis_sig_Target >= Target_thres[t]
		thresflag = tmp.astype(int)
		for tissue in Tissue_sig_MIMAT_algo.keys():
			tis_mirna = [mirna for mirna in Tissue_sig_MIMAT_algo[tissue] if mirna in thresflag.columns]
			tis_thresflag = thresflag.loc[:,tis_mirna]
			Count_tissue[tissue][t] = tis_thresflag.sum(axis=0).sum()
		thresflag2 = thresflag.as_matrix()
		get_id = np.where(thresflag2==1)
		for id_num in range(len(get_id[0])):
			symbol = thresflag.index[get_id[0][id_num]]
			mimat = thresflag.columns[get_id[1][id_num]]
			Top_assoc_pair[t].append((symbol,pre_hsa[mimat],mimat))
	for j in range(iteration):
		iter_P_RNA = np_random_shuffle(P_RNA)
		iter_P_Gene = np_random_shuffle(P_Gene)
		Dis_sig_mir = [z for z in Target_ana.columns if iter_P_RNA[pre_hsa[z]] < mi_alpha]
		Dis_sig_gene = [w for w in Target_ana.index if iter_P_Gene[w] < ge_alpha]
		Iter_Target = Target_ana.loc[Dis_sig_gene,Dis_sig_mir]
		for tissue in Tissue_sig_MIMAT_algo.keys():
			Count_i_tissue[tissue].append(np.zeros((len(Target_thres))))
		for t_j in range(len(Target_thres)):
			tmp2 = Iter_Target >= Target_thres[t_j]
			thresflag2 = tmp2.astype(int)
			for tissue in Tissue_sig_MIMAT_algo.keys():
				tis_mirna2 = [mirna for mirna in Tissue_sig_MIMAT_algo[tissue] if mirna in thresflag2.columns]
				tis_thresflag2 = thresflag2.loc[:,tis_mirna2]
				Count_i_tissue[tissue][j][t_j] = tis_thresflag2.sum(axis=0).sum()
		for t in range(len(Target_thres)):
			for tissue in Tissue_sig_MIMAT_algo.keys():
				SUM_iter_count_tissue[tissue][t] += Count_i_tissue[tissue][j][t]
				if Count_i_tissue[tissue][j][t] >= Count_tissue[tissue][t]:
					P_count_tissue[tissue][t] += 1
	for t in range(len(Target_thres)):
		for tissue in Tissue_sig_MIMAT_algo.keys():
			Mean_tissue[tissue][t] = SUM_iter_count_tissue[tissue][t]/iteration
			Fold_tissue[tissue][t] = Count_tissue[tissue][t]/(Mean_tissue[tissue][t]+0.000000001)
			Pval_tissue[tissue][t] = P_count_tissue[tissue][t]/iteration
	return Pval_tissue, Fold_tissue, Mean_tissue, Top_assoc_pair

# main
if __name__ == '__main__':
	# exec main function with multiple cores
	pool = Pool(num_threads)
	results = pool.map(migwas_func, Algo_list)
	# wrap up the results from each algorithm
	Top_assoc_pair_dic = {}
	P_t = {}
	F_t = {}
	M_t = {}
	all_pval_t = {}
	all_fold_t = {}
	all_z_t = {}
	sumP_t = {}
	sumF_t = {}
	for i in range(len(Algo_list)):
		P_t[Algo_list[i]], F_t[Algo_list[i]], M_t[Algo_list[i]],Top_assoc_pair_dic[Algo_list[i]] = results[i]
		all_pval_t[Algo_list[i]] ={}
		all_fold_t[Algo_list[i]] ={}
		all_z_t[Algo_list[i]] ={}
		sumP_t[Algo_list[i]] ={}
		sumF_t[Algo_list[i]] ={}
		for tissue in P_t[Algo_list[i]].keys():
			all_pval_t[Algo_list[i]][tissue] = []
			all_fold_t[Algo_list[i]][tissue] = []
			all_z_t[Algo_list[i]][tissue] = []
		for t in range(len(THRES_list[Algo_list[i]])):
			for tissue in P_t[Algo_list[i]].keys():
				if M_t[Algo_list[i]][tissue][t] >= thresflag:
					all_pval_t[Algo_list[i]][tissue].append(P_t[Algo_list[i]][tissue][t])
					all_fold_t[Algo_list[i]][tissue].append(F_t[Algo_list[i]][tissue][t])
		for tissue in P_t[Algo_list[i]].keys():
			for x in range(len(all_pval_t[Algo_list[i]][tissue])):
				if all_pval_t[Algo_list[i]][tissue][x] == 1:
					all_pval_t[Algo_list[i]][tissue][x] = iteration/(iteration+1)
				elif all_pval_t[Algo_list[i]][tissue][x] == 0:
					all_pval_t[Algo_list[i]][tissue][x] = 1/(iteration+1) # permutation_P should never be zero
				all_z_t[Algo_list[i]][tissue].append(norm.ppf(all_pval_t[Algo_list[i]][tissue][x]))
			sumP_t[Algo_list[i]][tissue] = norm.cdf(np.sum(all_z_t[Algo_list[i]][tissue])/len(all_z_t[Algo_list[i]][tissue]))
			all_fold_t[Algo_list[i]][tissue] = np.array(all_fold_t[Algo_list[i]][tissue]) + 0.000000001
			sumF_t[Algo_list[i]][tissue] = np.exp(np.mean(np.log(all_fold_t[Algo_list[i]][tissue])))
	# output
	OUT = open(os.path.join(BASEDIR,args.out+"_migwas_result.txt"), 'w')
	print("#tissue\tP_value\tFold_change", file = OUT)
	Total_P_t = {}
	Total_F_t = {}
	final_P_t = {}
	final_F_t = {}
	for tissue in P_t[Algo_list[0]].keys():
		Total_P_t[tissue] = []
		Total_F_t[tissue] = []
		for algorithm in Algo_list:
			Total_P_t[tissue].append(sumP_t[algorithm][tissue])
			Total_F_t[tissue].append(sumF_t[algorithm][tissue])
		Total_P_t[tissue] = np.array(Total_P_t[tissue])
		index_nan_p = np.where(Total_P_t[tissue] != Total_P_t[tissue]) # get nan
		Total_P_t[tissue] = np.delete(Total_P_t[tissue], index_nan_p, axis = 0)
		final_P_t[tissue] = norm.cdf(np.sum(norm.ppf(Total_P_t[tissue]))/(len(Total_P_t[tissue]))**0.5)
	
		Total_F_t[tissue] = np.array(Total_F_t[tissue])
		index_nan_f = np.where(Total_F_t[tissue] != Total_F_t[tissue]) # get nan
		Total_F_t[tissue] = np.delete(Total_F_t[tissue], index_nan_f, axis = 0)
		final_F_t[tissue] = np.exp(np.mean(np.log(Total_F_t[tissue])))
		print("{0}\t{1}\t{2}".format(tissue, final_P_t[tissue], final_F_t[tissue]), file = OUT)
	# top output (if ordered)
	if args.output_candidate:
		TOP_OUT = open(os.path.join(BASEDIR,args.out+"_candidates.txt"), 'w')
		len_thres = [len(THRES_list[algo]) for algo in Algo_list]
		for m in range(min(len_thres)):
			thres_percentile = percentile_list[m]
			pair_list_by_per = []
			for algo in Algo_list:
				for pair in Top_assoc_pair_dic[algo][m]:
					pair_list_by_per.append(pair)
			find_dup = set()
			dup_pair_list_by_per = [xx for xx in pair_list_by_per if xx in find_dup or find_dup.add(xx)]
			for out_pair in dup_pair_list_by_per:
				print("{0}\t{1}\t{2}\t{3}".format(thres_percentile,out_pair[0],out_pair[1],out_pair[2]),file = TOP_OUT)

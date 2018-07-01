# MIGWAS
This software evaluates enrichment of genome-wide association study (GWAS) signals on miRNA-target gene networks (MIGWAS) and partition them into various human tissues with the help of tissue specific miRNA expression data.

## Overview

## Requirements
- python 3.X
- scipy
- pandas
- six
- argparse
- math
- multiprocessing
- futures

## Installation
In order to get started with *MIGWAS*, you can just clone this repo as follows;
```bash
git clone https://github.com/saorisakaue/MIGWAS
cd ./MIGWAS
```

## Usage
### Step 0: Prepare your input
All you need is a text file with GWAS summary statistics.

| Column | Descriptions |
|:-----------:|:------------|
|1|rsID (*optional* with `--no-rsid` flag below.)|
|2|chromosome|
|3|BP position|
|4|z-score (*optional*)|
|5|p-value|

Please have a look at an example input at `./example/RA_trans.chr12.pos.P.txt` (RA GWAS result at chr12).

### Step 1: GWAS summary to gene- and miRNA- level *P* values
This part is based on the excelent work by Masahiro Kanai, which was implemented to calculate the corrected "gene association score" from a GWAS result, according to [MAGENTA](https://www.broadinstitute.org/mpg/magenta/)'s method. For detailed explanations, please visit [the original repository](https://github.com/mkanai/minimgnt).
The example command is as follows;
```bash
$ python3 ./minimgnt.py score_filename --out RA_trans [--cpus 4] [--not-remove-HLA] [--remove-NA] --no-rsid
```
#### Arguments and options
* *`score_filename`* : GWAS summary statistics.

| Option name | Descriptions | Required | Default |
|:-----------:|:------------|:------------:|:------------|
|`--out`, `-o`| An output prefix. This should preferably be a phenotype name. | Yes | "your_phenotype" |
| `--cpus`, `-j` | a number of cpus used for computation. | No | 1 |
| `--not-remove-HLA` | do not remove genes in HLA region from a result. | No | False |
| `--remove-NA` | remove genes with NA score from the output. | No | False |
| `--no-rsid` | use this flag when a score file doesn't contain a rsID column. | No | False　|

Output files will be generated at `./miRNA_P/` and `./gene_P/`.

### Step 2: MIGWAS analysis for all tissues and specific tissues

The example command is as follows;
```bash
$ python3 ./migwas.py --phenotype RA_trans --out miRA_RA [--cpus 4] --iterations 20000 [--output-candidate]
```

#### Options
| Option name | Descriptions | Required | Default |
|:-----------:|:------------|:------------:|:------------|
| `--phenotype`, `-p` | Name of the phenotype of interest (file name prefix from minimgnt output). | Yes | None |
| `--out`, `-o` | Output file prefix. | Yes | "your_migwas" |
| `--cpus`, `-j` | Number of CPUs to be used. | No | 1 |
| `--iterations`, `-i` | Number of iterations to simulate null distributions.| Yes | 20000 |
| `--output-candidate`, `-c` | If you want to output a list of candidate miRNAs and genes associated with the trait, set this flag.| No | False |

## Output
The example output is as follows;
```bash
$ head miRA_RA_migwas_result.txt
#tissue	P_value	Fold_change
endothelial_cell_of_hepatic_sinusoid	0.0734928970649366	1.3943485234788573
epithelial_cell_of_proximal_tubule	0.22659192222025776	0.9410164034752201
keratinocyte	0.3550971896625686	0.7721248194464695
```
Each cell's partitioned enrichment *P* value and fold change, as well as miRNA-gene enrichment for all tissues will be described.

## Acknowledgements
* Tissue specific miRNA-gene enrichment analysis was made possible by the awesome work from [FANTOM5](http://fantom.gsc.riken.jp/5/), a comprehensive expression catalog of miRNA expression in varitous human cells. The original data can be found [here](http://fantom.gsc.riken.jp/5/suppl/De_Rie_et_al_2017/vis_viewer/#/human#srna;miRNA;hsa-miR-6859-5p).
* The original [MAGENTA](https://www.broadinstitute.org/mpg/magenta/) was written by Ayellet Segre, Mark Daly, and David Altshuler of The Broad Institute.
    * Ayellet V. Segrè, DIAGRAM Consortium, MAGIC investigators, Leif Groop, Vamsi K. Mootha, Mark J. Daly, and David Altshuler (2010). **Common Inherited Variation in Mitochondrial Genes is not Enriched for Associations with Type 2 Diabetes or Related Glycemic Traits.** [PLoS Genetics Aug 12;6(8). pii: e1001058.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001058)
* Minimgnt (miniMAGENTA) part was written by [Masahiro Kanai](http://mkanai.github.io/), reimplementing the calculation of "gene associatino score" feature in Python.


## Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed.
Please refer to the [LICENCE](https://github.com/saorisakaue/MIGWAS/blob/master/LICENSE.md) page.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/3.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/3.0/88x31.png" /></a><br />

# Knockoff-based conditional differential gene expression (CDGE) analysis

This is a repository for R package performing knockoff-based CDGE analysis to identify direct effect genes. Here, we implement the CDGE analysis via the recently developed Ghostknockoff with lasso regression approach ([He et al., 2024](https://www.biorxiv.org/content/10.1101/2024.02.28.582621v3); [Chen et al., 2024](https://arxiv.org/abs/2402.12724)). 

## Dependences

### R packages on GitHub

- [R Package "knockoffsr"](https://github.com/biona001/knockoffsr)
- [R Package "ghostbasil"](https://github.com/JamesYang007/ghostbasil)

For the installation of the `ghostbasil` package, you can try the following commands.

```
git clone https://github.com/JamesYang007/ghostbasil
cd "…/ghostbasil/R"
ml mpfr
ml gcc
ml R
ml libgit2
ml system harfbuzz fribidi
git pull 
R
```

```R
install.packages('devtools') 
library(devtools)
install()
```

### R packages on CRAN

- [R Package "Matrix"](https://cran.r-project.org/web/packages/Matrix/index.html)
- [R Package "data.table"](https://cran.r-project.org/web/packages/data.table/index.html)
- [R Package "RSpectra"](https://cran.r-project.org/web/packages/RSpectra/index.html)
- [R Package "rsvd"](https://cran.r-project.org/web/packages/rsvd/index.html)
- [R Package "Rdsdp"](https://cran.r-project.org/web/packages/data.table/index.html)

## Installation

`KnockoffCDGE` is not on CRAN yet. You can install this package using devtools:
```R
library(devtools)
install_github("GuJQ5/KnockoffCDGE")
```
Then, you can load the package within R:
```R
library(KnockoffCDGE)
```

## Usage

`KnockoffCDGE` provides a one-stop inference of direct effect genes via the knockoff-based CDGE analysis. To begin with, you need inference results of existing differential gene expression (DGE) analysis and the expression level correlation matrix among genes of interest. In this package, a sample dataset of 682 genes in 23010 T-cells of 3 psoriasis patients is provided. This sample data is derived from the single-cell RNA-seq dataset collected by [Reynolds et al. (2021)](https://www.science.org/doi/10.1126/science.aba6500) in the study of inflammatory skin diseases.

```R
data("DGE_result")
data("Sigma")

head(DGE_result)
     Gene  baseMean log2FoldChange     lfcSE      pvalue       padj
1   BANF1 186.26196    -0.50524646 0.2282036 0.005705662 0.04675831
2  NDUFB4 626.66869     0.01660672 0.1266814 0.889103693 0.94995882
3 NDUFA11 421.51932    -0.40426853 0.2084349 0.019826898 0.10722996
4  NDUFV3  58.26491    -0.07503540 0.2138109 0.669875211 0.82787117
5   FOXO3  34.58603     0.17248214 0.2464750 0.358840484 0.60047479
6    HSF2  24.01166     0.48504321 0.2540897 0.011867007 0.07663629

Sigma[1:6,1:6]
             BANF1     NDUFB4    NDUFA11     NDUFV3      FOXO3       HSF2
BANF1   1.00000000 0.14262821 0.18673149 0.07655627 0.04272639 0.01529431
NDUFB4  0.14262821 1.00000000 0.19119376 0.08778896 0.04842344 0.02653903
NDUFA11 0.18673149 0.19119376 1.00000000 0.08088342 0.04469559 0.02334994
NDUFV3  0.07655627 0.08778896 0.08088342 1.00000000 0.01844828 0.01193820
FOXO3   0.04272639 0.04842344 0.04469559 0.01844828 1.00000000 0.02012618
HSF2    0.01529431 0.02653903 0.02334994 0.01193820 0.02012618 1.00000000
```

Then you should translate DGE analysis result into DGE Z-scores using p-values and log2FoldChange.

```R
Z<-Z_calculation(DGE_result$pvalue,DGE_result$log2FoldChange)
Z[1:10]
[1] -2.7642339  0.1394384 -2.3296076 -0.4263193  0.9175768  2.5160746  0.4755489  2.1628696 -1.8848544  2.8180229
```

Finally, you can perform knockoff-based CDGE analysis with DGE Z-scores, correlation matrix as follows.

```R
set.seed(433)
CDGE_result<-KnockoffCDGE(Z,Sigma,M=5,n=23010,method="ME",Gene_info=DGE_result,verbose=TRUE,tol=0.001)

682 representatives for 682 variables, 682 optimization variables
Iter 1: δ = 0.35852169562751923
Iter 2: δ = 0.3912557343340361
Iter 3: δ = 0.016588363288380092
Iter 4: δ = 0.0004856336999177757

CDGE_result[1:20,]
       Gene   baseMean log2FoldChange     lfcSE       pvalue         padj          Z         S  q_value_ko
1     BANF1 186.261961    -0.50524646 0.2282036 5.705662e-03 4.675831e-02 -2.7642339 0.4583510         Inf
2    NDUFB4 626.668685     0.01660672 0.1266814 8.891037e-01 9.499588e-01  0.1394384 0.4298007         Inf
3   NDUFA11 421.519323    -0.40426853 0.2084349 1.982690e-02 1.072300e-01 -2.3296076 0.4124828         Inf
4    NDUFV3  58.264909    -0.07503540 0.2138109 6.698752e-01 8.278712e-01 -0.4263193 0.5383378         Inf
5     FOXO3  34.586035     0.17248214 0.2464750 3.588405e-01 6.004748e-01  0.9175768 0.5259153         Inf
6      HSF2  24.011662     0.48504321 0.2540897 1.186701e-02 7.663629e-02  2.5160746 0.5569751         Inf
7     PYGO2  41.733923     0.08348481 0.2128035 6.343958e-01 8.056436e-01  0.4755489 0.5438025         Inf
8   TBL1XR1  93.597518     0.36966907 0.2032614 3.055121e-02 1.408559e-01  2.1628696 0.5243315         Inf
9    GTF2H5  61.264928    -0.37910451 0.3090825 5.944949e-02 2.145234e-01 -1.8848544 0.5066855         Inf
10   ADAM19  48.140643     0.62940215 0.4078989 4.832036e-03 4.258618e-02  2.8180229 0.5146926         Inf
11   ALKBH1  17.375140     0.17795360 0.2540359 3.486060e-01 5.899006e-01  0.9372966 0.5469624         Inf
12  METTL16  56.079196     0.61844229 0.2131691 4.184633e-04 7.767141e-03  3.5281574 0.5415033 0.001724138
13    INTS5  20.143503    -0.31148375 0.2670612 1.085424e-01 3.049661e-01 -1.6047794 0.5409285         Inf
14     GLCE  16.644493    -0.42156159 0.2580063 2.926158e-02 1.368207e-01 -2.1799444 0.5516036         Inf
15    SNX30   4.975468    -0.10581563 0.3258309 5.673302e-01 7.627043e-01 -0.5719879 0.5385734         Inf
16  PRPSAP2  49.154502    -0.12288147 0.1750012 4.267907e-01 6.609708e-01 -0.7946954 0.5307836         Inf
17    UBE4B  44.417665    -0.01907139 0.2119424 9.132875e-01 9.614455e-01 -0.1088928 0.5496289         Inf
18   RPRD1B  25.630895     0.11790345 0.2278416 5.167191e-01 7.286133e-01  0.6484111 0.5485474         Inf
19   PCSK1N  74.377222    -1.30504254 0.2687645 2.386519e-09 4.608676e-07 -5.9690412 0.5079039 0.001724138
20 FRA10AC1  32.771177    -0.25062612 0.2547872 1.886541e-01 4.224341e-01 -1.3145711 0.5411116         Inf
```

Here, argument `M` characterizes the number of knockoff copies you want to generate under the multiple knockoff framework, `n`  depicts the sample size of the analysis, `method` specifies the criterion to generate knockoffs, `Gene_info` includes supplementary information of genes of interest (if available), `verbose` indicates whether progress of computation is plotted and `tol` specifies the convergence tolerance in computing S-matrix for knockoff generation. The current version of this package provides four options for the knockoff criterion.

- `"ME"`: **Maximizing the entropy criterion** via interface to Julia ([Gimenez and Zou, 2019](https://proceedings.mlr.press/v89/gimenez19b.html); [Chu et al., 2024](https://doi.org/10.1093/bioinformatics/btae580));
- `"SDP"`: **Semidefinite program criterion** via interface to Julia ([Barber and Candès, 2015](https://projecteuclid.org/journals/annals-of-statistics/volume-43/issue-5/Controlling-the-false-discovery-rate-via-knockoffs/10.1214/15-AOS1337.full); [Candès et al., 2018](https://academic.oup.com/jrsssb/article/80/3/551/7048447));
- `"SDP_no_julia"`: **Semidefinite program criterion** without the need of interface to Julia.
- `"PCA"`: **PCA knockoffs** ([Fan et al., 2020](https://doi.org/10.1080/01621459.2019.1654878)).

Specifically, you can use the latter two if the installation of Julia fails.

For FDR-control inference, you can use knockoff q-values of all features in the column `"q_value_ko"`. Here, the q-value characterizes the smallest target FDR level for each feature to be selected.

## Citation

If you use this software in a research paper, please cite our paper.

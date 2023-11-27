# MAAT
## Workflow of MAAT

![image](https://github.com/wanghanmath/MAAT/assets/33062920/4679523d-a688-4fbb-b3d3-fa9bfc4f261b)

## Data Preparation

We need to load a genotype matrix, a gene expression vector and an annotation information matrix. A pre-specified sparsity level for the effect size vector is also needed.

The SNPs in the annotation information matrix are required to match SNPs in the genotype matrix.

## Gene Expression Imputation by incorporating multiple annotations
`maat_effect_size.R` is the main function of MAAT, which utilizes a product partition model with covariates (PPMx) to incorporate multifaceted annotation information into the TWAS imputation step.
We need to specify four arguments in `maat_effect_size.R`, i.e., the working directory, the name of genotype matrix, the name of gene expression vector, the name of annotation information matrix, and a sparsity level.

### Example Run
We can train a per-gene expression imputation model using the following code

```Rscript maat_effect_size.R my-folder genotype.csv expression.csv annotation.csv 0.1```
## Association Test
We use the test statistics in FUSION to perform gene-based association test
$$Z=\frac{Z_{SG}^{\top}\beta}{\beta^{\top} V\beta}$$
here $\beta$ is the effect size estimated from the first imputation step, $Z_{SG}$ is the $Z$-score from GWAS summary-level statistics, $V$ is the the LD matrix of analyzed SNPs.
## Annotation Assignment
We use a metric based on cosine distance to define the important annotation label for each significant gene-trait association
![image](https://github.com/wanghanmath/MAAT/assets/33062920/3c378c39-51cf-49ab-91ef-14bdbb82c5d4)


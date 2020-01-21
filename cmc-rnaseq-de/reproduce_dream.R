devtools::install_github("GabrielHoffman/variancePartition@3287d7b")
library(variancePartition)
library(tibble)
library(synapser)
synLogin()
param = BiocParallel::SnowParam(4, "SOCK", progressbar=TRUE)
BiocParallel::register(param)

COVARIATES <- readr::read_csv(synGet("syn21488871")$path) %>% 
  tibble::column_to_rownames(var = "SampleID")

NEW.COUNTS <- readRDS(synGet("syn21535143")$path)

CQN.GENE_EXPRESSION <- readRDS(synGet("syn21535146")$path)

# Post adjusted formula
random_effect <- "Individual_ID"
adjust.covars <- c("Dx.Tissue", "RIN")
formula <- glue::glue("~ ", glue::glue_collapse(adjust.covars, sep = " + ")) %>% 
  glue::glue(" + (1|{random_effect})")

tmp = NEW.COUNTS
tmp[is.na(tmp)] = 0
geneExpr = DGEList(tmp)
geneExpr = edgeR::calcNormFactors(geneExpr)

VOOM.GENE_EXPRESSION = variancePartition::voomWithDreamWeights(counts = geneExpr, 
                                                               formula = formula,
                                                               data = COVARIATES)

# Fit linear model using new weights and new design
# CQN NA counts to 0 counts for input to dream()
tmp = CQN.GENE_EXPRESSION$E
tmp[is.na(tmp)] = 0
VOOM.GENE_EXPRESSION$E = tmp
ADJUSTED.FIT = variancePartition::dream(exprObj = VOOM.GENE_EXPRESSION, 
                                        formula = formula,
                                        data = COVARIATES,
                                        computeResiduals = TRUE)

RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(ADJUSTED.FIT, CQN.GENE_EXPRESSION$E)




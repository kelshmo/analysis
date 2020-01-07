library("tibble")
synLogin()
param = BiocParallel::SnowParam(4, "SOCK", progressbar=TRUE)
BiocParallel::register(param)

metadata <- readr::read_csv(synGet("syn21488871")$path) %>% 
  tibble::column_to_rownames(var = "SampleID")

counts <- readr::read_csv(synGet("syn21488882")$path) %>% 
  tibble::column_to_rownames(var = "geneID") %>% 
  as.matrix()

tmp = counts
tmp[is.na(tmp)] = 0

dream_formula <- ~ Dx.Tissue + (1|Individual_ID)

vobjDream <- variancePartition::voomWithDreamWeights(counts = tmp, formula = dream_formula, data = metadata)



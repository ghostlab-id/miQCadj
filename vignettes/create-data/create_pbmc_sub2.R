# create dataset for miQCadj vignette
library(Seurat)
library(SeuratData)
library(SeuratObject)
pbmcsca <- SeuratData::LoadData("pbmcsca")
pbmc <- SeuratObject::UpdateSeuratObject(pbmcsca)
pbmcsca <- NULL
Seurat::DefaultAssay(pbmc) <- "RNA"

# adjust meta data column names
# Method -> tech
colnames(pbmc@meta.data)[colnames(pbmc@meta.data)=="Method"] <- "tech"
# remove "-" in tech entries
pbmc@meta.data$tech <- gsub("\\-", "", pbmc@meta.data$tech)
# Experiment -> sample
colnames(pbmc@meta.data)[colnames(pbmc@meta.data)=="Experiment"] <- "sample"
# combine sample & tech
pbmc[["sample_tech"]] <- paste(pbmc$sample, pbmc$tech, sep="_")

# remove duplicated columns
pbmc@meta.data[c("nUMI", "nGene")] <- NULL

# convert percent.mito into numeric
pbmc@meta.data$percent.mito <- as.numeric(pbmc@meta.data$percent.mito)

# select sample_tech objects for vignette
pbmc_sub2 <- subset(pbmc, subset=sample_tech %in% c("pbmc1_Smartseq2", "pbmc2_CELSeq2"))

# The stored percent.mito shows small deviations (~0.3%), so it is recomputed for consistency.
colnames(pbmc_sub2)[which(colnames(pbmc_sub2)=="percent.mito")] <- "percent.mt"

## calculate percent feature using seurat
pbmc_sub2[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")/100

counts <- pbmc_sub2@assays$RNA@layers$counts
meta <- pbmc_sub2@meta.data
pbmc_sub2 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)

utils::str(pbmc_sub2)
# Formal class 'Seurat' [package "SeuratObject"] with 13 slots
# ..@ assays      :List of 1
# .. ..$ RNA:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
# .. .. .. ..@ layers    :List of 1
# .. .. .. .. ..$ counts:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# .. .. .. .. .. .. ..@ i       : int [1:1313914] 2 5 7 11 12 18 21 22 38 44 ...
# .. .. .. .. .. .. ..@ p       : int [1:527] 0 2200 4638 6512 8992 11188 13404 15735 17810 20640 ...
# .. .. .. .. .. .. ..@ Dim     : int [1:2] 33694 526
# .. .. .. .. .. .. ..@ Dimnames:List of 2
# .. .. .. .. .. .. .. ..$ : NULL
# .. .. .. .. .. .. .. ..$ : NULL
# .. .. .. .. .. .. ..@ x       : num [1:1313914] 366 1 35 1 388 6 441 154 232 375 ...
# .. .. .. .. .. .. ..@ factors : list()
# .. .. .. ..@ cells     :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
# .. .. .. .. .. ..@ .Data: logi [1:526, 1] TRUE TRUE TRUE TRUE TRUE TRUE ...
# .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. .. .. .. .. ..$ : chr [1:526] "Cell_1" "Cell_2" "Cell_3" "Cell_4" ...
# .. .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. .. .. ..$ dim     : int [1:2] 526 1
# .. .. .. .. .. ..$ dimnames:List of 2
# .. .. .. .. .. .. ..$ : chr [1:526] "Cell_1" "Cell_2" "Cell_3" "Cell_4" ...
# .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. ..@ features  :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
# .. .. .. .. .. ..@ .Data: logi [1:33694, 1] TRUE TRUE TRUE TRUE TRUE TRUE ...
# .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. .. .. .. .. ..$ : chr [1:33694] "Feature1" "Feature2" "Feature3" "Feature4" ...
# .. .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. .. .. ..$ dim     : int [1:2] 33694 1
# .. .. .. .. .. ..$ dimnames:List of 2
# .. .. .. .. .. .. ..$ : chr [1:33694] "Feature1" "Feature2" "Feature3" "Feature4" ...
# .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. ..@ default   : int 1
# .. .. .. ..@ assay.orig: chr(0) 
# .. .. .. ..@ meta.data :'data.frame':	33694 obs. of  0 variables
# .. .. .. ..@ misc      :List of 1
# .. .. .. .. ..$ calcN: logi TRUE
# .. .. .. ..@ key       : chr "rna_"
# ..@ meta.data   :'data.frame':	526 obs. of  10 variables:
#   .. ..$ orig.ident  : Factor w/ 2 levels "pbmc1","pbmc2": 1 1 1 1 1 1 1 1 1 1 ...
# .. ..$ nCount_RNA  : num [1:526] 437125 335596 302204 377420 385514 ...
# .. ..$ nFeature_RNA: int [1:526] 2200 2438 1874 2480 2196 2216 2331 2075 2830 2263 ...
# .. ..$ percent.mito: num [1:526] 0.0297 0.0312 0.0431 0.026 0.0405 ...
# .. ..$ Cluster     : chr [1:526] "0" "0" "0" "0" ...
# .. ..$ CellType    : chr [1:526] "Cytotoxic T cell" "Cytotoxic T cell" "Cytotoxic T cell" "Cytotoxic T cell" ...
# .. ..$ sample      : chr [1:526] "pbmc1" "pbmc1" "pbmc1" "pbmc1" ...
# .. ..$ tech        : chr [1:526] "Smartseq2" "Smartseq2" "Smartseq2" "Smartseq2" ...
# .. ..$ sample_tech : chr [1:526] "pbmc1_Smartseq2" "pbmc1_Smartseq2" "pbmc1_Smartseq2" "pbmc1_Smartseq2" ...
# .. ..$ percent.mt  : num [1:526] 0.0425 0.0396 0.0639 0.0344 0.0516 ...
# ..@ active.assay: chr "RNA"
# ..@ active.ident: Factor w/ 1 level "Cell": 1 1 1 1 1 1 1 1 1 1 ...
# .. ..- attr(*, "names")= chr [1:526] "Cell_1" "Cell_2" "Cell_3" "Cell_4" ...
# ..@ graphs      : list()
# ..@ neighbors   : list()
# ..@ reductions  : list()
# ..@ images      : list()
# ..@ project.name: chr "SeuratProject"
# ..@ misc        : list()
# ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
# .. ..$ : int [1:3] 5 0 2
# ..@ commands    : list()
# ..@ tools       : list()

utils::object.size(pbmc_sub2)
# 18497096 bytes

# Rename the final object
usethis::use_data(pbmc_sub2, compress = "xz", overwrite = TRUE)
# ✔ Setting active project to 'C:/Users/Eli/Desktop/Patho/temp_git/miQCadj'
# ✔ Setting LazyData to 'true' in 'DESCRIPTION'
# ✔ Saving 'pbmc_sub2' to 'data/pbmc_sub2.rda'
# • Document your data (see 'https://r-pkgs.org/data.html')

# cleanup
rm(pbmc, pbmc_sub2)

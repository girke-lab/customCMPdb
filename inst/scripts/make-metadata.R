library(AnnotationHub)
library(AnnotationHubData)
meta <- data.frame(
  Title = c("annot_0.1", 
            "drugage_build2", "drugbank_5.1.5", "cmap02", "lincs_pilot1"),
  Description = c("SQLite database containing compound annotations from four resources",
                  "SDF file of compounds in DrugAge database",
                  "SDF file of compounds in DrugBank database",
                  "SDF file of compounds in CMAP02 database",
                  "SDF file of compounds in LINCS database"),
  BiocVersion = rep("3.11", 5),
  Genome = rep("GRCh38", 5),
  # getValidSourceTypes()
  SourceType = c("TSV", "XLS/XLSX", "TXT",
                 "TSV", "TSV"),
  SourceUrl = c("https://bit.ly/3dCWKWo, https://github.com/yduan004/drugbankR, https://github.com/tgirke/longevityDrugs, ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
                "https://bit.ly/3dCWKWo",
                "https://www.drugbank.ca/releases/latest#structures",
                "https://github.com/tgirke/longevityDrugs",
                "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz"),
  SourceVersion = c("Build 2; 5.1.5; Build02; Pilot Phase 1", 
                    "Build 2", "5.1.5",
                    "Build02", "Pilot Phase 1"),
  Species = rep("Homo sapiens", 5),
  TaxonomyId = rep(9606, 5),
  Coordinate_1_based = rep(TRUE, 5),
  DataProvider = c("DrugAge, DrugBank, Broad Institute, Broad Institute",
                   "DrugAge", "DrugBank",
                   "Broad Institute", "Broad Institute"),
  Maintainer = rep("Yuzhu Duan <yduan004@ucr.edu>", 5),
  RDataClass = c("character", 
                 "character", "character",
                 "character", "character"),
  # AnnotationHub::DispatchClassList()
  DispatchClass = c("FilePath", 
                    "FilePath", "FilePath",
                    "FilePath", "FilePath"),
  RDataPath = c("compoundCollectionData/compoundCollection_0.1.db",
                "compoundCollectionData/drugage_build2.sdf",
                "compoundCollectionData/drugbank_5.1.5.sdf",
                "compoundCollectionData/cmap02.sdf",
                "compoundCollectionData/lincs_pilot1.sdf"),
  Tags = c("annot", "drugage", "drugbank",
           "cmap", "lincs")
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)

AnnotationHubData::makeAnnotationHubMetadata("../compoundCollectionData")

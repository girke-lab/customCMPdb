#' Compound Annotation and Structure Datasets
#'
#' @name compoundCollectionData-package
#' @aliases compoundCollectionData-package compoundCollectionData
#' @docType package
#' @description This package contains annotation and 
#' structure datasets for compounds in DrugAge, DrugBank, CMAP02 and LINCS 
#' databases. The actual datasets are stored in 
#' \code{\link[AnnotationHub]{AnnotationHub}}
#' 
#' @details
#' The description of the 5 datasets in this package is as follows:
#' 
#' \strong{Annotation SQLite database:}
#' 
#' It is a SQLite database storing compound annotation tables for DrugAge, 
#' DrugBank, CMAP02 and LINCS, respectively. It also contains an ID mapping 
#' table of ChEMBL ID to IDs of individual databases. The annotation database
#' was stored in the 'compoundCollection_0.1.db' file
#' 
#' \strong{DrugAge SDF:}
#' 
#' It is an SDF (Structure-Data File) file storing molecular structures of 
#' DrugAge comopunds. The source DrugAge annotation file was downloaded from
#' \href{http://genomics.senescence.info/drugs/dataset.zip}{here}. The extracted csv 
#' file only contains drug names, without id mappings to external resources 
#' such as PubChem or ChEMBL. This 'drugage.csv' file was further processed by the 
#' \code{processDrugage} function in this package. The result DrugAge annotation table
#' as well as the id-mapping table (DrugAge internal id to ChEMBL ID) were then
#' stored in an SQLite database named as 'compoundCollection'. 
#' The drug structures were obtained from PubChem CIDs by \code{getIds} function 
#' from \pkg{ChemmineR} package. The 'SDFset' object was then written to the 
#' 'drugage_build2.sdf' file
#' 
#' \strong{DrugBank SDF:}
#' 
#' This SDF file stores structures of compounds in 
#' \href{https://www.drugbank.ca/}{DrugBank} database. The full DrugBank xml 
#' file was downloaded from \url{https://www.drugbank.ca/releases/latest}.
#' The most recent release version at the time of writing this document is 5.1.5.  
#' The extracted xml file was procesed by the \code{dbxml2df} function in this package.
#' The result DrugBank annotation table was then stored in the compoundCollection 
#' SQLite database. The DrugBank to ChEMBL id mappings were obtained from 
#' \href{ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz}{UniChem}.
#' The DrugBank SDF file was downloaded from 
#' \url{https://www.drugbank.ca/releases/latest#structures}
#' and made some validity checks and modifications via utilities in the 
#' \pkg{ChemmineR} package. The results were written to the 'drugbank_5.1.5.sdf' file
#' 
#' \strong{CMAP SDF:}
#' 
#' The CMAP compound instance table was downloaded from 
#' \href{http://www.broadinstitute.org/cmap/cmap_instances_02.xls}{CMAP02} 
#' website and processed by the \code{buildCMAPdb} function
#' in this package. The result 'cmap.db' contains both compound annotation and 
#' structure information. 
#' Since the annotation table only contains PubChem CID, the ChEMBL ids were added 
#' via PubChem CID to ChEMBL id mappings from 
#' \href{ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz}{UniChem}.
#' The CMAP internal IDs were made for ChEMBL id to CMAP id mappings. The 
#' structures was written to the 'cmap02.sdf' file
#' 
#' \strong{LINCS SDF:}
#' 
#' The LINCS compound annotation table was downloaded from
#' \href{ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz}{GEO}.
#' where only compounds type were selected.
#' The LINCS ids were mapped to ChEMBL ids via inchi key. The LINCS compounds 
#' structures were obtained from PubChem CIDs via \code{getIDs} function from
#' \pkg{ChemmineR} package. The structures were written to the 'lincs_pilot1.sdf' file
#' 
#' The R script of generating the above 5 datasets is available at the 
#' 'inst/scripts/make-data.R' file in this package.
#' 
#' @import AnnotationHub
#' @author
#' \itemize{
#'   \item Yuzhu Duan (yduan004@ucr.edu)
#'   \item Thomas Girke (thomas.girke@ucr.edu)
#' } 
#' @examples 
#' library(AnnotationHub)
#' \dontrun{
#'     ah <- AnnotationHub()
#'     
#'     ## Load compoundCollection annotation SQLite database
#'     query(ah, c("compoundCollectionData", "annot_0.1"))
#'     annot_path <- ah[["AH79563"]]
#'     library(RSQLite)
#'     conn <- dbConnect(SQLite(), annot_path)
#'     dbListTables(conn)
#'     dbDisconnect(conn)
#'     
#'     ## Load DrugAge SDF file
#'     query(ah, c("compoundCollectionData", "drugage_build2"))
#'     da_path <- ah[["AH79564"]]
#'     da_sdfset <- ChemmineR::read.SDFset(da_path)
#'     
#'     ## Load DrugBank SDF file
#'     query(ah, c("compoundCollectionData", "drugbank_5.1.5"))
#'     db_path <- ah[["AH79565"]]
#'     db_sdfset <- ChemmineR::read.SDFset(db_path)
#'     
#'     ## Load CMAP SDF file
#'     query(ah, c("compoundCollectionData", "cmap02"))
#'     cmap_path <- ah[["AH79566"]]
#'     cmap_sdfset <- ChemmineR::read.SDFset(cmap_path)
#'     
#'     ## Load LINCS SDF file
#'     query(ah, c("compoundCollectionData", "lincs_pilot1"))
#'     lincs_path <- ah[["AH79567"]]
#'     lincs_sdfset <- ChemmineR::read.SDFset(lincs_path)
#' }
#' 

NULL
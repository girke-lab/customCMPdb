#' Customize and Query Compound Annotation Database
#'
#' @name customCMPdb-package
#' @aliases customCMPdb-package customCMPdb
#' @docType package
#' @description 
#' This package is served as the query and customization interface for compound 
#' annotations from \href{https://genomics.senescence.info/drugs/}{DrugAge}, 
#' \href{https://www.drugbank.ca/}{DrugBank}, 
#' \href{https://portals.broadinstitute.org/cmap/}{CMAP02}
#' and \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742}{LINCS} databases.
#' It also stores the structure SDF datasets for compounds in the above four databases.
#' 
#' Specifically, the annotation database created by this package is an SQLite database
#' containing 5 tables, including 4 compound annotation tables from DrugAge, 
#' DrugBank, CMAP02 and LINCS databases, respectively. The other one is an ID 
#' mapping table of ChEMBL IDs to IDs of individual databases. The other 4 datasets 
#' stores the structures of compounds in the DrugAge, DrugBank, CMAP02 and LINCS 
#' databases in SDF files. For detailed description of the 5 datasets generated
#' by this package, please consult to the vignette of this package by running
#' \code{browseVignettes("customCMPdb")} The actual datasets are hosted in 
#' \code{\link[AnnotationHub]{AnnotationHub}}.
#' 
#' This package also provides functionalities to customize and query the compound
#' annotation SQLite database. Users could add their customized compound annotation
#' tables to the SQLite database and query both the default (DrugAge, DrugBank, CMAP02,
#' LINCS) and customized annotations by providing ChEMBL ids of the query compounds.
#' The customization and query functions are available at \code{\link{customAnnot}}
#' and \code{\link{queryAnnotDB}}, respectively. 
#'  
#' @details
#' The description of the 5 datasets in this package is as follows.
#' 
#' \strong{Annotation SQLite database:}
#' 
#' It is a SQLite database storing compound annotation tables for DrugAge, 
#' DrugBank, CMAP02 and LINCS, respectively. It also contains an ID mapping 
#' table of ChEMBL ID to IDs of individual databases. 
#' 
#' \strong{DrugAge SDF:}
#' 
#' It is an SDF (Structure-Data File) file storing molecular structures of 
#' DrugAge compounds. The source DrugAge annotation file was downloaded from
#' \href{http://genomics.senescence.info/drugs/dataset.zip}{here}. The extracted csv 
#' file only contains drug names, without id mappings to external resources 
#' such as PubChem or ChEMBL. The extracted 'drugage.csv' file was further processed by the 
#' \code{\link{processDrugage}} function in this package. The result DrugAge annotation table
#' as well as the id-mapping table (DrugAge internal id to ChEMBL ID) were then
#' stored in the SQLite annotation database named as 'compoundCollection'. 
#' The drug structures were obtained from PubChem CIDs by \code{\link[ChemmineR]{getIds}} 
#' function from \pkg{ChemmineR} package. The \code{SDFset} object was then 
#' written to the \code{drugage_build2.sdf} file
#' 
#' \strong{DrugBank SDF:}
#' 
#' This SDF file stores structures of compounds in 
#' \href{https://www.drugbank.ca/}{DrugBank} database. The full DrugBank xml 
#' file was downloaded from \url{https://www.drugbank.ca/releases/latest}.
#' The most recent release version at the time of writing this document is 5.1.5.  
#' The extracted xml file was processed by the \code{\link{dbxml2df}} function in this package.
#' The result DrugBank annotation table was then stored in the \code{compoundCollection}
#' SQLite database. The DrugBank to ChEMBL id mappings were obtained from 
#' \href{ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz}{UniChem}.
#' The DrugBank SDF file was downloaded from 
#' \url{https://www.drugbank.ca/releases/latest#structures}.
#' Some validity checks and modifications were made via utilities in the 
#' \pkg{ChemmineR} package. The results were written to the \code{drugbank_5.1.5.sdf} file
#' 
#' \strong{CMAP SDF:}
#' 
#' The CMAP compound instance table was downloaded from 
#' \href{http://www.broadinstitute.org/cmap/cmap_instances_02.xls}{CMAP02} 
#' website and processed by the \code{\link{buildCMAPdb}} function
#' in this package. The result 'cmap.db' contains both compound annotation and 
#' structure information. 
#' Since the annotation table only contains PubChem CID, the ChEMBL ids were added 
#' via PubChem CID to ChEMBL id mappings from 
#' \href{ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz}{UniChem}.
#' The CMAP internal IDs were made for ChEMBL id to CMAP id mappings. The 
#' structures were written to the \code{cmap02.sdf} file
#' 
#' \strong{LINCS SDF:}
#' 
#' The LINCS compound annotation table was downloaded from
#' \href{ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz}{GEO}.
#' where only compounds type were selected.
#' The LINCS ids were mapped to ChEMBL ids via inchi key. The LINCS compounds 
#' structures were obtained from PubChem CIDs via \code{\link[ChemmineR]{getIds}} function from
#' \pkg{ChemmineR} package. The structures were written to the \code{lincs_pilot1.sdf} file
#' 
#' The R script of generating the above 5 datasets is available at the 
#' 'inst/scripts/make-data.R' file in this package.  The file location can
#' be found by running \code{system.file("scripts/make-data.R",package="compoundCollectionData")}
#' in user's R session or from the 
#' \href{https://github.com/yduan004/customCMPdb/blob/master/inst/scripts/make-data.R}{GitHub repository}
#' of this package.
#' 
#' @import AnnotationHub
#' 
#' @seealso
#' \code{\link{customAnnot}}, \code{\link{queryAnnotDB}}
#' 
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
#'     drugAgeAnnot <- dbReadTable(conn, "drugAgeAnnot")
#'     head(drugAgeAnnot)
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

NULL
#' Build DrugAge Annotation Database
#' 
#' This function builds the DrugAge annotation SQLite database from the 
#' 'drugage_id_mapping' table stored in the 'inst/extdata' directory of this 
#' package. The 'drugage_id_mapping.tsv' table contains the DrugAge compounds 
#' annotation information (such as species, avg_lifespan_change etc) as well 
#' as the compound name to ChEMBL id and PubChem id mappings. 
#' 
#' Part of the id mappings in the 'drugage_id_mapping.tsv' table is generated 
#' by the \code{processDrugage} function for compound names that have ChEMBL 
#' ids from the ChEMBL database (version 24). The missing IDs were added 
#' manually. A semi-manual approach was to use this web service: 
#' \url{https://cts.fiehnlab.ucdavis.edu/batch}. After the semi-manual process,
#' the left ones were manually mapped to ChEMBL, PubChem and DrugBank ids. 
#' The mixed items were commented. 
#' @param dest_path character(1), destination path of the result DrugAge 
#' annotation SQLite database
#' @return DrugAge annotation SQLite database 
#' @examples 
#'  buildDrugAgeDB(dest_path=tempfile(fileext="_drugage.db"))
#' @export
#'  
buildDrugAgeDB <- function(dest_path){
    ## Read DrugAge Mapping from inst/extdata
    
    # require(gsheet)
    # url<-"https://bit.ly/3dCWKWo"
    # drugAge_mapping <- suppressWarnings(gsheet2tbl(url))
    # ## Add columns names
    # names(drugAge_mapping) <- as.character(unlist(drugAge_mapping[1,]))
    # drugAge_mapping <- drugAge_mapping[-1,]
    
    da_path <- system.file("extdata/drugage_id_mapping.tsv", 
                           package="compoundCollectionData")
    drugAge_mapping <- read.delim(da_path)
    ## Create internal DrugAge_id named ida000xxx
    drugAge_ids <- paste0("ida",sprintf("%05d",1:nrow(drugAge_mapping)))
    drugAge_mapping2 <- data.frame(drugage_id=drugAge_ids, drugAge_mapping)
    
    ## Create DrugAge Database
    da <- dbConnect(SQLite(), dest_path)
    ## write id_mapping table.
    id_maptb <- na.omit(data.frame(drugage_id=drugAge_mapping2$drugage_id,
                                   chembl_id=drugAge_mapping2$chembl_id, 
                                   stringsAsFactors=FALSE)) 
    dbWriteTable(da,'id_mapping', id_maptb) 
    
    ## write DrugAge Annotation table
    drugAgeAnnot <- data.frame(drugAge_mapping2[1:13], drugAge_mapping2[15:16])
    dbWriteTable(da, 'drugAgeAnnot', drugAgeAnnot)
    dbDisconnect(da)
}

#' Process Source DrugAge Dataset
#' 
#' This function processes the source DrugAge datasets by adding the ChEMBL, 
#' PubChem and DrugBank id mapping information to the source DrugAge table 
#' which only has comopund names without id mapping information. Source file of 
#' DrugAge is linked here: \url{http://genomics.senescence.info/drugs/dataset.zip}
#' 
#' This function only annotates compound names that have ChEMBL 
#' ids from the ChEMBL database (version 24). The missing IDs were added 
#' manually. A semi-manual approach was to use this web service: 
#' \url{https://cts.fiehnlab.ucdavis.edu/batch}. After the semi-manual process,
#' the left ones were manually mapped to ChEMBL, PubChem and DrugBank ids. 
#' The mixed items were commented. 
#' @param drugagefile character(1), path to the destination DrugAge annotation 
#' file with id mappings
#' @param redownloaddrugage TRUE or FALSE indicating whether to redownload the 
#' DrugAge dataset if it has been downloaded before
#' @return write the 'drugage_id_mapping.tsv' table
#' @examples 
#' library(ChemmineR)
#' \dontrun{
#' processDrugage(drugagefile="drugage_id_mapping.tsv", redownloaddrugage=FALSE)
#' # Now the missing IDs need to be added manually. A semi-manual approach is to 
#' # use this web service: https://cts.fiehnlab.ucdavis.edu/batch
#' }
#' @export
#' 
processDrugage <- function(drugagefile="drugage_id_mapping.tsv", redownloaddrugage=FALSE) {
    if(redownloaddrugage==TRUE) {
        download.file("http://genomics.senescence.info/drugs/dataset.zip", "downloads/dataset.zip")
        unzip("downloads/dataset.zip", exdir="downloads/")
    }
    drugage <- read.csv("downloads/drugage.csv")
    drugage <- drugage[!duplicated(drugage$compound_name),]
    ## Get mol_dict table from chembl_db
    chembldb <- "/bigdata/girkelab/shared/lcshared/chemoinformatics/compoundDBs/chembl_24/chembl_24_sqlite/chembl_24.db"
    mydb <- dbConnect(SQLite(), chembldb)
    mol_dict <- dbGetQuery(mydb, 'SELECT * FROM molecule_dictionary')
    mol_dict_sub <- mol_dict[!is.na(mol_dict$pref_name),] # mol_dict from below
    prefname <- as.character(mol_dict_sub$pref_name)
    chemblid <- as.character(mol_dict_sub$chembl_id)
    fact <- tapply(chemblid, factor(prefname), paste, collapse=", ")
    ## Load ChEMBL to PubChem CID and DrugBank ID mappings (generated with downloadUniChem Fct)
    chembl2pubchem <- read.delim(gzfile("downloads/src1src22.txt.gz"))
    chembl2pubchem_vec <- tapply(as.character(chembl2pubchem[,2]), 
                                 factor(chembl2pubchem[,1]), paste, collapse=", ")
    chembl2drugbank <- read.delim(gzfile("downloads/src1src2.txt.gz"))
    chembl2drugbank_vec <- tapply(as.character(chembl2drugbank[,2]), 
                                  factor(chembl2drugbank[,1]), paste, collapse=", ")
    ## Assemble results
    drugage <- cbind(drugage, pref_name=names(fact[toupper(drugage$compound_name)]), 
                     chembl_id=fact[toupper(drugage$compound_name)])
    drugage <- cbind(drugage, pubchem_cid=chembl2pubchem_vec[as.character(drugage$chembl_id)], 
                     drugbank_id=chembl2drugbank_vec[as.character(drugage$chembl_id)])
    write.table(drugage, drugagefile, row.names=FALSE, quote=FALSE, sep="\t")
}

#############
## UniChem ##
#############
## UniChem CMP ID mappings from here: 
##      https://www.ebi.ac.uk/unichem/ucquery/listSources
## Note: above html table gives numbering to select proper src_id for ftp downloads, 
## e.g. DrugBank is src2
##      ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/
## Examples:
downloadUniChem <- function(rerun) {
    if(rerun==TRUE) {
        ## ChEMBL to DrugBank mapping in src1src2.txt.gz: 
        download.file(
            "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz", 
            "downloads/src1src2.txt.gz")
        ## ChEMBL to PubChem CID mapping in src1src22.txt.gz
        download.file(
            "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz", 
            "downloads/src1src22.txt.gz")
        ## ChEMBL to ChEBI mapping in src1src7.txt.gz
        download.file(
            "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src7.txt.gz", 
            "downloads/src1src7.txt.gz")
    }
}
## Usage:
# downloadUniChem(rerun=FALSE)

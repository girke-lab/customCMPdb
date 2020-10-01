#' Annotation Queries with Compound IDs
#'
#' This function can be used to query compound annotations from the default
#' resources as well as the custom resources stored in the SQLite annotation
#' database. The default annotation resources are DrugAge, DrugBank, CMAP02 and
#' LINCS. The customized compound annotations could be added/deleted by the
#' \code{\link{customAnnot}} utilities.
#' 
#' @details 
#' The input of this query function could be a set of ChEMBL IDs, it returns a
#' data.frame storing annotations of the input compounds from selected
#' annotation resources defined by the \code{annot} argument.
#' 
#' Since in the SQLite annotation database, ID identifiers from different ID systems, 
#' such as DrugBank and LINCS, are connected by ChEMBL IDs, it is hard to tell 
#' whether two IDs, such as DB00341, BRD-A42571354, refer to the same compound if 
#' either of them lack ID mappings to ChEMBL. So for querying compounds that don't 
#' have ChEMBL IDs, only one isolated database where the compounds belong to are 
#' supported. For example, a compound with LINCS id as "BRD-A00150179" doesn't have 
#' the ChEMBL ID mapping, when it is passed to the `chembl_id` argument, 
#' the `annot` need only to be set as `lincsAnnot` and the result
#' will be the compound annotation table from the LINCS annotation.
#'
#' @param chembl_id character vector of ChEMBL IDs or compound ids from other 
#' annotation system..
#' @param annot character vector of annotation resources, such as
#' \code{drugAgeAnnot}, \code{DrugBankAnnot}, \code{cmapAnnot}, \code{lincsAnnot} 
#' or names of the annotation tables added by users. The \code{\link{listAnnot}}
#' function lists the available options for the \code{annot} argument names.
#' @return data.frame of annotation result
#' @import RSQLite
#' @examples
#' query_id <- c("CHEMBL1000309", "CHEMBL100014", "CHEMBL100109",
#'                "CHEMBL100", "CHEMBL1000", "CHEMBL10")
#' qres <- queryAnnotDB(query_id, annot=c("drugAgeAnnot", "lincsAnnot"))
#'
#' # Add a custom compound annotation table
#' chembl_id <- c("CHEMBL1000309", "CHEMBL100014", "CHEMBL10",
#'                "CHEMBL100", "CHEMBL1000", NA)
#' annot_tb <- data.frame(compound_name=paste0("name", 1:6),
#'         chembl_id=chembl_id,
#'         feature1=paste0("f", 1:6),
#'         feature2=rnorm(6))
#' addCustomAnnot(annot_tb, annot_name="myCustom2")
#'
#' # query custom annotation
#' qres2 <- queryAnnotDB(query_id, annot=c("lincsAnnot", "myCustom2"))
#' 
#' # query compounds that don't have ChEMBL IDs
#' query_id <- c("BRD-A00474148", "BRD-A00150179", "BRD-A00763758", "BRD-A00267231")
#' qres3 <- queryAnnotDB(chembl_id=query_id, annot=c("lincsAnnot"))
#' qres3
#' @export
queryAnnotDB <- function(chembl_id,
        annot=c("drugAgeAnnot", "DrugBankAnnot", "cmapAnnot", "lincsAnnot")){
    ah <- AnnotationHub()
    dbpath <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), dbpath)
    
    da_cols <- c("species", "strain", "dosage", "avg_lifespan_change",
                 "max_lifespan_change", "gender", "significance")
    db_cols <- c("drugbank_id", "name", "cas.number", "unii", "state",
                 "groups", "fda.label")
    cmap_cols <- c("instance_id", "concentration..M.", "duration..h.", "cell2",
                   "array3", "perturbation_scan_id", "vehicle_scan_id4",
                   "scanner", "catalog_number")
    lincs_cols <- c("lincs_id", "pert_iname", "is_touchstone", "inchi_key",
                    "pubchem_cid")
    
    if(isCHEMBL(chembl_id)){
        chembl_id <- paste0("(\"", paste(chembl_id, collapse="\", \""), "\")")
        query <- dbSendQuery(conn, paste("SELECT a.chembl_id, a.drugbank_id,
            a.lincs_id, b.species, b.strain, b.dosage, b.avg_lifespan_change,
            b.max_lifespan_change, b.gender, b.significance,
            c.name, c.'cas-number', c.unii, c.state, c.groups, c.'fda-label',
            d.instance_id, d.'concentration..M.', d.'duration..h.',
            d.cell2, d.array3, d.perturbation_scan_id, d.vehicle_scan_id4,
            d.scanner, d.catalog_number,
            e.pert_iname, e.is_touchstone, e.inchi_key, e.pubchem_cid
            FROM id_mapping AS a
            LEFT JOIN drugAgeAnnot AS b ON a.drugage_id = b.drugage_id
            LEFT JOIN DrugBankAnnot AS c ON a.drugbank_id = c.drugbank_id
            LEFT JOIN cmapAnnot AS d ON a.cmap_id = d.cmap_id
            LEFT JOIN lincsAnnot AS e ON a.lincs_id = e.lincs_id
            WHERE a.chembl_id IN", chembl_id,
            "GROUP BY a.chembl_id
            ORDER BY a.chembl_id"))
        assays <- dbFetch(query)
        dbClearResult(query)
        result <- data.frame(assays)
        
        res <- result[, "chembl_id", drop=FALSE]
        for(x in annot){
            if(x == "drugAgeAnnot") res <- cbind(res, result[,da_cols])
            if(x == "DrugBankAnnot") res <- cbind(res, result[,db_cols])
            if(x == "cmapAnnot") res <- cbind(res, result[,cmap_cols])
            if(x == "lincsAnnot") res <- cbind(res, result[,lincs_cols])
            if(! x %in% c("drugAgeAnnot", "DrugBankAnnot", "cmapAnnot", "lincsAnnot")){
                idcol <- paste0(x, "_id")
                cust_annot <- dbGetQuery(conn, paste0("SELECT a.chembl_id, a.",
                                  idcol, ", b.* FROM id_mapping AS a LEFT JOIN ", x,
                                  " AS b ON a.", idcol, " = b.", idcol,
                                  " WHERE a.chembl_id IN", chembl_id,
                                  " GROUP BY a.chembl_id
                                  ORDER BY a.chembl_id"))
                cust_annot <- cust_annot[, !colnames(cust_annot) %in%
                                             c("chembl_id", idcol)]
                res <- cbind(res, cust_annot)
            }
        }
    } else {
        if(length(annot) > 1) 
            stop("The query ids passed to 'chembl_id' argument are not identified as ChEMBL IDs,
                 so only one annotation resource where the query id come frome is supported.")
        annot_tb <- dbReadTable(conn, annot)
        index <- apply(annot_tb, 1, function(r) any(chembl_id %in% r))
        res <- annot_tb[index,]
        if(annot == "drugAgeAnnot") res <- res[,da_cols]
        if(annot == "DrugBankAnnot") res <- res[,db_cols]
        if(annot == "cmapAnnot") res <- res[,cmap_cols]
        if(annot == "lincsAnnot") res <- res[,lincs_cols]
    }
    dbDisconnect(conn)
    return(res)
}

isCHEMBL <- function(id){
    return(all(grepl("CHEMBL", id)))
}

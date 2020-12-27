#‘ 
#' Get Compound Structures from Four Resources
#'
#' This function could be used to get SDFset of compounds in CMAP2, LINCS, DrugAge
#' or DrugBank databases. The \code{cid} of the SDFset are compound names instead
#' of their internal IDs.
#' 
#' @param source character(1), one of "CMAP2", "LINCS", "DrugBank", "DrugAge"
#' @return SDFset object of compounds in the \code{source} database, the \code{cid}
#' of the SDFset are compound names.
#' @seealso 
#' \code{\link[ChemmineR]{SDFset}}
#' @examples
#' da_sdf <- getSDFwithName(source="DrugAge")
#' @export
getSDFwithName <- function(source="LINCS"){
    ah <- AnnotationHub()
    annot_path <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), annot_path)
    if(source == "LINCS"){
        lincs_annot <- dbReadTable(conn, "lincsAnnot")
        lincs_sdf_path <- ah[["AH79567"]]
        lincs_sdf <-ChemmineR::read.SDFset(lincs_sdf_path)
        cid(lincs_sdf) <- sdfid(lincs_sdf)
        ## make lincs cids as pert names
        brd_ids <- lincs_annot$lincs_id
        names(brd_ids) <- lincs_annot$pert_iname
        brd_uniq <- brd_ids[!duplicated(names(brd_ids))] # 19811
        brd_common <- brd_uniq[brd_uniq %in% cid(lincs_sdf)] # 19758
        res_sdf <- lincs_sdf[brd_common]
        cid(res_sdf) <- names(brd_common)
    }
    if(source == "CMAP2"){
        cmap_annot <- dbReadTable(conn, "cmapAnnot")
        cmap_sdf_path <- ah[["AH79566"]]
        cmap_sdf <-ChemmineR::read.SDFset(cmap_sdf_path)
        cid(cmap_sdf) <- sdfid(cmap_sdf)
        cmap_ids <- cmap_annot$cmap_id
        names(cmap_ids) <- cmap_annot$cmap_name
        res_sdf <- cmap_sdf[cmap_ids]
        cid(res_sdf) <- names(cmap_ids)
    }
    if(source == "DrugBank"){
        db_annot <- dbReadTable(conn, "drugBankAnnot")
        db_sdf_path <- ah[["AH79565"]]
        db_sdf <-ChemmineR::read.SDFset(db_sdf_path)
        cid(db_sdf) <- sdfid(db_sdf)
        db_ids <- db_annot$drugbank_id
        names(db_ids) <- db_annot$name
        db_common <- db_ids[db_ids %in% cid(db_sdf)] # 10569
        res_sdf <- db_sdf[db_common]
        cid(res_sdf) <- names(db_common)
    }
    if(source == "DrugAge"){
        da_annot <- dbReadTable(conn, "drugAgeAnnot")
        da_sdf_path <- ah[["AH79564"]]
        da_sdf <-ChemmineR::read.SDFset(da_sdf_path)
        cid(da_sdf) <- sdfid(da_sdf)
        da_ids <- da_annot$drugage_id
        names(da_ids) <- da_annot$compound_name
        da_uniq <- da_ids[!duplicated(names(da_ids))] # 420
        da_common <- da_uniq[da_uniq %in% cid(da_sdf)] # 223
        res_sdf <- da_sdf[da_common]
        cid(res_sdf) <- names(da_common)
    }
    dbDisconnect(conn)
    return(res_sdf)
}

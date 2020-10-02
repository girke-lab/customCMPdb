#' Add/Delete Custom Annotation
#'
#' Functions could be used to add/delete user's custom compound annotations
#' to the annotation SQLite database in the
#' \code{\link[compoundCollectionData]{compoundCollectionData}} package.
#' The added custom compound annotation table should contains a column named as
#' \code{chembl_id} that represents the ChEMBL ids of the added compounds.
#'
#' @rdname customAnnot
#' @aliases addCustomAnnot customAnnot deleteAnnot listAnnot defaultAnnot
#' @param annot_tb data.frame representing the custom annotation table,
#' Note, it should contains a 'chembl_id' column representing the compound
#' ChEMBL ids
#' @param annot_name character(1), user defined name of the annotation table
#' @importFrom stats na.omit
#' @examples
#' chembl_id <- c("CHEMBL1000309", "CHEMBL100014", "CHEMBL10",
#'                "CHEMBL100", "CHEMBL1000", NA)
#' annot_tb <- data.frame(compound_name=paste0("name", 1:6),
#'         chembl_id=chembl_id,
#'         feature1=paste0("f", 1:6),
#'         feature2=rnorm(6))
#' addCustomAnnot(annot_tb, annot_name="mycustom3")
#' @export
addCustomAnnot <- function(annot_tb, annot_name){
    # check validity of annot_tb
    if(! "chembl_id" %in% colnames(annot_tb)){
        stop("The input annot_tb does not contain a 'chembl_id' column,
             Please make sure that your custom comppounds have ChEMBL ids!")
    }
    if(! any(grepl("CHEMBL", na.omit(annot_tb$chembl_id)))){
        stop("The ChEMBL ids are not in the correct format, please check!")
    }
    # check whether annot_name already exists in sqlite
    ah <- AnnotationHub()
    annot_path <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), annot_path)
    annot_names <- dbListTables(conn)
    if(tolower(annot_name) %in% tolower(annot_names)){
        ans <- readline(paste("This name has existed in the SQLite database,",
                            "do you want to overwrite it? (yes/no)"))
        if(tolower(ans)=="no" | ans==""){
            dbDisconnect(conn)
            return()
        }
    }
    ndigit <- nchar(as.character(nrow(annot_tb)))
    annot_tb2 <- data.frame(
        internal_id=paste0(toupper(annot_name),
            sprintf(paste0("%0", ndigit, "d"), seq_len(nrow(annot_tb)))),
        annot_tb
    )
    iid_name <- paste0(annot_name, "_id")
    colnames(annot_tb2)[1] <- iid_name
    chem2in <- na.omit(annot_tb2[ ,c("chembl_id", iid_name)])

    # write annotation table and id_mapping table to SQLite db
    id_mapping <- dbReadTable(conn, "id_mapping")
    id_mapping <- id_mapping[,colnames(id_mapping) != iid_name]
    id_mapping2 <- merge(id_mapping, chem2in, by="chembl_id",
                         all.x=TRUE, all.y=TRUE)
    dbWriteTable(conn, "id_mapping", id_mapping2, overwrite=TRUE)

    dbWriteTable(conn, annot_name,
                 annot_tb2,
                 overwrite=TRUE)
    message("The SQLite database now contains the following tables:\n",
            paste(dbListTables(conn), collapse=" "))
    dbDisconnect(conn)
}

#' @rdname customAnnot
#' @examples
#' deleteAnnot("mycustom3")
#' @export
deleteAnnot <- function(annot_name){
    if(tolower(annot_name) %in% c("cmapannot", "drugageannot",
                                  "drugbankannot", "lincsannot")){
        stop("The default annotation resources could not be deleted!")
    }
    ah <- AnnotationHub()
    annot_path <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), annot_path)
    annot_names <- dbListTables(conn)
    if(! annot_name %in% annot_names){
        dbDisconnect(conn)
        stop("The 'annot_name' does not exist in the SQLite database!")
    }
    id_map <- dbReadTable(conn, "id_mapping")
    # delete annot_name column in id_mapping table
    id_map <- id_map[,colnames(id_map) != paste(annot_name, "id", sep="_")]
    # remove rows that are all NAs except for chembl_id column
    tmp <- id_map[,colnames(id_map) != "chembl_id"]
    del_rows <- which(apply(is.na(tmp), 1, all))
    if(length(del_rows) > 0){
        id_map <- id_map[-del_rows, ]
    }
    dbWriteTable(conn, "id_mapping", id_map, overwrite=TRUE)
    # delete annot_name table
    dbRemoveTable(conn, annot_name)
    message("The SQLite database now contains the following tables:\n",
            paste(dbListTables(conn), collapse=" "))
    dbDisconnect(conn)
}

#' @description The \code{listAnnot} function lists the available annotation
#' resources in the SQLite annotation database.
#' @importFrom AnnotationHub AnnotationHub
#' @return character vector of names of the annotation tables in the SQLite DB
#' @rdname customAnnot
#' @examples
#' annot_names <- listAnnot()
#' @export
listAnnot <- function(){
    ah <- AnnotationHub()
    annot_path <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), annot_path)
    tb_names <- dbListTables(conn)
    annot_names <- tb_names[tb_names != "id_mapping"]
    dbDisconnect(conn)
    print(annot_names)
}

#' @rdname customAnnot
#' @description The \code{defaultAnnot} function sets the annotation SQLite
#' database to the default one by deleting the existing one and re-downloading
#' from AnnotationHub.
#' @return character(1), path to the annotation SQLite database
#' @examples
#' # defaultAnnot()
#' @export
defaultAnnot <- function(){
    ah <- AnnotationHub()
    annot_path <- ah[["AH79563"]]
    unlink(annot_path)
    annot_path <- ah[["AH79563"]]
    return(annot_path)
}

getidmap <- function(){
    ah <- AnnotationHub()
    annot_path <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), annot_path)
    id_map <- dbReadTable(conn, "id_mapping")
    dbDisconnect(conn)
    return(id_map)
}

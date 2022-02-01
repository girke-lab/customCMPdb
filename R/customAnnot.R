#' Load Compound Annotation Database
#' 
#' The compound annotation tables from different databases/sources are stored
#' in one SQLite database. This function can be used to load the SQLite 
#' annotation database
#' @examples 
#' conn <- loadAnnot()
#' @export
loadAnnot <- function(){
    ah <- AnnotationHub()
    # query(ah, c("customCMPdb", "annot_0.1"))
    annot_path <- ah[["AH79563"]]
    conn <- dbConnect(SQLite(), annot_path)
    return(conn)
}

#' Add/Delete Custom Annotation
#'
#' Functions could be used to add/delete user's custom compound annotations
#' from the annotation SQLite database.
#' The added custom compound annotation table should contains a column named as
#' \code{chembl_id} that represents the ChEMBL ids of the added compounds.
#'
#' @rdname customAnnot
#' @aliases addCustomAnnot customAnnot deleteAnnot listAnnot defaultAnnot
#' @param annot_tb data.frame representing the custom annotation table,
#' Note, it should contains a 'chembl_id' column representing the compound
#' ChEMBL ids
#' @param id_col column name in \code{annot_tb} that is used as ID column, 
#' this column must contain unique identifiers. If not defined, an automatically
#' generated ID column will be appended.
#' @param annot_name character(1), user defined name of the annotation table
#' @param overwrite a logical specifying whether to overwrite an existing table 
#' or not. Its default is FALSE.
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
addCustomAnnot <- function(annot_tb, id_col=NULL, annot_name, overwrite=FALSE){
    # check validity of annot_tb
    if(! "chembl_id" %in% colnames(annot_tb)){
        stop("The input annot_tb does not contain a 'chembl_id' column,
             Please make sure that your custom comppounds have ChEMBL ids!")
    }
    if(! isCHEMBL(annot_tb$chembl_id)){
        stop("The ChEMBL ids are not in the correct format, please check!")
    }
    # check whether annot_name already exists in SQLite
    conn <- loadAnnot()
    annot_names <- dbListTables(conn)
    if(tolower(annot_name) %in% tolower(annot_names) & !overwrite){
        # ans <- readline(paste("This name exists in the SQLite database,",
        #                     "do you want to overwrite it? (yes/no)"))
        # if(tolower(ans)=="no" | ans==""){
        #     dbDisconnect(conn)
        #     return()
        # }
        warning("The annot_name exists in the SQLite database, the old table is not overwritten, set 'overwrite=TRUE' to overwrite the existing one")
        dbDisconnect(conn)
        return()
    }
    ndigit <- nchar(as.character(nrow(annot_tb)))
    iid_name <- paste0(annot_name, "_id")
    if(is.null(id_col)){
        annot_tb <- data.frame(
            internal_id=paste0(toupper(annot_name),
                               sprintf(paste0("%0", ndigit, "d"), seq_len(nrow(annot_tb)))),
            annot_tb
        )
        colnames(annot_tb)[1] <- iid_name
    } else {
        if(sum(duplicated(annot_tb[[id_col]])) > 0 | sum(is.na(annot_tb[[id_col]])) > 0){
            stop("The id_col of annot_tb need to be unique and doesn't contain NA values!")
        }
        annot_tb <- data.frame(
            internal_id=annot_tb[[id_col]],
            annot_tb
        )
        colnames(annot_tb)[1] <- iid_name
    }
    
    chem2in <- na.omit(annot_tb[ ,c("chembl_id", iid_name)])

    # write annotation table and id_mapping table to SQLite db
    id_mapping <- dbReadTable(conn, "id_mapping")
    id_mapping <- id_mapping[,colnames(id_mapping) != iid_name]
    id_mapping2 <- merge(id_mapping, chem2in, by="chembl_id",
                         all.x=TRUE, all.y=TRUE)
    dbWriteTable(conn, "id_mapping", id_mapping2, overwrite=TRUE)

    dbWriteTable(conn, annot_name,
                 annot_tb,
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
                                  "drugbankannot", "lincsannot", "drugage4", "lincs2")){
        stop("The default annotation resources could not be deleted!")
    }
    conn <- loadAnnot()
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
    conn <- loadAnnot()
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
    conn <- loadAnnot()
    id_map <- dbReadTable(conn, "id_mapping")
    dbDisconnect(conn)
    return(id_map)
}

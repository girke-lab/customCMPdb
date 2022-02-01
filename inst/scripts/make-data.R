##########################################
## Create SQLite database and SDF files ##
##########################################
## Author: Yuzhu Duan
## Last update: 12-April-2020

# The SQLite database containing id mappings and drug annotation tables in
# DrugAge, DrugBank, CMAP2 and LINCS databases
getwd() # get user's current working directory. In my case, named as "compoundCollection"

##############
## DrugBank ##
##############

# get drugbank annotation table
## go to https://www.drugbank.ca/releases/latest download the drugbank xml 
## database as 'drugbank_5.1.5.xml.zip' file to `inst/scripts` directory
unzip("inst/scripts/drugbank_5.1.5.xml.zip", exdir="inst/scripts")
file.rename("inst/scripts/full database.xml", "inst/scripts/drugbank_5.1.5.xml")
library(customCMPdb)
drugbank_dataframe <- dbxml2df(xmlfile="inst/scripts/drugbank_5.1.5.xml", version="5.1.5")
df2SQLite(dbdf=drugbank_dataframe, version="5.1.5", dest_dir="./inst/scripts")
library(RSQLite)
dbpath <- "inst/scripts/drugbank_5.1.5.db"
conn <- dbConnect(SQLite(), dbpath)
dbListTables(conn)
dbdf <- dbGetQuery(conn, 'SELECT * FROM dbdf')
pryr::object_size(dbdf)
colnames(dbdf)

# get drugbank id to ChEMBL id mapping table
download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz",
              "inst/scripts/chembl2drugbank.txt.gz")
library(readr)
chem2db <- read_tsv("inst/scripts/chembl2drugbank.txt.gz") # 7290 X 2
colnames(chem2db) <- c("chembl_id", "drugbank_id")
length(intersect(dbdf$`drugbank-id`, chem2db$drugbank_id))
## 7290 / 13475 (54.1%) drugs in drugbank have chembl id
library(dplyr)
dbdf <- dbdf %>% rename("drugbank_id"="drugbank-id")

# download drugbank sdf file
## downlaod 'drugbank_all_structures.sdf.zip' file from 
## https://www.drugbank.ca/releases/latest#structures
## to the 'inst/scripts' directory
unzip("inst/scripts/drugbank_all_structures.sdf.zip", exdir="inst/scripts")
file.rename("inst/scripts/structures.sdf", "inst/scripts/drugbank_5.1.5.sdf")
db_stru <- read.SDFset("inst/scripts/drugbank_5.1.5.sdf") # 10,695

# Try to write drugbank id to Molecule_Name
dbids <- sapply(seq_along(db_stru), function(i){
    unlist(datablock(db_stru[[i]]))[1]
})
cid(db_stru) <- dbids
db_valid <- db_stru[validSDF(db_stru)] # 10,569
write.SDF(db_valid, "inst/scripts/drugbank_5.1.5.sdf", cid=TRUE)

#############
## DrugAge ##
#############

# Generate drugage.db by using 'buildDrugAgeDB' function
buildDrugAgeDB(dest_path="inst/scripts/drugage.db")
# Rename 'drugage.db' as 'compoundCollection.db' and add DrugBank compound 
# annotation and id_mappings to here
file.rename("inst/scripts/drugage.db", "inst/scripts/compoundCollection.db")
cc_path <- "inst/scripts/compoundCollection.db"
cc_conn <- dbConnect(SQLite(), cc_path)
dbListTables(cc_conn)
id_mapping <- dbGetQuery(cc_conn, "SELECT * FROM id_mapping")
## the 'chembl_id' column in DrugAge id mapping table contains several ChEMBL ids 
## in one row separated by comma
chem_vec <- id_mapping$chembl_id
names(chem_vec) <- id_mapping$drugage_id
chem_list <- lapply(chem_vec, function(x) unlist(strsplit(x, split=",")))
ida_dup <- rep(names(chem_list), times=sapply(chem_list, length))
id_mapping2 <- unique(data.frame(drugage_id=ida_dup, chembl_id=unlist(chem_list)))
# Add drugbank id to id mapping table
id_mapping3 <- merge(chem2db, id_mapping2, by="chembl_id", all.x=TRUE, all.y=TRUE) # 7448 X 3
dbWriteTable(cc_conn, "id_mapping", id_mapping3, overwrite=TRUE)
dbWriteTable(cc_conn, "DrugBankAnnot", dbdf)
dbListTables(cc_conn)
dbDisconnect(cc_conn)

# Get DrugAge sdf file
da_annot <- dbGetQuery(cc_conn, "SELECT * FROM drugAgeAnnot")
ida2pubchem <- na.omit(da_annot[,c("drugage_id", "pubchem_cid")])
pid_unique <- gsub(",.*", "", ida2pubchem$pubchem_cid)
library(ChemmineR)
da_cmp <- getIds(as.numeric(pid_unique))
cid(da_cmp) <- ida2pubchem$drugage_id
write.SDF(da_cmp, "inst/scripts/drugage.sdf", cid=TRUE)
## double check that the cids of the imported sdfset are drugage internal ids
da_cmp2 <- read.SDFset("inst/scripts/drugage.sdf")
head(cid(da_cmp2))
head(sdfid(da_cmp2))

###########
## CMAP2 ##
###########

# Get CMAP instance annotation table
download.file("https://portals.broadinstitute.org/cmap/cmap_instances_02.xls",
              "inst/scripts/cmap_instances_02.xls")
## Note, this file required some cleaning in LibreOffice (Excel would work for this too).
## The cleaning process removes the last lines that contain cell line info.
## After this it was saved as tab delimited txt file named 'cmap_instances_02.txt'
inst_path <- system.file("extdata/cmap_instances_02.txt", package="customCMPdb")
cmap_inst <- read.delim(inst_path, check.names=FALSE)
## The useful information in this cmap instance table are only compound names, 
## catelog number and catelog name

# Build 'cmap.db' 
## Since cmap only provides inconsistently formatted compound names and order 
## numbers. Obtaining the structures and PubChem IDs for CMAP was much harder 
## than expected. The 'buildCMAPdb' function builds 'cmap.db' that contains 
## id mapping and compound structure information.

## For about 2/3 of the CMAP drugs, one can obtain their PubChem/DrugBank IDs from 
## the DMAP site here: http://bio.informatics.iupui.edu/cmaps. 
## The SMILES strings for CMAP entries were obtained from ChemBank. Compounds 
## were matched by names using the 'stringdist' library where cmap_name from 
## CMAP were mapped to the closest name in ChemBank.
library(ChemmineR)
buildCMAPdb(dest_dir="inst/scripts")
conn <- dbConnect(SQLite(), "inst/scripts/cmap.db") # or conn <- initDb(mypath)
results <- getAllCompoundIds(conn)
sdfset <- getCompounds(conn, results, keepOrder=TRUE)
sdfset # An instance of "SDFset" with 1309 molecules
plot(sdfset[1:4], print=FALSE)

myfeat <- listFeatures(conn)
feat <- getCompoundFeatures(conn, results, myfeat)
feat[1:4,]

cmp_tb <- dbGetQuery(conn, 'SELECT * FROM compounds')
cmapAnnot <- as.data.frame(datablock2ma(datablock(sdfset))) # 1309 X 27
cmap_id <- paste0("CMAP", sprintf("%04d", as.integer(rownames(cmapAnnot))))
cmapAnnot2 <- data.frame(cmap_id=cmap_id, cmapAnnot, stringsAsFactors=FALSE)
# write CMAP SDF file
cid(sdfset) <- cmap_id
write.SDF(sdfset, file="inst/scripts/cmap.sdf", cid=TRUE)
cmap_sdf <- read.SDFset("inst/scripts/cmap.sdf")

# get pubchem cid to ChEMBL id mapping table
download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz",
              "inst/scripts/chembl2pubchemcid.txt.gz")
library(readr)
chem2pub <- read_tsv("inst/scripts/chembl2pubchemcid.txt.gz") # 1837162 X 2
colnames(chem2pub) <- c("chembl_id", "pubchem_cid")
chem2pub$pubchem_cid <- as.character(chem2pub$pubchem_cid)
# get cmap_id (created internal) to pubchem cid mappings
cmap2pub <- cmapAnnot2[,c("cmap_id", "PUBCHEM_ID")]
sum(cmap2pub$PUBCHEM_ID=="NA") # 443
## remove "NA"
cmap2pub <- cmap2pub[cmap2pub$PUBCHEM_ID != "NA",] # 866 X 2
dup <- cmap2pub[cmap2pub$PUBCHEM_ID %in% cmap2pub$PUBCHEM_ID[duplicated(cmap2pub$PUBCHEM_ID)],]
cmapAnnot2[cmapAnnot2$cmap_id %in% dup$cmap_id,]
### two cmap drugs with different names have the same pubchem cid, there are three this situation
# cmap_id to chembl id mapping
cmap2pub$PUBCHEM_ID <- gsub("CID", "", cmap2pub$PUBCHEM_ID)
library(dplyr)
cmap2chem <- cmap2pub %>% left_join(chem2pub, by=c("PUBCHEM_ID"="pubchem_cid"))
cmap2chem <- unique(na.omit(cmap2chem[,c("cmap_id", "chembl_id")])) # 800 X 2

cc_path <- "inst/scripts/compoundCollection.db"
cc_conn <- dbConnect(SQLite(), cc_path)
dbListTables(cc_conn)
id_mapping <- dbGetQuery(cc_conn, "SELECT * FROM id_mapping")
id_mapping2 <- merge(id_mapping, cmap2chem, by="chembl_id", all.x=TRUE, all.y=TRUE) # 7605 X 4
dbWriteTable(cc_conn, "id_mapping", id_mapping2, overwrite=TRUE)
dbWriteTable(cc_conn, "cmapAnnot", cmapAnnot2)
dbListTables(cc_conn)
dbDisconnect(cc_conn)

###########
## LINCS ##
###########

# get LINCS annotation table
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz"
download.file(url, "./inst/scripts/GSE92742_Broad_LINCS_pert_info.txt.gz")
pertIDs <- read.delim("inst/scripts/GSE92742_Broad_LINCS_pert_info.txt.gz")
pertIDs <- pertIDs[pertIDs$pert_type=="trt_cp", ] # 20,413 X 8
lincsAnnot <- pertIDs
colnames(lincsAnnot)[1] <- "lincs_id"

# get lincs internal id (e.g. BRD-A00100033) to chembl id mapping via inchi_key
## Download and extract the ChEMBL SQLite database from 
## ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_25/chembl_25_sqlite.tar.gz
## with 3.3GB. Replace the following path to the path of your downloaded ChEMBL db 
chembldb <- "~/insync/project/ChEMBL_data/chembl_25/chembl_25_sqlite/chembl_25.db"
library(RSQLite)
conn <- dbConnect(SQLite(), chembldb)
cmp_str <- dbGetQuery(conn, 'SELECT chembl_id, molregno, standard_inchi_key
                             FROM compound_structures, chembl_id_lookup
                             WHERE compound_structures.molregno = chembl_id_lookup.entity_id') # 3,099,738 X 3
lincs_inchi <- pertIDs[,c("pert_id", "inchi_key")]
lincs_inchi <- lincs_inchi[lincs_inchi$inchi_key != "-666",] # 20350 X 2
length(unique(lincs_inchi$inchi_key)) # 20307
library(dplyr)
lincs_inchi2 <- lincs_inchi %>% left_join(cmp_str[,c("chembl_id", "standard_inchi_key")],
                                          by=c("inchi_key"="standard_inchi_key"))
lincs2chem <- unique(na.omit(lincs_inchi2[,c("pert_id", "chembl_id")])) # 10,930 X 2
colnames(lincs2chem) <- c("lincs_id", "chembl_id")
length(unique(lincs2chem$lincs_id)) # 5774 unique lincs internal id
length(unique(lincs2chem$chembl_id)) # 10,878 unique chembl id
write_tsv(lincs2chem, "inst/scripts/lincs2chem.tsv")

# write annotation table and id_mapping table to SQLite db
cc_path <- "inst/scripts/compoundCollection.db"
cc_conn <- dbConnect(SQLite(), cc_path)
dbListTables(cc_conn)
id_mapping <- dbGetQuery(cc_conn, "SELECT * FROM id_mapping")
id_mapping2 <- merge(id_mapping, lincs2chem, by="chembl_id", all.x=TRUE, all.y=TRUE) # 17,052 X 5
dbWriteTable(cc_conn, "id_mapping", id_mapping2, overwrite=TRUE)
dbWriteTable(cc_conn, "lincsAnnot", lincsAnnot)
dbListTables(cc_conn)
dbDisconnect(cc_conn)

# write LINCS sdf file
library(ChemmineR)
library(ChemmineOB) # make sure openbabel module is loaded!!
pertIDs <- pertIDs[pertIDs$canonical_smiles!=-666, ]
smiset <- as.character(pertIDs$canonical_smiles)
pertids <- as.character(pertIDs$pert_id)
names(smiset) <- pertids
lincs_sdfset <- smiles2sdf(smiset)
valid <- validSDF(lincs_sdfset); lincs_sdfset <- lincs_sdfset[valid]
## note, original smiset contained 20350 and final valid SDFset 20333
saveRDS(lincs_sdfset, "inst/scripts/lincs_sdfset.rds")
lincs_sdfset <- readRDS("inst/scripts/lincs_sdfset.rds")
write.SDF(lincs_sdfset, file="inst/scripts/lincs.sdf")

# Save the 'cmap.sdf', 'lincs.sdf', 'drugage.sdf', 'drugbank.sdf', 'compoundCollection.db'
# files to '../ccdata/customCMPdb' directory and renamed as
# 'cmap02.sdf', 'lincs_pilot1.sdf', 'drugage_build2.sdf', 'drugbank_5.1.5.sdf'
# and 'compoundCollection_0.1.db'
# Upload the above files to AnnotationHub

###################################################
## Update the SQLite database to include compound annotation tables of DrugAge build 4 and LINCS beta 2020
###################################################

#### Update lincsAnnot to include more annotations from lincs_pert_info 2017 ####
#################################################################################

conn <- loadAnnot()
dbListTables(conn)
idmap <- dbReadTable(conn, "id_mapping")
cmpdb_lincs <- dbReadTable(conn, "lincsAnnot")
data("lincs_pert_info")
cmpdb_lincs2 <- tibble(left_join(cmpdb_lincs, select(lincs_pert_info, pert_id, MOA:selectivity_comment),
                          by=c(lincs_id="pert_id")))
dbWriteTable(conn, "lincsAnnot", cmpdb_lincs2, overwrite=TRUE)

#### Add DrugAge build 4 annotation table ####
##############################################

library(dplyr); library(readr); library(stringr); library(magrittr)
## clean DrugAge build 4 table
download.file("https://genomics.senescence.info/drugs/dataset.zip", 
              "~/insync/project/DrugAge/dataset_build4.zip")
unzip("~/insync/project/DrugAge/dataset_build4.zip", 
      exdir="~/insync/project/DrugAge/dataset_build4")
drugage <- read_csv("~/insync/project/DrugAge/dataset_build4/drugage.csv")

drugage_clean <- drugage %>% 
    select(-11, -12) %>%
    group_by(compound_name) %>% 
    summarise_all(paste, collapse="; ")

drugage_scrape <- lapply(seq_along(drugage_clean$compound_name), function(n) {
    cat("doing", n, "\n")
    system(paste0(
        'curl "https://genomics.senescence.info/drugs/drug_details.php?compound_name=',
        str_replace_all(drugage_clean$compound_name[n], " ", "%20"),
        '"'
    ), intern=TRUE) %>% 
        {
            cid <- if(length(cid_line <- str_which(., "<dt>PubChem CID")) > 0) {
                .[cid_line + 1] %>% str_extract("[0-9]+")
            } else NA
            iupac <- if(length(iupac_line <- str_which(., "<dt>IUPAC Name")) > 0) {
                .[iupac_line + 1] %>% str_remove("^[ ]+<dd>") %>% str_remove("&nbsp;</dd>\\\r$") 
            } else NA
            list(
                cid = cid,
                iupac = iupac
            )
        }
})

drugage_clean$pubchem_cid <- sapply(drugage_scrape, `[[`, "cid")
drugage_clean$iupac <- sapply(drugage_scrape, `[[`, "iupac")
drugage_clean$pubchem_cid[drugage_clean$compound_name == "1,2,4-triazolo[1,5-a]pyridine"] <- "67508"
drugage_clean$pubchem_cid[drugage_clean$compound_name == "7-cyclopentyl-5-(4-phenoxy)phenyl-7H-pyrrolo[2,3-d]pyrimidin-4-ylamine"] <- "6603792"
drugage_clean$pubchem_cid[drugage_clean$compound_name == "N-(2-[4-(4-chlorophenyl)piperazin-1-yl]ethyl)-3-methoxybenzamide"] <- "3626837"

## get pubchem cid to ChEMBL id mapping table
download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz",
              "~/insync/project/UniChem/chembl2pubchemcid.txt.gz") # download at 2022/01/31
chem2pub <- read_tsv("~/insync/project/UniChem/chembl2pubchemcid.txt.gz") # 2027837 X 2
colnames(chem2pub) <- c("chembl_id", "pubchem_cid")
chem2pub$pubchem_cid <- as.character(chem2pub$pubchem_cid)
drugage_clean %<>% left_join(chem2pub, by = "pubchem_cid")
drugage_clean %>% print(width=Inf)
write_tsv(drugage_clean, "~/insync/project/DrugAge/data/drugage4_annot.tsv")
## Add to SQLite
addCustomAnnot(drugage_clean, id_col="compound_name", annot_name="drugAge4")

da4 <- dbReadTable(conn, "drugAge4")
idmap <- getidmap()

#### Add curated LINCS beta 2020 compound annotation table ####
###############################################################
data("lincs_pert_info2")
addCustomAnnot(lincs_pert_info2, id_col="pert_id", annot_name="lincs2")

cmplincs2 <- dbReadTable(conn, "lincs2")
idmap <- getidmap()
sample_n(idmap, 50)
## copy and rename the cached SQLite
cache_path <- ah[["AH79563"]]
file.copy(cache_path, "~/insync/project/ccdata/customCMPdb/compoundCollection_0.1.db", 
          overwrite=TRUE)
dbDisconnect(conn)
## Double check
conn2 <- dbConnect(SQLite(), "~/insync/project/ccdata/customCMPdb/compoundCollection_0.1.db")
dbListTables(conn2)
idmap <- dbReadTable(conn2, "id_mapping")
sample_n(idmap, 50)
lincs_update <- dbReadTable(conn2, "lincsAnnot")
lincs2 <- dbReadTable(conn2, "lincs2")
da4 <- dbReadTable(conn2, "drugAge4")

## Update inst/extdata/metadata.csv, ask for Bioconductor core team to update the compoundCollection_0.1.db file.
## And upload the above files to AnnotationHub after their reply.
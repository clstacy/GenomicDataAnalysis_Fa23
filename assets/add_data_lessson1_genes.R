library(tidyverse)
library(here)

here::dr_here()

# pull original file with edgeR outputs
url <- "https://github.com/clstacy/GenomicDataAnalysis_Fa23/raw/main/data/12864_2020_6673_MOESM3_ESM.xlsx"

# create a new variable with the desired destination and name
destfile <- here("assets/12864_2020_6673_MOESM3_ESM.xlsx")

# use the curl package to download this file
curl::curl_download(url, destfile)

# read in the file
# gene_counts_rsem <- readxl::read_excel(destfile)
DE_results_heatshock <- readxl::read_excel(destfile, sheet = "edgeR QL",.name_repair = "universal")


# load in the GO term ontologies from geneontology.org for S. cerevesiae
sgd_gaf <- read_delim("http://current.geneontology.org/annotations/sgd.gaf.gz", 
                      delim = "\t", escape_double = FALSE, comment = "!", trim_ws = TRUE, 
                      col_names = c("DB", "Object ID", "Object Symbol",
                                    "Qualifier", "GO ID", "DB:Reference",
                                    "Evidence Code", "With (or) From",
                                    "Aspect", "Object Name", "Object Synonym",
                                    "Object Type", "Taxon", "Date",
                                    "Assigned By", "Annotation Extension",
                                    "Gene Product Form ID"),
                      name_repair= "universal", col_types=cols(Date = col_date("%Y%m%d"))
                      ) %>%
  mutate(ID = str_split_i(sgd_gaf$Object.Synonym, "\\|", 1))


# combine the two for student analysis
DE_results_heatshock_annotated <- sgd_gaf %>%
  dplyr::select(-c(Qualifier, DB.Reference, Evidence.Code, With..or..From, 
                   Aspect, Assigned.By, Annotation.Extension, 
                   Gene.Product.Form.ID, Date)) %>%
  group_by(ID, DB, Object.ID, Object.Symbol, Object.Name, Object.Synonym,
           Object.Type, Taxon) %>% 
  summarise(GO.ID = str_c(sort(GO.ID), collapse = "|"), 
            .groups = 'drop') %>%
  right_join(DE_results_heatshock, by = join_by("ID" == "ID")) %>%
  #rearrange columns for easy analysis
  dplyr::select(
    ID, Object.Symbol, Object.Name,
    Phenol.HS.logFC, Phenol.HS.p.value, Phenol.HS.FDR,
    RNeasy.HS.logFC, RNeasy.HS.p.value, RNeasy.HS.FDR,
    Direct.zol.HS.logFC, Direct.zol.HS.p.value, Direct.zol.HS.FDR,
    DB, Object.ID, Object.Synonym, Object.Type, Taxon, GO.ID
  )

# write file into github data folder
write_tsv(DE_results_heatshock_annotated,
          file = here("data/DE_results_heatshock_annotated.tsv"))



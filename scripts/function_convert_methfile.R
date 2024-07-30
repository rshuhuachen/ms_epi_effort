### This is a function to convert a united methylkit file 
### to a long file, where each line is a CpG site per individual
### and in the columns we have CpG site, nC, nT, cov, %meth

convert_meth <- function(methfile, novar){
  #extract data
  pacman::p_load(dplyr, methylKit)
  data <- methylKit::getData(methfile) #get methdata to dataframe
  
  if (novar == "remove"){
    data <- data[apply(data[,grep("numC", names(data), value = T)], 1, var,na.rm=T) != 0,]#no var in numC = all 0 = 0% meth
    data <- data[apply(data[,grep("numT", names(data), value = T)], 1, var,na.rm=T) != 0,]#no var in numT = all 0 = 100% meth
  }
  
  message(paste0("Out of ", nrow(methfile), " CpG sites, kept ", nrow(data), " which is ", round(nrow(data)/nrow(methfile), 2), "% removed"))
  
  data <- data %>% mutate(chr_pos = paste0(data$chr, "_", data$start), .before=chr) #create chr_pos column
  
  #wide to long df, one row per id per cpg site, seperate for C and T and then merge
  data_c <- cbind(data[,c(1:5)], data[,grepl("numC",names(data))]) #take only C cols
  data_t <- cbind(data[,c(1:5)], data[,grepl("numT",names(data))]) #take only T cols
  names(data_c)[6:ncol(data_c)] <- methfile@sample.ids #rename cols according to sample ID
  names(data_t)[6:ncol(data_c)] <- methfile@sample.ids#rename cols according to sample ID
  long_c <- gather(data_c, lib_id, numC, methfile@sample.ids[1]:methfile@sample.ids[length(methfile@sample.ids)], factor_key=TRUE) #wide to long, each row is one ID
  long_t <- gather(data_t, lib_id, numT,  methfile@sample.ids[1]:methfile@sample.ids[length(methfile@sample.ids)], factor_key=TRUE)
  
  #combine two long files
  long <- left_join(long_c, long_t, by = join_by(chr_pos, chr, start, end, strand, lib_id))
  ### add cols/rem cols
  long <- long %>% dplyr::select(-c(end, strand))
  long$cov <- long$numC + long$numT
  long$methperc <- long$numC / long$cov
  
  long <- long %>% mutate(epi_nr = gsub(".*_", "", lib_id),
                        lib = gsub("_.*", "", lib_id), .after =lib_id) # add epinr and library (splitting up lib_id)
  
  samplenr <- long %>% dplyr::select(c(epi_nr, lib)) %>% unique() %>% group_by(epi_nr) %>% mutate(n_sample = row_number())
  #in case there are repeats of the same sample
  
  long <- left_join(long, samplenr, by = c("epi_nr", "lib"))
  long$epi_nr <- as.integer(long$epi_nr)
  
  #add some other data from pheno file
  load("data/phenotypes/fulldata_complete_epi_withdates.RData")
  long <- left_join(long, all_pheno_epi[,c("epi_nr", "id", "year", "fulldate")], by = "epi_nr")
  
  return(long)}
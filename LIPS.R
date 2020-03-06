library('BiocManager')
library('genbankr')
library('seqinr')
library('ape')
library('rentrez')
library('GenomicRanges')
library('rtracklayer')
library('tidyverse')
library('ggplot2')


getwd()
setwd('/Users/david/Google Drive/PhD/LIPS/LIPS (R-code)')
set_entrez_key('2940c8fce95de8bb652be512772faf78a009')
accession_lassa <- as.vector(read_table('/Users/david/Google Drive/PhD/LIPS/Sequences/Lassa_complete_NP_accessionID.seq'))
accession_lassa_gb <- GBAccession(accession_lassa$Accession_ID)

# Accession Codes -------------------------------------------------------
Lassa_test_1 <- c('AF181853.1','KM821861.1','KM822082.1','MH887803.1',
                  'MK107932.1','MH053484.1','KM821847.1','AY628205.1')

Lassa_accession  <- c('AF181853.1','AF181854.1','AF246121.2','AF333969.1',
                     'AY179173.1','AY628201.1','AY628203.1','AY628205.1','AY628206.1',
                     'AY628207.1','AY628208.1','FR832711.1','GU481068.1','GU481070.1',
                     'GU481072.1','GU481074.1','GU481076.1','GU481078.1','GU830839.1',
                     'HQ688672.1','HQ688673.1','J04324.1','KF478765.1','KF478766.1',
                     'KM821773.1','KM821775.1','KM821777.1','KM821779.1','KM821781.1',
                     'KM821783.1','KM821785.1','KM821787.1','KM821789.1','KM821790.1',
                     'KM821792.1','KM821794.1','KM821796.1','KM821798.1','KM821800.1',
                     'KM821802.1','KM821804.1','KM821806.1','KM821808.1','KM821810.1',
                     'KM821812.1','KM821814.1','KM821815.1','KM821816.1','KM821818.1',
                     'KM821820.1','KM821822.1','KM821824.1','KM821826.1','KM821828.1',
                     'KM821830.1','KM821832.1','KM821833.1','KM821835.1','KM821837.1',
                     'KM821839.1','KM821841.1','KM821843.1','KM821845.1','KM821847.1',
                     'KM821848.1','KM821850.1','KM821852.1','KM821854.1','KM821856.1',
                     'KM821857.1','KM821859.1','KM821861.1','KM821863.1','KM821865.1',
                     'KM821867.1','KM821869.1','KM821871.1','KM821872.1','KM821874.1',
                     'KM821876.1','KM821877.1','KM821878.1','KM821880.1','KM821882.1',
                     'KM821884.1','KM821885.1','KM821886.1','KM821887.1','KM821889.1',
                     'KM821890.1','KM821892.1','KM821894.1','KM821895.1','KM821897.1',
                     'KM821899.1','KM821900.1','KM821902.1','KM821904.1','KM821906.1',
                     'KM821907.1','KM821909.1','KM821910.1','KM821912.1','KM821914.1',
                     'KM821916.1','KM821918.1','KM821920.1','KM821922.1','KM821924.1',
                     'KM821925.1','KM821927.1','KM821929.1','KM821931.1','KM821933.1',
                     'KM821935.1','KM821937.1','KM821939.1','KM821941.1','KM821943.1',
                     'KM821945.1','KM821947.1','KM821949.1','KM821951.1','KM821953.1',
                     'KM821955.1','KM821957.1','KM821959.1','KM821961.1','KM821963.1',
                     'KM821965.1','KM821967.1','KM821969.1','KM821971.1','KM821972.1',
                     'KM821974.1','KM821976.1','M821977.1','KM821979.1','KM821980.1',
                     'KM821982.1','KM821984.1','KM821986.1','KM821988.1','KM821990.1',
                     'KM821992.1','KM821994.1','KM821996.1','KM821998.1','KM822000.1',
                     'KM822002.1','KM822004.1','KM822006.1','KM822008.1','KM822010.1',
                     'KM822012.1','KM822014.1','KM822016.1','KM822018.1','KM822020.1',
                     'KM822022.1','KM822023.1','KM822025.1','KM822027.1','KM822029.1',
                     'KM822031.1','KM822033.1','KM822035.1','KM822037.1','KM822039.1',
                     'KM822041.1','KM822043.1','KM822045.1','KM822047.1','KM822049.1',
                     'KM822051.1','KM822053.1','KM822055.1','KM822057.1','KM822059.1',
                     'KM822060.1','KM822062.1','KM822063.1','KM822065.1','KM822066.1',
                     'KM822067.1','KM822069.1','KM822071.1','KM822072.1','KM822074.1',
                     'KM822076.1','KM822078.1','KM822080.1','KM822082.1','KM822084.1',
                     'KM822086.1','KM822088.1','KM822089.1','KM822091.1','KM822093.1',
                     'KM822095.1','KM822097.1','KM822099.1','KM822101.1','KM822103.1',
                     'KM822105.1','KM822107.1','KM822109.1','KM822111.1','KM822113.1',
                     'KM822115.1','KM822117.1','KM822118.1','KM822120.1','KM822122.1',
                     'KM822124.1','KM822126.1','KM822128.1','KM822130.1','KM822132.1',
                     'KT992416.1','KT992417.1','KT992418.1','KT992419.1','KT992420.1',
                     'KT992421.1','KT992422.1','KT992423.1','KT992424.1','KT992425.1',
                     'KU961971.1','KU978807.1','KU978810.1','KU978811.1','LC387468.1',
                     'LC387470.1','LC387472.1','LC387474.1','LC387476.1','LC387478.1',
                     'LC387480.1','LC387482.1','LC387484.1','LC387486.1','LC387488.1',
                     'MF990889.1','MG812631.1','MG812632.1','MG812635.1','MG812636.1',
                     'MG812639.1','MG812641.1','MG812643.1','MG812645.1','MG812647.1',
                     'MG812648.1','MG812650.1','MG812653.1','MG812655.1','MG812657.1',
                     'MG812658.1','MG812661.1','MG812663.1','MG812665.1','MG812667.1',
                     'MG812669.1','MG812670.1','MG812673.1','MG812675.1','MG812677.1',
                     'MG812679.1','MG812681.1','MG812683.1','MG812685.1','MH053463.1',
                     'MH053468.1','MH053470.1','MH053472.1','MH053474.1','MH053475.1',
                     'MH053476.1','MH053478.1','MH053480.1','MH053482.1','MH053483.1',
                     'MH053484.1','MH053485.1','MH053487.1','MH053488.1','MH053491.1',
                     'MH053493.1','MH053497.1','MH053498.1','MH053501.1','MH053502.1',
                     'MH053506.1','MH053511.1','MH053513.1','MH053514.1','MH053518.1',
                     'MH053522.1','MH053526.1','MH053527.1','MH053528.1','MH053531.1',
                     'MH053533.1','MH053536.1','MH053537.1','MH053538.1','MH053544.1',
                     'MH053545.1','MH053551.1','MH053552.1','MH053553.1','MH053554.1',
                     'MH053557.1','MH053559.1','MH053560.1','MH053561.1','MH053562.1',
                     'MH053564.1','MH053567.1','MH053572.1','MH053573.1','MH053574.1',
                     'MH053575.1','MH053576.1','MH053577.1','MH053578.1','MH053580.1',
                     'MH053582.1','MH053586.1','MH053588.1','MH053589.1','MH157028.1',
                     'MH157032.1','MH157034.1','MH157035.1','MH157037.1','MH157038.1',
                     'MH157041.1','MH157044.1','MH157045.1','MH157046.1','MH157050.1',
                     'MH215279.1','MH215281.1','MH215283.1','MH215285.1','MH215287.1',
                     'MH215289.1','MH887755.1','MH887757.1','MH887758.1','MH887759.1',
                     'MH887760.1','MH887761.1','MH887762.1','MH887763.1','MH887764.1',
                     'MH887765.1','MH887766.1','MH887767.1','MH887768.1','MH887769.1',
                     'MH887770.1','MH887771.1','MH887772.1','MH887773.1','MH887774.1',
                     'MH887775.1','MH887776.1','MH887777.1','MH887778.1','MH887779.1',
                     'MH887780.1','MH887781.1','MH887783.1','MH887784.1','MH887785.1',
                     'MH887786.1','MH887787.1','MH887788.1','MH887789.1','MH887790.1',
                     'MH887791.1','MH887792.1','MH887793.1','MH887795.1','MH887796.1',
                     'MH887797.1','MH887798.1','MH887799.1','MH887800.1','MH887801.1',
                     'MH887802.1','MH887803.1','MH887804.1','MH887805.1','MH887806.1',
                     'MH887807.1','MH887808.1','MH887809.1','MH887810.1','MH887811.1',
                     'MH887812.1','MH887813.1','MH887815.1','MH887816.1','MH887817.1',
                     'MH887818.1','MH887819.1','MH887820.1','MH887821.1','MH887822.1',
                     'MH887823.1','MH887824.1','MH887825.1','MH887826.1','MH887827.1',
                     'MH887828.1','MH887829.1','MH887830.1','MH887831.1','MH887832.1',
                     'MH887833.1','MH887834.1','MH887835.1','MH887836.1','MH887837.1',
                     'MH887838.1','MH887839.1','MH887840.1','MH887841.1','MH887842.1',
                     'MH887843.1','MH887844.1','MH887845.1','MH887846.1','MH887847.1',
                     'MH887848.1','MH887849.1','MH887851.1','MH887852.1','MH887853.1',
                     'MH887854.1','MH887855.1','MH887856.1','MH887857.1','MH887858.1',
                     'MH887860.1','MH887861.1','MH887863.1','MH887864.1','MH887865.1',
                     'MH887866.1','MH887867.1','MH887868.1','MH887869.1','MH887870.1',
                     'MH887871.1','MH887872.1','MH887873.1','MH887874.1','MH887875.1',
                     'MH887877.1','MH887878.1','MH887879.1','MH887880.1','MH887881.1',
                     'MH887882.1','MH887883.1','MH887884.1','MH887885.1','MH887886.1',
                     'MH887887.1','MH887888.1','MH887889.1','MH887890.1','MH887893.1',
                     'MH887894.1','MH887895.1','MH887896.1','MH887897.1','MK107897.1',
                     'MK107898.1','MK107899.1','MK107900.1','MK107901.1','MK107902.1',
                     'MK107903.1','MK107904.1','MK107905.1','MK107906.1','MK107907.1',
                     'MK107908.1','MK107909.1','MK107910.1','MK107911.1','MK107912.1',
                     'MK107913.1','MK107914.1','MK107915.1','MK107916.1','MK107917.1',
                     'MK107918.1','MK107919.1','MK107920.1','MK107921.1','MK107922.1',
                     'MK107923.1','MK107924.1','MK107925.1','MK107926.1','MK107927.1',
                     'MK107928.1','MK107929.1','MK107930.1','MK107931.1','MK107932.1',
                     'MK107933.1','MK107934.1','MK107935.1','MK107936.1','MK107937.1',
                     'MK107938.1','MK107939.1','MK107940.1','MK107941.1','MK107942.1',
                     'MK107943.1','MK107944.1','MK107945.1','MK107946.1','MK107947.1',
                     'MK107948.1','MK107949.1','MK107950.1','MK107951.1','MK107952.1',
                     'MK107953.1','MK107954.1','MK107955.1','MK107956.1','MK107957.1',
                     'MK107958.1','MK107959.1','MK107960.1','MK107961.1','MK107962.1',
                     'MK107963.1','MK107964.1','MK107965.1','MK107966.1','MK107967.1',
                     'MK107968.1','MK107969.1','MK107970.1','MK117846.1','MK117847.1',
                     'MK117848.1','MK117849.1','MK117851.1','MK117852.1','MK117853.1',
                     'MK117854.1','MK117855.1','MK117856.1','MK117857.1','MK117858.1',
                     'MK117859.1','MK117860.1','MK117861.1','MK117862.1','MK117863.1',
                     'MK117864.1','MK117865.1','MK117948.1','MK117949.1','MK117950.1',
                     'MK117951.1','MK117952.1','MK117953.1','MK117954.1','MK117956.1',
                     'MK117957.1','MK117958.1','MK117959.1','MK117960.1','MK117961.1',
                     'MK117962.1','MK117963.1','MK117964.1','MK117965.1','MK117966.1',
                     'MK117967.1','MK117968.1','MK117969.1','MK117970.1','MK117971.1',
                     'MK117972.1','MK117973.1','MK117974.1','MK117975.1','MK117976.1',
                     'MK117977.1','MK117978.1','MK117979.1','MK117980.1','MK117981.1',
                     'MK117982.1','MK117983.1','MK117984.1','MK117985.1','MK117986.1',
                     'MK117987.1','MK117988.1','MK117989.1','MK117990.1','MK117991.1',
                     'MK117992.1','MK117993.1','MK117994.1','MK117995.1','MK117996.1',
                     'MK117997.1','MK117998.1','MK117999.1','MK118000.1','MK118001.1',
                     'MK118002.1','MK118003.1','MK118005.1','MK118006.1','MK118007.1',
                     'MK118008.1','MK118009.1','MK118010.1','MK118011.1','MK118012.1',
                     'MK118013.1','MK118014.1','MK118015.1','MK118016.1','MK118017.1',
                     'MK118018.1','MK118019.1','MK118020.1','MK118021.1','MK118022.1',
                     'MK118023.1','MK118024.1','MK118025.1','MK118026.1','MK118027.1',
                     'MK118028.1','MK118029.1','MK118030.1','MK118031.1','MK118032.1',
                     'MK118033.1','MK118034.1','MK118035.1','MK118036.1','MK118037.1',
                     'MK118038.1','MK118039.1','MK281624.1','MK281625.1','MK281627.1',
                     'MK281631.1','MK291248.1','MK291249.1','MK291250.1','MN275173.1',
                     'NC_004296.1','X52400.1'
)


# Download Genbank based on accession numbers -----------------------------
Lassa_list <- list() #An empty list to be populated
for(entry in 1:length(accession_lassa_gb)){
  entry.name <- accession_lassa_gb[entry]
  Lassa_list[entry.name] <- readGenBank(GBAccession(accession_lassa_gb[entry]), partial = T)
}

# Create a list of identified sequences
Lassa_sequences <- data.frame('accession_lassa_gb'=accession_lassa_gb)
Lassa_sequences$host <- as.character(lapply(Lassa_list, function(x) x@sources@elementMetadata@listData$host))

# Append where and when they were obtained from
Lassa_sequences$country <- as.character(lapply(Lassa_list, function(x) x@sources@elementMetadata@listData$country))
Lassa_sequences$country <- gsub(':.*','\\',Lassa_sequences$country)
Lassa_sequences$year <- as.character(lapply(Lassa_list, function(x) x@sources@elementMetadata@listData$collection_date))
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}
Lassa_sequences$year <- RIGHT(Lassa_sequences$year, 4)

# Append the AA sequence for Nucleoprotein gene
Lassa_sequences$myaa <- Lassa_sequences$lmyse <- NA
for(entry in names(Lassa_list)) {
  if(any(Lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "NP")) {
    myposition <- which(Lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "NP")
    myaa <- as.character(Lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb==entry] <- myaa
    Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb==entry] <- nchar(myaa)
  } else {
    if(any(Lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "nucleoprotein")) {
      myposition <- which(Lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "nucleoprotein")
      myaa <- as.character(Lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb==entry] <- myaa
      Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb==entry] <- nchar(myaa)
    } else {
      if(any(Lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "nucleoprotein")) {
        myposition <- which(Lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "nucleoprotein")
        myaa <- as.character(Lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
        Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb==entry] <- myaa
        Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb==entry] <- nchar(myaa)
      } else {
        if(any(Lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "nucleocapsid protein")) {
          myposition <- which(Lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "nucleocapsid protein")
          myaa <- as.character(Lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
          Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb==entry] <- myaa
          Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb==entry] <- nchar(myaa)
        } else {
        if(any(Lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "hypothetical protein")) {
          myposition <- which(Lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "hypothetical protein")
          myaa <- as.character(Lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
          Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb==entry] <- myaa
          Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb==entry] <- nchar(myaa)
        } else { 
            myaa <- as.character(Lassa_list[["X52400.1"]]@cds@elementMetadata@listData$translation[2])
            Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb=="X52400.1"] <- myaa
            Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb=="X52400.1"] <- nchar(myaa)
          }
        }
      }
    }
  }
}

Lassa_sequences$n_a <- is.na(Lassa_sequences$myaa)
  
Lassa_sequences_NA <- subset(Lassa_sequences, country != "NULL" & year != "NA")
# Analysis ----------------------------------------------------------------

table(Lassa_sequences_NA$host, Lassa_sequences_NA$country)

ggplot(Lassa_sequences_NA, aes(x = year, fill = country))+
  geom_bar(stat = "count")+
  labs(title = "Number of LASV sequences and country of origin in NCBI GenBank",
       x = "Year", y = "Number of NCBI sequences")+
  theme(axis.text.x = element_text(angle = 90))

ggplot(Lassa_sequences_NA, aes(x = year, fill = host))+
  geom_bar(stat = "count")+
  labs(title = "Number of LASV sequences and host of viral sequence in NCBI GenBank",
       x = "Year", y = "Number of NCBI sequences")

# Identify source of rodent hosts ---------------------------------------
rodent_hosts <- as.data.frame(c(1:25))
if (Lassa_sequences$host != 'NULL' && Lassa_sequences$host != 'Homo sapiens'){
  rodent_hosts$accession <- Lassa_sequences$Accession_ID
}


# Interrogate the list of genbank files -----------------------------------
for(entry in 1:length(rodent_hosts$accession)){
  entry.name <- rodent_hosts[entry]
  rodent_lassa[entry.name] <- cds(Lassa_list[[entry.name]])
}



# Append the AA sequence for Nucleoprotein gene


rodent_hosts <- 1:25
if (Lassa_sequences$host != 'NULL' & Lassa_sequences$host != 'Homo sapiens')
  rodent_hosts$accession <- Lassa_sequences$Accession_ID
  

Rodent_lassa_list <- list() #An empty list to be populated
rodent_accession <- rodent_hosts$accession
for(entry in 1:length(rodent_accession)){
  entry.name <- rodent_accession[entry]
  Rodent_lassa_list[entry.name] <- readGenBank(GBAccession(rodent_accession[entry]), partial = T)
}


# Alternative method ------------------------------------------------------

# Use downloaded GFF file and import to R
file_lassa = file.path('/Users/david/Google Drive/Phd/Lips/Sequences/Lassa_complete_NP.gff3')
lseq = readGFF(fil)

# Drop data for the GP and non NP regions
lseq_np <- lseq[ which(lseq$type == 'CDS' & lseq$product == 'nucleoprotein'),]

# Identify the sequences that are derived from non-humans
lseq_non_h <- for(entry in lseq) {
  non_h_accession <- which(lseq(lseq$`nat-host` != 'NA' 'Homo sapiens')),
}

# Test code
Lassa_test <- list() #An empty list to be populated
for(entry in 1:length(Lassa_test_1)){
  entry.name <- Lassa_test_1[entry]
  Lassa_test[entry.name] <- readGenBank(GBAccession(Lassa_test_1[entry]), partial = T)
}

for(entry in 1:length(Lassa_test)){
  print(cds(Lassa_test[[entry]]))
}






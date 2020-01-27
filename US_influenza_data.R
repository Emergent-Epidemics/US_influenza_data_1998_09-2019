#SV Scarpino
#National ILI, flu lab tests, and total influenza cases
#Sept 2019

###########
#libraries#
###########
library(ggplot2)

###########
#acc funcs#
###########


###############
#Global Params#
###############
save_results <- FALSE
make_plots <- FALSE
time_stamp <- as.numeric(Sys.time())

######
#data#
######
ili <- read.csv("Data/ILINet.csv", skip = 1, stringsAsFactors = FALSE)
lab_old <- read.csv("Data/WHO_NREVSS_Combined_prior_to_2015_16.csv", skip = 1, stringsAsFactors = FALSE)
lab_new <- read.csv("Data/WHO_NREVSS_Clinical_Labs.csv", skip = 1, stringsAsFactors = FALSE)
pop <- read.csv("Data/us_pops.csv", skip = 1, stringsAsFactors = FALSE)
pops <- pop$Value[match(ili$YEAR, pop$Date)]*1000000 #US pop is in millions 

ili$X.UNWEIGHTED.ILI <- as.numeric(ili$X.UNWEIGHTED.ILI)
lab_old$PERCENT.POSITIVE <- as.numeric(lab_old$PERCENT.POSITIVE)
lab_new$PERCENT.POSITIVE <- as.numeric(lab_new$PERCENT.POSITIVE)
ili$WEEK <- as.numeric(ili$WEEK)

############
#Build data#
############
#1. Multiply ili by percent pos. from WHO labs
ili_raw_raw <- ili$X.UNWEIGHTED.ILI
ili_raw <- ili_raw_raw/100
lab_pos_raw <- c(lab_old$PERCENT.POSITIVE, lab_new$PERCENT.POSITIVE) #NOTE should add a code check here that the dates line up (checked manually on Oct. 7th 2019)
lab_pos <- lab_pos_raw/100
flu_ili_who_raw <- ili_raw * lab_pos
flu_ili_who <- flu_ili_who_raw * pops
flu_ili_who_corrected <- flu_ili_who * 0.7 #Rolfes et al. 2019, see ciz075.pdf in the lit folder, estimated that the 2017-2018 influenza season was ~70% lower than what one obtains from ili_pos * who_pos. NOTE - should return and use the actual values from Rolfes with confidence intervals

#2. build season indicator
season <- rep(NA, nrow(ili))
start <- 1997:2018
stop <- 1998:2019
for(i in 1:length(start)){
  sea.1.i <- which(ili$WEEK >= 40 & ili$YEAR == start[i])
  sea.2.i <- which(ili$WEEK <= 20 & ili$YEAR == stop[i])
  season[c(sea.1.i, sea.2.i)] <- paste0(start[i], "-", stop[i])
}

pan1 <- which(ili$WEEK >= 16 & ili$WEEK <= 33 & ili$YEAR == 2009) #note, should come back and check these against the literature
season[pan1] <- "2009 Pandemic First Wave"

pan2 <- which(ili$WEEK > 33 & ili$YEAR == 2009 | ili$WEEK <= 20 & ili$YEAR == 2010) #note, should come back and check these against the literature
season[pan2] <- "2009 Pandemic Second Wave"

#3. building data file
out <- data.frame(season, ili$YEAR, ili$WEEK, flu_ili_who, flu_ili_who_corrected, ili$X.UNWEIGHTED.ILI, lab_pos, pops)
colnames(out) <- c("Season", "Year", "Week", "Estimated_total_flu_US_ILI_WHO",  "Estimated_total_flu_US_ILI_WHO_RolfesCorrected", "ILI_per_pop_weighted", "Percent_Pos_Labs", "US_pop")

if(save_results == TRUE){
  #this is currently influenza_usa_sept2019.csv in the Results directory
  filename <- paste0(time_stamp, "influenza_usa.csv")
  write.csv(out, file = filename)
}

#4. Aggregate metrics
max_sea <- by(out$Estimated_total_flu_US_ILI_WHO_RolfesCorrected, out$Season, max, na.rm = T)
tot_sea <- by(out$Estimated_total_flu_US_ILI_WHO_RolfesCorrected, out$Season, sum, na.rm = T)

pop_out <- c(pop$Value[23:12], pop$Value[11], pop$Value[11], pop$Value[10:2])*1000000 #this is obviously bad and needs to be fixed
out2 <- data.frame(names(max_sea), as.numeric(max_sea), as.numeric(tot_sea), pop_out, as.numeric(tot_sea)/pop_out)
colnames(out2) <- c("season", "peak_cases", "total_cases", "population", "total_prop")

if(save_results == TRUE){
  #this is currently flu_season_summary_USA_sept2019.csv in the Results directory
  filename_season <- paste0(time_stamp, "flu_season_summary_usa.csv")
  write.csv(out2, file = filename_season, row.names = FALSE, quote = FALSE)
}

#######
#Plots#
#######
if(make_plots == TRUE){
  #I. Time series
  #quartz()
  y <- out$Estimated_total_flu_US_ILI_WHO_RolfesCorrected
  x <- 1:length(y) #should use package to transform CDC weeks in dates.
  
  plot(x, y, type = "l", main = "Influenza in United States 1998 - 2019", xlab = " ", ylab = "Estimated influenza cases", xaxt = "n", col = "#000000", lty = 1, lwd = 3, bty = "n")
  at.x <- (1:nrow(out))[which(out[,"Week"] == 1)]
  lab.x <- out[,"Year"][which(out[,"Week"] == 1)]
  axis(1, at = at.x, labels = lab.x)
  
  #II. Peak cases by season
  #quartz()
  par(mar = c(15, 4, 4, 2))
  bplt <- barplot(out2$peak_cases, xaxt = "n", ylab = "Peak number of influenza cases", main = "Peak influenza cases in the United States")
  labs <- as.character(out2$season)
  text(bplt, par("usr")[3]-0.25, labels = labs, srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex = 1) #this isn't quite right in terms of the label position, note the pandemics are too low.
  
  #III. Total prop infected by season
  #quartz()
  par(mar = c(15, 4, 4, 2))
  bplt2 <- barplot(out2$total_prop, xaxt = "n", ylab = "Total proportion of US infected with influenza", main = "Total proportion infected with influenza in the United States")
  labs <- as.character(out2$season)
  text(bplt2, 0, labels = labs, srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex = 1) #this isn't quite right in terms of the label position, note the pandemics are too low.
  
  #IV. Total cases by season
  #quartz()
  hist(out2$total_cases, col = "gray", xlab = "Total US influenza cases", main = "Total number of influenza cases in the US by season")
  abline(v = mean(out2$total_cases), col = "red", lwd = 3, lty = 3)
}
#SV Scarpino
#Flu lab tests and percent positive
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
lab_new <- read.csv("Data/WHO_NREVSS_Clinical_Labs.csv", skip = 1, stringsAsFactors = FALSE)

##########
#Analysis#
##########
states <- unique(lab_new$REGION)
years <- 2016:2019
week1 <- 50:52
week2 <- 1:3
slope1 <- matrix(NA, ncol = length(states), nrow = length(years))
colnames(slope1) <- states
row.names(slope1) <- years
slope1 <- as.data.frame(slope1)
slope2 <- slope1
ratio <- slope1
for(i in 1:length(states)){
  for(j in 1:length(years)){
    #slope for weeks 50:52
    use.i.j.1 <- which(lab_new$REGION == states[i] & lab_new$YEAR == years[j] & lab_new$WEEK %in% week1)
    if(length(use.i.j.1) != 3){
      stop()
    }
    ord.i.j.1 <- order(lab_new$WEEK[use.i.j.1])
    use.i.j.1.ord <- use.i.j.1[ord.i.j.1]
    slope1.i.j <- try(lm(as.numeric(lab_new$PERCENT.POSITIVE[use.i.j.1.ord]) ~ as.numeric(lab_new$TOTAL.SPECIMENS[use.i.j.1.ord])), silent = TRUE)
    if(is(slope1.i.j)[1] == "try-error"){
      slope1[j,i] <- NA
    }else{
      slope1[j,i] <- slope1.i.j$coefficients[2]
    }
    
    #slope for weeks 1:3
    use.i.j.2 <- which(lab_new$REGION == states[i] & lab_new$YEAR == (years[j]+1) & lab_new$WEEK %in% week2)
    if(length(use.i.j.2) != 3){
      stop()
    }
    ord.i.j.2 <- order(lab_new$WEEK[use.i.j.2])
    use.i.j.2.ord <- use.i.j.2[ord.i.j.2]
    slope2.i.j <- try(lm(as.numeric(lab_new$PERCENT.POSITIVE[use.i.j.2.ord]) ~ as.numeric(lab_new$TOTAL.SPECIMENS[use.i.j.2.ord])), silent = TRUE)
    
    if(is(slope2.i.j)[1] == "try-error"){
      slope2[j,i] <- NA
    }else{
      slope2[j,i] <- slope2.i.j$coefficients[2]
    }
  }
}

for(i in 1:nrow(ratio)){
  ratio[i,] <- as.numeric(slope1[i,])/as.numeric(slope2[i,])
}

wash <- which(colnames(ratio) == "Washington")
length(which(as.matrix(ratio) >= ratio["2019", wash]))/length(which(is.na(as.matrix(ratio)) == FALSE))
rowMeans(ratio, na.rm = TRUE)

hist(log(as.matrix(ratio)), col = "gray", main = "Change in slope % pos flu vs. testing WHO labs", xlab = "Ratio of slopes (weeks 50:52/1:3) (log scale)")
abline(v = ratio["2019", wash], col = "red", lwd = 3, lty = 3)
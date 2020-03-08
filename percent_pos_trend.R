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
lab_new$TOTAL.SPECIMENS <- as.numeric(lab_new$TOTAL.SPECIMENS)
lab_new$PERCENT.POSITIVE <- as.numeric(lab_new$PERCENT.POSITIVE)

##########
#Analysis#
##########
states <- unique(lab_new$REGION)
years <- 2017:2020
week1 <- 1:4
week2 <- 5:8
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
    if(length(use.i.j.1) != 4){
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
    use.i.j.2 <- which(lab_new$REGION == states[i] & lab_new$YEAR == (years[j]) & lab_new$WEEK %in% week2)
    if(length(use.i.j.2) != 4){
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

cali <- which(colnames(ratio) == "California")
length(which(as.matrix(ratio) >= ratio["2020", cali]))/length(which(is.na(as.matrix(ratio)) == FALSE))
rowMeans(ratio, na.rm = TRUE)

hist(log(as.matrix(ratio)), col = "gray", main = "Change in slope % pos flu vs. testing WHO labs", xlab = "Ratio of slopes (weeks 50:52/1:3) (log scale)")
abline(v = ratio["2020", cali], col = "red", lwd = 3, lty = 3)

use <- which(lab_new$WEEK %in% c(1:8) & lab_new$REGION %in% c("New York", "Massachusetts", "California", "Washington", "Oregon", "Texas"))
quartz()
ggplot(lab_new[use,], aes(as.numeric(TOTAL.SPECIMENS), as.numeric(PERCENT.POSITIVE), color = as.factor(YEAR))) + geom_point() + facet_wrap(~REGION)  + scale_color_brewer(palette = "Dark2") + ylab("Percent positive for influenza") + xlab("Epi Week") + theme(legend.position = "right", legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01))+labs(color = "Year")

use <- which(lab_new$WEEK %in% c(1:8) & lab_new$REGION %in% c("Washington") & lab_new$YEAR %in% c(2017, 2018, 2019, 2020))
quartz()
ggplot(lab_new[use,], aes(as.numeric(TOTAL.SPECIMENS), as.numeric(PERCENT.POSITIVE), color = as.factor(YEAR))) +  facet_wrap(~REGION) + geom_point() + scale_color_brewer(palette = "Dark2") + ylab("Percent positive for influenza") + xlab("Number of tests") + theme(legend.position = "right", legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01))+labs(color = "Year") + geom_smooth(method = "lm")

use <- which(lab_new$WEEK %in% c(1:8) & lab_new$REGION %in% c("Washington", "Oregon", "Idaho") & lab_new$YEAR %in% c(2017, 2018, 2019, 2020))
ggplot(lab_new[use,], aes(as.factor(WEEK), as.numeric(PERCENT.POSITIVE), fill = as.factor(YEAR))) + geom_bar(stat = "identity", position = "dodge") + facet_wrap(~REGION)  + scale_fill_brewer(palette = "Dark2") + ylab("Percent positive for influenza") + xlab("Epi Week") + theme(legend.position = "right", legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_y_continuous(expand = c(0.001,0.001))+ labs(fill = "Year")

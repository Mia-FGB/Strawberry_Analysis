library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

textsize <- c(8)
pointsize <- c(2)
pointalpha <-c(0.4)
pointshape <- c(1)
pointwidth <- c(1)
xvjust <- c(0.2)
yvjust <- c(0.8)

# Plot coverage vs position
airseq_data = read.table("error_bar_sequence_data.txt", header=TRUE)
airseq_data$Date <- as.Date(airseq_data$Date, format="%d/%m/%y")


plotit <- function(index, titletext, filename) {
    mini <- index+1;
    maxi <- index+2;
    xdata=airseq_data$Date
    ydata=airseq_data[,index]
    minerror=airseq_data[,mini];
    maxerror=airseq_data[,maxi];
    pdf(filename, width=3.5, height=2.5)
    gp <- qplot(xdata, ydata, geom="line") + geom_point(shape=1, size=2) +
      scale_x_date(limits = as.Date(c('2015-06-01','2015-07-31'))) + geom_line(color="black", size=0.1) +
      ggtitle(titletext) + theme(text = element_text(size=textsize)) +
      xlab("Date") + ylab("Hits per million reads") +
      theme(plot.margin = unit(c(0.02,0.02,0.04,0.02), "npc")) +
      theme(axis.title.y=element_text(vjust=0.2)) + theme(axis.title.x=element_text(vjust=-0.2)) +
      theme(plot.margin = unit(c(0.02,0.02,0.04,0.02), "npc")) + theme(axis.title.x=element_text(vjust=-xvjust)) +
      theme(axis.title.y=element_text(vjust=yvjust)) + geom_errorbar(aes(ymin=minerror,ymax=maxerror))
    print(gp)
    garbage <- dev.off()
    return(gp)
}


weather_data = read.table("/Users/leggettr/Desktop/Airseq/weather2.txt", header=TRUE)
weather_data$Date <- as.Date(weather_data$Date, format="%d/%m/%y")
pdf("/Users/leggettr/Desktop/temperature.pdf", width=6, height=1.5)
print(ggplot(weather_data, aes(x=weather_data$Date, y=weather_data$OutdoorTemperature)) + scale_x_date(limits = as.Date(c('2015-06-01','2015-07-31'))) + geom_line(color="red") + scale_y_continuous(limits=c(0, 40)) + theme(text = element_text(size=textsize)) + xlab("Date") + ylab("Temp/C") + theme(axis.title.y=element_text(vjust=0.2)) + theme(axis.title.x=element_text(vjust=-0.2)) + theme(plot.margin = unit(c(0.02,0.02,0.04,0.02), "npc")) + theme(axis.title.x=element_text(vjust=-xvjust)) + theme(axis.title.y=element_text(vjust=yvjust)))
garbage <- dev.off()

pdf("/Users/leggettr/Desktop/humidity.pdf", width=6, height=1.5)
print(ggplot(weather_data, aes(x=weather_data$Date, y=weather_data$OutdoorHumidity)) + scale_x_date(limits = as.Date(c('2015-06-01','2015-07-31'))) + geom_line(color="blue") + scale_y_continuous(limits=c(0, 100)) + theme(text = element_text(size=textsize)) + xlab("Date") + ylab("Humidity %") + theme(axis.title.y=element_text(vjust=0.2)) + theme(axis.title.x=element_text(vjust=-0.2)) + theme(plot.margin = unit(c(0.02,0.02,0.04,0.02), "npc")) + theme(axis.title.x=element_text(vjust=-xvjust)) + theme(axis.title.y=element_text(vjust=yvjust)))
garbage <- dev.off()

pdf("/Users/leggettr/Desktop/rainfall.pdf", width=6, height=1.5)
print(ggplot(weather_data, aes(x=weather_data$Date, y=weather_data$DayRainfall)) + scale_x_date(limits = as.Date(c('2015-06-01','2015-07-31'))) + geom_line(color="forestgreen") + theme(text = element_text(size=textsize)) + xlab("Date") + ylab("Rainfall/mm") + theme(axis.title.y=element_text(vjust=0.2)) + theme(axis.title.x=element_text(vjust=-0.2)) + theme(plot.margin = unit(c(0.02,0.02,0.04,0.02), "npc")) + theme(axis.title.x=element_text(vjust=-xvjust)) + theme(axis.title.y=element_text(vjust=yvjust)))
garbage <- dev.off()

ydata=airseq_data$Fusarium_culmorum
#pno <- plotit(21, "Phaeosphaeria nodorum", "/Users/leggettr/Desktop/Phaeosphaeria_nodorum.pdf")
#pst <- plotit(22, "Puccinia striiformis", "/Users/leggettr/Desktop/Puccinia_striiformis.pdf")
#pgr <- plotit(23, "Puccinia graminis", "/Users/leggettr/Desktop/Puccinia_graminis.pdf")
#psy <- plotit(24, "Pseudomonas syringae", "/Users/leggettr/Desktop/Pseudomonas_syringae.pdf")
#bgh <- plotit(25, "bgh_dh14", "/Users/leggettr/Desktop/bgh_dh14.pdf")
#bci <- plotit(26, "Botrytis cinerea", "/Users/leggettr/Desktop/Botrytis_cinerea.pdf")
#lam <- plotit(27, "Lambda", "/Users/leggettr/Desktop/Lambda.pdf")
#ztr <- plotit(28, "Zymoseptoria tritici", "/Users/leggettr/Desktop/Zymoseptoria_tritici.pdf")
#pre <- plotit(29, "Puccinia recondita", "/Users/leggettr/Desktop/Puccinia_recondita.pdf")
#ptr <- plotit(30, "Pyrenophora tritici", "/Users/leggettr/Desktop/Pyrenophora_tritici.pdf")
#fgr <- plotit(31, "Fusarium graminearum", "/Users/leggettr/Desktop/Fusarium_graminearum.pdf")
#bth <- plotit(32, "Bacillus thuringiensis", "/Users/leggettr/Desktop/Bacillus_thuringiensis.pdf")
fcu <- plotit(20, "Fusarium culmorum", "/Users/leggettr/Desktop/Fusarium_culmorum.pdf")
pno <- plotit(23, "Phaeosphaeria nodorum", "/Users/leggettr/Desktop/Phaeosphaeria_nodorum.pdf")
pst <- plotit(26, "Puccinia striiformis", "/Users/leggettr/Desktop/Puccinia_striiformis.pdf")
pgr <- plotit(29, "Puccinia graminis", "/Users/leggettr/Desktop/Puccinia_graminis.pdf")
psy <- plotit(32, "Pseudomonas syringae", "/Users/leggettr/Desktop/Pseudomonas_syringae.pdf")
bgh <- plotit(35, "bgh_dh14", "/Users/leggettr/Desktop/bgh_dh14.pdf")
bci <- plotit(38, "Botrytis cinerea", "/Users/leggettr/Desktop/Botrytis_cinerea.pdf")
lam <- plotit(41, "Lambda", "/Users/leggettr/Desktop/Lambda.pdf")
ztr <- plotit(44, "Zymoseptoria tritici", "/Users/leggettr/Desktop/Zymoseptoria_tritici.pdf")
pre <- plotit(47, "Puccinia recondita", "/Users/leggettr/Desktop/Puccinia_recondita.pdf")
ptr <- plotit(50, "Pyrenophora tritici", "/Users/leggettr/Desktop/Pyrenophora_tritici.pdf")
fgr <- plotit(53, "Fusarium graminearum", "/Users/leggettr/Desktop/Fusarium_graminearum.pdf")
bth <- plotit(56, "Bacillus thuringiensis", "/Users/leggettr/Desktop/Bacillus_thuringiensis.pdf")

pdf("/Users/leggettr/Desktop/all.pdf")
grid.arrange(fcu, pno, pst, pgr, psy, bgh, bci, lam, ztr, pre, ptr, fgr, bth, ncol = 3, main = "Pathogen hits")
garbage <- dev.off()

library(meta)      # Used for function metagen
library(gemtc)     # Used for function blobbogram
library(grid)

data <- read.table(textConnection('
  Study,      Year, OR,   CI.l, CI.u, N,      Interval,                                          group, location,            style
  "Ball",     2014, 1.07, 0.86, 1.34, "40441",  "IPI:  0 to 5 months (matched analysis)",           1,   "Perth",            "normal"
  "Ball",     2014, 1.04, 0.87, 1.23, "40441",  "IPI:  6 to 11 months (matched analysis)",          2,   "Perth",            "normal"
  "Ball",     2014, 0.87, 0.73, 1.03, "40441",  "IPI: 12 to 17 months (matched analysis)",          3,   "Perth",            "normal"
  "Schachar", 2016, 1.20, 1.13, 1.27, "302706", "IPI:  0 to 5 months (matched analysis)",           1,   "California",       "normal"
  "Schachar", 2016, 1.14, 1.08, 1.21, "302706", "IPI:  6 to 11 months (matched analysis)",          2,   "California",       "normal"
  "Schachar", 2016, 1.13, 1.07, 1.19, "302706", "IPI: 12 to 17 months (matched analysis)",          3,   "California",       "normal"
  "Hanley",   2017, 0.85, 0.71, 1.02, "38178",  "IPI:  0 to 5 months (matched analysis)",           1,   "British Columbia", "normal"
  "Hanley",   2017, 0.91, 0.79, 1.04, "38178",  "IPI:  6 to 11 months (matched analysis)",          2,   "British Columbia", "normal"
  "Hanley",   2017, 0.97, 0.84, 1.11, "38178",  "IPI: 12 to 17 months (matched analysis)",          3,   "British Columbia", "normal"
  "Howard",   2013, 1.15, 1.09, 1.20, "Not stated",  "IPI: not stratified (propensity analysis)",   4,   "Louisiana",        "normal"
'), header=TRUE, sep=",",strip.white=TRUE)

data$Study <- paste(data$Study, ", ",data$Year, sep="")

#Ln transformations
meta <- data
meta$OR.Ln <- log(meta$OR)
meta$CI.l.Ln <- log(meta$CI.l)
meta$CI.u.Ln <- log(meta$CI.u)

# OR conversion to effect size
# From:
# Chinn, 2000. From http://pubmed.gov/11113947
# Replicating Chinn using her data
effect.size <- log(1.32)/1.81
SE <-(log(5.68) - log(0.31))/2/1.96/1.81
# Now using our data
meta$CI.Ln_Width <- meta$CI.u.Ln - meta$CI.l.Ln
meta$SE.Ln <- meta$CI.Ln_Width/2/1.96 
meta$effect.size <- meta$OR.Ln/1.81
meta$effect.size.SE <- meta$SE.Ln/1.81
#SD from SE: http://handbook.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
# Not needed in this analysis
# meta$SD.Ln <- meta$SE.Ln * meta$N^0.5  #Varies with sample size
# meta$variance= meta$SD.Ln^2
# SMD = MEAN.DIFF / SD

#-------------------------
# Meta-analysis from library meta
meta1 <- metagen(meta$effect.size, meta$effect.size.SE, sm="OR", byvar=Interval, data=meta, comb.fixed=FALSE, backtransf = TRUE, studlab=meta$Study)
# WARNING: after backtransf for display in forest plots, note that point estimates and I.I.s do not exactly match the data inputted in the table 'data' above
forest(meta1,print.p=FALSE, xlim=c(0.5,2), leftcols=c("studlab"), print.tau2=FALSE,col.diamond="blue", col.diamond.lines="blue", print.I2.ci=TRUE,overall=FALSE)
#grid.text("Association of IPI and PTB in matched analyses", 0.5, 0.97, gp = gpar(fontsize = 14, fontface = "bold"))

#-------------------------
# Added meta-analysis results to dataframe for display in blobbogram
# The data below was manually built from the results of meta1 above
summary <- read.table(textConnection('
  Study      Year OR   CI.l CI.u N      Interval          group location style
  "Summary, I2 = 85% (55%-95%)"  NA 1.0220 0.9035 1.1559 NA  "IPI: 12 to 17 months (matched analysis)" 1 " " "pooled"
  "Summary, I2 = 78% (30%-93%)"  NA 1.0193 0.9392 1.1061 NA  "IPI:  6 to 11 months (matched analysis)" 2 " " "pooled"
  "Summary, I2 = 82% (44%-94%)"  NA 0.9999 0.9147 1.1094 NA  "IPI: 12 to 17 months (matched analysis)" 3 " " "pooled"
'), header=TRUE)

myframe <- NULL
myframe <- rbind(data,summary)
myframe$id <- myframe$Study #paste(myframe$intervention,', ', myframe$study, ' (', myframe$year, ')', sep='')
myframe$groupname <- myframe$Interval
myframe <- myframe[order(myframe$group),]
# Transforming data back from their natural logs
myframe$pe <- log(myframe$OR)
myframe$ci.l <- log(myframe$CI.l)
myframe$ci.u <- log(myframe$CI.u)

#-------------------------
#Blobbogram from library gemtc
blobbogram(myframe, group.labels=c('IPI: 0 to 5 months (matched analysis)', 'IPI: 6 to 11 months (matched analysis)', 'IPI: 12 to 17 months (matched analysis)', 'IPI: not stratified (propensity analysis)'),
           grouped=TRUE,
           columns=c('location', 'N'), column.labels=c('Location', 'Mothers'),
           id.label="Study", ci.label="OR (95% CI)", log.scale=TRUE, xlim=c(log(0.5),log(2)), digits=3)

mtext("Association of IPI and PTB in\nmatched or propensity analyses", side=3, line = -0.1, font=2, cex=1.75)

###########################################

# ALEDX2 #
# adapted from Christina Kumpitsch/ AG Moissl-Eichinger

###########################################


#install.packages("devtools")
#devtools::install_github("ggloor/ALDEx_bioc")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ALDEx2")


library(ALDEx2)

#################

## import data ##

#################

#prepare your otu table: remove all not-needed columns, arrange by groups 
#Import tab-del. file by the Import Dataset function; make sure column and row names are recognized as such

#you don't import an extra metadata file. Therefore you have to sort your data already in your otu.file. - in my case first 33 columns normosmic samples following 20 dysosmic
conds <- c(rep("tp1_complaints", 10), rep("tp2_complaints", 10))

#We recommend 128 or more mc.samples for the t-test, 1000 for a rigorous effect size calculation, and at least 16 for ANOVA.^[in fact we recommend that the number of samples in the smallest group multiplied by the number of DMC be equal at least 1000 in order to generate a reasonably stable estimate of the posterior distribution]
x.all <- aldex(tp1.2_genus_complaints, conds, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")

#The left panel is an Bland-Altman or MA plot that shows the relationship between (relative) Abundance and Difference. The right panel is an effect plot that shows the relationship between Difference and Dispersion. In both plots features that are not significant are in grey or black. Features that are statistically significant are in red. The Log-ratio abundance axis is the clr value for the feature.

###############################
# ALDEX clr module
#The workflow for the modular approach first generates random instances of the centred log-ratio transformed values. There are three inputs: counts table, a vector of conditions, and the number of Monte-Carlo (DMC) instances; and several parameters: a string indicating if iqlr, zero or all feature are used as the denominator is required, and level of verbosity (TRUE or FALSE). We recommend 128 or more mc.samples for the t-test, 1000 for a rigorous effect size calculation, and at least 16 for ANOVA.14
# the output is in the S3 object 'x'
x <- aldex.clr(tp1.2_genus_complaints, conds, mc.samples=128, denom="all", verbose=F)

#if you want to save this data
write.csv(x@analysisData, file = "genus_tp1-2_complaints_clr.csv")

###############################
#ALDEX ttest module
#The next operation performs the Welch's t and Wilcoxon rank test for the situation when there are only two conditions. There are only two inputs: the aldex object from aldex.clr and whether a paired test should be conducted or not (TRUE or FALSE).
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
write.csv(x.tt, file = "genus_tp1-2_complaints_clr_xtt.csv")

###############################
#ALDEX kw module
#Alternatively to the t-test, the user can perform the glm and Kruskal Wallace tests for one-way ANOVA of two or more conditions. Here there are only two inputs: the aldex object from aldex.clr, and the vector of conditions. Note that this is slow! and is not evaluated for this documentation.
x.kw <- aldex.kw(x)
write.csv(x.kw, file = "genus_tp1-2_complaints_clr_kw.csv")

###############################
#ALDEX effect module
#Finally, we estimate effect size and the within and between condition values in the case of two conditions. This step is required for plotting, in our lab we base our conclusions primarily on the output of this function1515 Macklaim et al. (2013);McMurrough et al. (2014);Bian et al. (n.d.).. There is one input: the aldex object from aldex.clr, ; and several parameters: a flag as to whether to include values for all samples or not are used as the denominator, and the level of verbosity. It is also possible to include the 95% confidence interval information for the effect size estimate with the flag CI=TRUE. This can be helpful when deciding whether to include or exclude specific features from consideration. We find that a large effect but that is an outlier for the extremes of the effect distribution can be false positives.
x.effect <- aldex.effect(x, CI=T, verbose=FALSE)
write.csv(x.effect, file = "genus_tp1-2_complaints_clr_xeffect.csv")

###############################
#ALDEX plot module
#plot
x.all <- data.frame(x.tt,x.effect)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch")
aldex.plot(x.all, type="MW", test="welch")
#The left panel is the MA plot, the right is the MW (effect) plot. In both plots red represents features called as differentially abundant with q <0.1; grey are abundant, but not differentially abundant; black are rare, but not differentially abundant. This function uses the combined output from the aldex.ttest and aldex.effect functions


#??? we.ep - Expected P value of Welch's t test
#??? we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#??? wi.ep - Expected P value of Wilcoxon rank test
#??? wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#??? kw.ep - Expected P value of Kruskal-Wallace test
#??? kw.eBH - Expected Benjamini-Hochberg corrected P value of Kruskal-Wallace test
#??? glm.ep - Expected P value of glm test
#??? glm.eBH - Expected Benjamini-Hochberg corrected P value of glm test
#??? rab.all - median clr value for all samples in the feature
#??? rab.win.NS - median clr value for the NS group of samples
#??? rab.win.S - median clr value for the S group of samples
#??? rab.X1 BNS.q50 - median expression value of features in sample X1 BNS if [include.item.summary=TRUE]
#??? dif.btw - median difference in clr values between S and NS groups
#??? dif.win - median of the largest difference in clr values within S and NS groups
#??? effect - median effect size: diff.btw / max(dif.win) for all instances
#??? overlap - proportion of effect size that overlaps 0 (i.e. no effect)

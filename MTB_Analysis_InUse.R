library("dplyr")
library("ggplot2")
library("ggpubr")
 
froot = 'C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/'
mtb <- read.csv('C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/MTB_ANALYSISSSS.csv')

mtb<-mtb %>% mutate_at(vars(Age,X,X12_Onset,X12_AB,X12_BA,Coh12,
                            X20_Onset,X20_AB,X20_BA,Coh20,GFP.Diff,
                 ITC1_300ms,ITC2_300ms,ITC3_300ms,ITC1.Avg,ITC2.Avg,
                 ITC3.Avg, CMR, Jane, MRT, GDT,WB.MEMR,HP.MEMR), as.numeric)
apply(mtb[which(sapply(mtb,is.numeric))],2,shapiro.test)#Shapiro Wilks across all columns (2 signifies columns)

mtb$X<-mtb$X*1e6
mtb$GFP.Diff <-mtb$GFP.Diff*1e6
mtb$Coh12 <-mtb$Coh12*1e6
mtb$Coh20 <-mtb$Coh20*1e6

mtb["Ages"] = cut(mtb$Age, c(0, 35, Inf), c("<35 y", ">35 y"), include.lowest=TRUE)
mtb$Ages = as.factor(mtb$Ages)

#Checking normality with Shapiro Wilk test first
with(mtb, shapiro.test(GDT[Age_group == ">35"]))

#Scatter plots with Pearson Co-eff, CI, and regression line 
jpeg("C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/GFPCoh20vsMST.jpeg",
     width = 4, height = 3, units = 'in', res = 500)

p<-ggscatter(data=mtb, x = "MRT", y = "GFP.Diff",
          # color = "black", shape = 21, size = 2, 
          add = "reg.line", add.params = list(color = "black", fill = "lightgray"),
          conf.int =FALSE,conf.int.level=0.95,cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson",label.x =-1.5,label.y=1,
                                label.sep = "\n"),
          xlab = "MRT Threshold (dB SNR)", ylab = expression("EEG Coherent Steady-State("~mu~"V)"), 
          use = "complete.obs", color="Ages",palette = c("#228B22", "#800080"),
          shape="Ages", size = 3)+
      theme(legend.position="top",
      legend.justification="right",
      legend.box.spacing=unit(0,"pt"),
      legend.margin=margin(0,0,0,0))

#dev.off()
ggsave(p, file="C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/GFPvsMRT4.4.jpeg", 
       width = 4, height = 3, units = "in", dpi=500)

#expression("EEG Coherence metric("~mu~"V)")
#expression("Baselined Steady-State("~mu~"V)")
#expression("EEG Coherence metric("~mu~"V)
#color="Age",show.legend.text = FALSE -- To show age on the plot
#ggplot(mtb,aes(Jane,GFP.Diff,color=Age))+
 # geom_point(size=2)+
  #geom_smooth(method=lm)+theme_classic()


#group_by(mtb, Age_group) %>%
#  summarise(
#    count = n(),
 #   mean = mean(ITC1.Avg, na.rm = TRUE),
#    sd = sd(ITC1.Avg,  na.rm = TRUE))
# ITC2.Avg, ITC3.Avg, GFP.Diff,CMR,
# GDT, Jane, MRT,MEMR,
# ITC2.Avg, ITC3.Avg, GFP.Diff,CMR,
# GDT, Jane, MRT,MEMR,

### F test
#var.test(ITC2.Avg ~ Age_group, data = mtb)

#dev.off()
#while (!is.null(dev.list()))  dev.off()
#print(plot(1))

library(ggpubr)

#t tests with box plots, p-value on it 
x <- which(names(mtb) == "Age_group") # name of grouping variable
y <- which(names(mtb) == "ITC1.Avg" # names of variables to test
           | names(mtb) == "ITC2.Avg" |
             names(mtb) == "ITC3.Avg" |
             names(mtb) == "GFP.Diff" |
             names(mtb) == "CMR"|
             names(mtb) == "GDT"|
             names(mtb) == "Jane"|
             names(mtb) == "MRT"|
             names(mtb) == "WB.MEMR"|
             names(mtb) == "HP.MEMR")
method <- "t.test" # one of "wilcox.test" or "t.test"
paired <- FALSE # if paired make sure that in the dataframe you have first all individuals at T1, then all individuals again at T2

for (i in y) {
  for (j in x) {
    ifelse(paired == TRUE,
           p <- ggpaired(mtb,
                         x = colnames(mtb[j]), y = colnames(mtb[i]),
                         color = colnames(mtb[j]), line.color = "gray", line.size = 0.4,
                         palette = "npg",
                         legend = "none",
                         xlab = colnames(mtb[j]),
                         ylab = colnames(mtb[i]),
                         add = "jitter"
           ),
           p <- ggboxplot(mtb,
                          x = colnames(mtb[j]), y = colnames(mtb[i]),
                          color = colnames(mtb[j]),
                          palette = "npg",
                          legend = "none",
                          add = "jitter"
           )
    )
    #  Add p-value
    print(p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                                 method = method,
                                 paired = paired,
                                 # group.by = NULL,
                                 ref.group = NULL
    ))
  }}
  # Save plots to jpeg Makes a separate file for each plot.
  ggsave(p, file=paste0("C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/plot_", i,".jpeg"), width = 4, height = 3, units = "in")

#Alternative box plots
##1
p<-ggboxplot(mtb,x = "Age", y = "Jane",
             color="Age",xlab="Age (in years)",
             ylab= "MST Threshold (dB SNR)",
             palette = c("#228B22", "#800080"),
             legend = "none",
             add = "jitter", shape="Age",
             add.params=list(size=3))
#+ geom_text(x=2, y=-25, label="p=8.5e-05")
ggsave(p, file="C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/JanevsAge4.3.1.jpeg", width = 3, height = 3, units = "in", dpi=500)

###2
ggplot(temp, aes(x=Ages, y=c(GDT)) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0)))

jpeg("C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/GDTvsAge1.jpeg",
     width = 4, height = 3, units = 'in', res = 500)

###3
bwplot(ITC~Ages|paste0("ITC",time),data=temp,
       layout=c(3,1),xlab='Age(in years)',
       ylab='Normalized ITC Values',strip=FALSE,
       par.settings=list(box.rectangle=list(col=c("#228B22", "#800080")),
                                       box.umbrella=list(col=c("#228B22", "#800080")),
                                       box.dot=list(col="black")))

dev.off()

####Multiple Linear Regression 
length(mtb$Jane) <- length(mtb$GFP.Diff)
length(mtb$ITC1.Avg) <- length(mtb$GFP.Diff)

model <- lm(mtb$Jane~mtb$GFP.Diff+mtb$ITC1.Avg, data=mtb)

plot("Jane", mtb$fitted.values,
     ylab = "Predicted Values",
     xlab = "Observed Values")
abline(a = 0, b = 1, lwd=2,
       col = "green")


ggscatter(x = "prediction",
          y = "actual",
          data = data.frame(prediction = predict(model),
                            actual = mtb$Jane), rm.na=TRUE) +
  geom_abline(intercept = 0,
              slope = 1) 


Y = a+XGFP.x1+XITC1.x2 #To predict using model by adding values at x1 and x2 


# Plotting multiple Regression Lines
ggplt+geom_smooth(method=lm,se=FALSE,fullrange=TRUE,
                  aes(color=Tree))
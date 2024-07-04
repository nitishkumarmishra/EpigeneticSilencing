library(ggpubr)
data("ToothGrowth")

#The most common methods for comparing means include:
  
#METHODS	R FUNCTION	DESCRIPTION
#T-test	          t.test()	          Compare two groups (parametric)
#Wilcoxon test	  wilcox.test()	      Compare two groups (non-parametric)
#ANOVA	          aov() or anova()  	Compare multiple groups (parametric)
#Kruskal-Wallis	  kruskal.test()	    Compare multiple groups (non-parametric)


#R functions to add p-values
#Here we present two new R functions in the ggpubr package:

#compare_means(): easy to use solution to performs one and multiple mean comparisons.
#stat_compare_means(): easy to use solution to automatically add p-values and significance levels to a ggplot.


##Create a box plot with p-values:
p <- ggboxplot(ToothGrowth, x = "supp", y = "len", color = "supp", palette = "jco", add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")
#Add p-values and significance levels to ggplots
p + stat_compare_means( aes(label = ..p.signif..), label.x = 1.5, label.y = 40) # This put na is pvalue >= 0.05
p + stat_compare_means( label = "p.signif", label.x = 1.5, label.y = 40) # specify the argument label as a character vector

### ******************************************** ###
#Compare more than two groups
# Global test
compare_means(len ~ dose,  data = ToothGrowth, method = "anova")
# Default method = "kruskal.test" for multiple groups
ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "jco")+
  stat_compare_means()
# Change method to anova
ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "jco")+
  stat_compare_means(method = "anova")


### ******************************************* ###
#Pairwise comparisons. If the grouping variable contains more than two levels, then pairwise tests will be performed automatically. The default method is "wilcox.test". You can change this to "t.test".
# Perorm pairwise comparisons
compare_means(len ~ dose,  data = ToothGrowth)
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

#Add p-values and significance levels to ggplots
#If you want to specify the precise y location of bars, use the argument label.y:
ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+ 
  stat_compare_means(label.y = 45)


### ******************************************** ###
#Multiple pairwise tests against a reference group:
# Pairwise comparison against reference
compare_means(len ~ dose,  data = ToothGrowth, ref.group = "0.5", method = "t.test")
# Visualize
ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5") # Pairwise comparison against reference


### ******************************************** ###
#Multiple pairwise tests against all (base-mean):
 # Comparison of each group against base-mean
compare_means(len ~ dose,  data = ToothGrowth, ref.group = ".all.", method = "t.test")
# Visualize
ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test", 
                     ref.group = ".all.") # Pairwise comparison against all

### ****************************************** ###
library(survminer)
data("myeloma", package = "survminer")
# Perform the test
compare_means(DEPDC1 ~ molecular_group,  data = myeloma, ref.group = ".all.", method = "t.test")

# Visualize the expression profile
ggboxplot(myeloma, x = "molecular_group", y = "DEPDC1", color = "molecular_group", add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(myeloma$DEPDC1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test", 
                     ref.group = ".all.") # Pairwise comparison against all

### ****************************************** ###
#Note that, if you want to hide the ns symbol, specify the argument hide.ns = TRUE.
# Visualize the expression profile
ggboxplot(myeloma, x = "molecular_group", y = "DEPDC1", color = "molecular_group", add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(myeloma$DEPDC1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test", 
                     ref.group = ".all.", hide.ns = TRUE) # Pairwise comparison against all

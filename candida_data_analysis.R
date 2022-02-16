library(tidyverse)
library(readxl)
library(ggpubr)
library(car)
library(emmeans)


candida <- read_excel("Experiment/candida_data.xlsx")

candida$Heat <- as.factor(candida$Heat)
candida$Metal <- as.factor(candida$Metal)

mod1 <- aov(Expression~Treatment, data=subset(candida, Target=="Vit"))
summary(mod1)
outlierTest(mod1)

plot(mod1)
hist(resid(mod1))

# Now perform post-hoc to determine which means have sig diff.
pairwise.t.test(subset(candida, Target=="Vit")$Expression, subset(candida, Target=="Vit")$Treatment,
                p.adjust.method="bonf")

TukeyHSD(mod1) # No significant differences in the mean vitrellogenin expression

# A further way of doing this is with the emmeans package which gives useful summaries.
vit_means <- emmeans(mod1, spec="Treatment")
summary(vit_means)

# Can choose other adjustment methods
pairs(vit_means, adjust="bonferroni")

# Plot the means
plot(vit_means, comparisons=TRUE)

library(ggsci) # colour palette library

p <- ggboxplot(subset(candida, Target=="Vit"), x="Treatment", y="Expression",
                   add="jitter", fill="Treatment", palette="npg", bxp.errorbar=TRUE,
               bxp.errorbar.width=0.2, show.legend=FALSE,
                   xlab="Treatment Group", ylab="Vg Expression Level")+
  stat_compare_means(method="anova",
                     label.x = 3, label.y = 1.5)
  stat_compare_means(aes(label=..p.signif..))
p

p + theme(legend.position = "none")

# Although it does appear that heat shocked individuals that had previously been exposed to metal 
# exposure did express lower levels of vitrellogenin, the difference in means was not significant.

# Now model differences in means of HsP expression across treatments

mod2 <- aov(log(Expression)~Treatment, data=subset(candida, Target=="HsP"))
summary(mod2)

plot(mod2)
hist(resid(mod2))

# Now perform post-hoc to determine which means have sig diff.
pairwise.t.test(subset(candida, Target=="HsP")$Expression, subset(candida, Target=="HsP")$Treatment,
                p.adjust.method="bonf")
TukeyHSD(mod2)

p2 <- ggboxplot(subset(candida, Target=="HsP"), x="Treatment", y="Expression",
          add="jitter", fill="Treatment", palette="npg", bxp.errorbar=TRUE,
          bxp.errorbar.width=0.2, show.legend=FALSE,
          xlab="Treatment", ylab="log(HSP70 Expression Level)")+
  stat_compare_means(method="anova",
                     label.x = 0.75, label.y = 1.4)
  stat_compare_means(aes(label=..p.signif..))
p2

p2+theme(legend.position = "none")

# Treatment groups 1 and 3 (that received heat shock) 
# had significantly higher HSP70 expression than the control and treatment group 2 (both no heat 
# shock) (all four p values < 0.01)
  
# multifactorial anova to look for interactions between heat and metal

# On vitellogenin expression
mod3 <- aov(Expression~Metal*Heat, data=subset(candida, Target=="Vit"))
summary(mod3)
plot(mod3)
hist(resid(mod3))
plot(resid(mod3) ~ fitted(mod3)) # Check homogeneity

TukeyHSD(mod3)


# On HsP expression
mod4 <- aov(log(Expression)~Metal*Heat, data=subset(candida, Target=="HsP"))
summary(mod4)

plot(mod4)
hist(resid(mod4))
plot(resid(mod4) ~ fitted(mod4))

TukeyHSD(mod4)

# correlation between HsP and vit

hsp_dat <- subset(candida, Target=="HsP")
vit_dat <- subset(candida, Target=="Vit")
treat_dat <- candida$Treatment[1:12]

new_df <- cbind(treat_dat, hsp_dat["Expression"], vit_dat["Expression"])
names(new_df)[1] <- "Treatment"
names(new_df)[2] <- "HsP"
names(new_df)[3] <- "Vit"
new_df

obs_6 <- new_df[6,]
obs_6

lm_plot <- ggplot(new_df, aes(x = HsP, y = Vit))+
  theme_bw()+
  geom_point(aes(color=Treatment))+
  geom_text(data=obs_6, label="T1.C", vjust=2, hjust=1, colour="red")+
  geom_smooth(method=lm, color="black")+
  labs(x="HSP70 Expression Level",
       y="Vg Expression Level")
lm_plot

lm_plot + theme(panel.grid.major=element_blank(), # Remove set of grid cells. Leaves large grid cells.
                panel.grid.minor=element_blank(), # With above, completely remove all grid cells.
                panel.border=element_blank(), # Remove borders of graph
                axis.line=element_line(),
                legend.title=element_blank())

express_mod <- lm(HsP~Vit, data=new_df)
summary(express_mod)

plot(express_mod)
hist(resid(express_mod))

express_mod2 <- lm(HsP~Vit, data=new_df[-6,])
summary(express_mod2)

plot(express_mod2)
hist(resid(express_mod2))

# linear regression model found HSP70 transcription positively correlated with Vg transcription 
# but the effect was not significant (t = 1.58, p = 0.15)

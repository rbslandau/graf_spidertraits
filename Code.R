# R Script Graf_et_al_2018_R_revised: R-Code
# Nadin Graf, Karina P. Battes, Mirela Cimpean, Pitt Dittrich,
# Martin H. Entling, Moritz Link, Andreas Scharmüller, Verena C. Schreiner,
# Eduard Szöcs, Ralf B. Schäfer
#----------------------------------------------------------------------------------------
# Do agricultural pesticides in streams affect riparian spiders?
#----------------------------------------------------------------------------------------
# submitted to 'Science of the total Environment'

# The code has been written by:
# Nadin Graf
# University of Koblenz-Landau
# Fortstrasse 7
# 76829 Landau
# GERMANY
# graf-nadin@uni-landau.de

# Revised by Moritz Link and Ralf B. Schäfer

# The code has been written to reproduce
# See 'attachment' for details

### Update March 27 2019 ##############################################################
# This script has been updated because of an error in the EPA Ecotox database         #
# The database incorrectly reported a value for atrazine a factor of 1000 too low     #
# EC50 for Deleatidium sp. is 0.6 mg/L not 0.6 µg/L as reported in the ECOTOX 		  #
# database (version December 2018) 													  #
# see Boonthai C., Scott R.R. & Chapman R.B. (2000). Acetylcholinesterase as a 		  #
# biomarker to assess the effect of chlorpyrifos and atrazine on some New Zealand 	  #
# aquatic invertebrates. Australasian Journal of Ecotoxicology 6, 59–64				  #
#######################################################################################

##### Changes in Graf et al. 2019  ###########################################
# The change has no impact on any of the results reported and conclusions 	 #
# drawn in our publication, except for negligible changes:					 #
# - slight changes in explained variance, p-values and correlations 		 #
# that did not impact any of the interpretations							 #
# - model for ballooning is now statistically significant, hence displays the#
# estimated relationship with toxicity										 #
# For reproducing the original version of the paper, use the variable 		 #
# max_sumTU_ms instead of max_sumTU_ms_corr in the models (needs slight 	 #
# adjustment of the code)													 #
##############################################################################

#######################
#      Structure      #
#######################

# Structure of the code:
# 1.  # ENVIRONMENTAL DATA SPCA ---------------
# 1.1 # proxies for agricultural land use -----
# 1.2 # sparse pca ----------------------------
# 2.  # RDA -----------------------------------
# 3.  # CWM -----------------------------------
# 4.  # CORRELATIONS --------------------------
# 4.1 # test for diversity --------------------
# 4.2 # test for abundance --------------------
# 4.3 # test for size -------------------------
# 4.4 # test for mobility ---------------------
# 4.5 # test for shading preference -----------
# 4.6 # test for moisture preference ----------
# 4.7 # test for specilisation ----------------
# 4.8 test for web  builder -------------------
# 4.9 # combine plots -------------------------
#----------------------------------------------

setwd("~/Gitprojects/Publications/graf_spidertraits")

#### 1. ENVIRONMENTAL DATA SPCA ####

lu_data <- read.csv("environmental_data.csv", header = TRUE, sep = ",")


# Aim of the SPCA is to create a sparse principal component
# that reflects the intensity of agricultural landuse at
# each sampling site and puts it into relation to the other
# sampling sites
# Decide, which variables are a proxy of agricultural land use.

head(lu_data)

#### 1.1 proxies for agricultural land use ####
# select variables that are suitable proxies for agricultural land use

proxy_vars <- c(
  "cl", "nh4", "no2", "no3", "po4", "so4", "shore_cover_shrubs",
  "shore_cover_meadow", "fid_factor", "land_use_meadows", "land_use_agri",
  "totalarea", "ext_agri", "int_agri"
)
## only keep the proxy variables
lu_dat <- lu_data[, proxy_vars]
str(lu_dat)

par(mfrow = c(4, 3))
for (nam in c(colnames(lu_dat))) {
  y <- lu_dat[, nam]
  plot(density(y), main = "", xlab = nam)
  plot(density(sqrt(y)), main = "", xlab = nam)
  plot(density(log10(y + 1)), main = "", xlab = nam)
}

#### 1.1 sparse pca ####
library(vegan)
library(pcaPP)

# set 2 as k.max
k.max <- 2
set.seed(1000)
oTPO_D1 <- opt.TPO(scale(lu_dat), k.max = k.max, method = "sd")
summary(oTPO_D1$pc)
oTPO_D1$pc$load

#  Tradeoff Curves: Explained Variance vs. sparseness
par(mfrow = c(1, k.max))
for (i in 1:k.max) plot(oTPO_D1, k = i)
# L0 gives number of variables with zero loadings on PC
# lambda opt is the estimated optimal penalty
# ECV1 is the empirical cumulated variance of the PCs

##  Tradeoff Curves: Explained Variance vs. lambda
par(mfrow = c(1, k.max))
for (i in 1:k.max) plot(oTPO_D1, k = i, f.x = "lambda")
# lambda represents penalty term

# Extract the spc scores
spc_D1 <- sPCAgrid(scale(lu_dat), k = k.max, lambda = oTPO_D1$pc.noord$lambda, method = "sd")
load_spcD1 <- scores(spc_D1, choices = 1, display = "species", scaling = 0)

pca_axes_D1 <- data.frame(vegan::scores(spc_D1, disp = "sites", choices = c(1:ncol(lu_dat)), scaling = "sites"))
pca_axes_D1$site <- lu_data$site
pca_axes_D1
# write.csv(pca_axes_D1, 'C:/Users/N/PhD/Digital/Rumänien/trait/txt/SPC_axis.csv', row.names = FALSE)
par(mfrow = c(1, 1))
biplot(spc_D1)


#### 2. RDA ########
#### 2.1 RDA all ####
library(ggplot2)
library(vegan)

data <- read.csv("stat_trait.csv", header = TRUE, sep = ",")
env <- data[, -c(1, 6:20)]
rownames(env) <- c("A", "B", "C", "E", "F", "G", "H", "I", "K", "L", "N", "O", "P", "Q", "R", "S", "T")
spe <- read.csv("allspiders.csv", header = TRUE, sep = ",")
spe <- spe[, -c(1)] # exclude name
rownames(spe) <- c("A", "B", "C", "E", "F", "G", "H", "I", "K", "L", "N", "O", "P", "Q", "R", "S", "T")
names(spe)

set.seed(1000)
expl <- rda(formula = spe ~ max_sumTU_ms_corr + shading + Comp.1 + Comp.2, data = env, scale = TRUE)
expl
summary(expl)
set.seed(222)
anova(expl, by = "terms")
anova(expl, by ='axis')

set.seed(1000)
specno <- specnumber(spe, MARGIN = 2)
svg(filename = "rdaplotall_review.svg", width = 10, height = 8, pointsize = 12)
plot(expl,
  type = "n", scaling = 2,
  xlab = "RDA1 (11% of total variance)", ylab = "RDA2 (7% of total variance)", las = 1, xlim = c(-1, 1), ylim = c(-1.5, 1)
)
orditorp(expl,
  scaling = 1, cex = 1, display = "species", col = "darkred", priority = specno,
  air = 2, pch = "+"
)
ef <- envfit(expl, env)
plot(ef, col = "cornflower blue")
dev.off()

#### 2.1 RDA tox ####
set.seed(1000)
tox <- rda(
  formula = spe ~ max_sumTU_ms_corr + Condition(shading + Comp.1 + Comp.2),
  data = env, scale = TRUE
)
summary(tox)
set.seed(1000)
anova(tox, by = "terms")
plot(tox)

sort(scores(tox)$species[, 1])

tox1 <- rda(
  formula = spe ~ max_sumTU_ms_corr + Condition(Comp.1 + Comp.2),
  data = env, scale = TRUE
)
set.seed(1000)
anova(tox1, by = "terms")


#### 3. CWM ####
library(FD)

trait <- read.csv("trait.csv", header = TRUE, row.names = 1, sep = ",")
names(trait)
identical(row.names(trait), row.names(t(spe))) # identical

ex1 <- dbFD(trait, spe, corr = "lingoes") # lingoes
ex1 # warnings due to rows with the same values
#write.csv(ex1, file = "cwm.csv")

#### 4. CORRELATIONS ####

library(vegan)

# correlaton between explanatories
cor.test(data$Comp.1, data$Comp.2)
cor.test(data$Comp.1, data$shading)
cor.test(data$Comp.1, data$max_sumTU_ms_corr)
cor.test(data$Comp.2, data$shading)
cor.test(data$Comp.2, data$max_sumTU_ms_corr)
cor.test(data$shading, data$max_sumTU_ms_corr)

#### 4.1 test for number of spiders #####
sp2.1 <- glm(spe ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "poisson")
summary(sp2.1)
drop1(sp2.1, test = "Chi") #-> exclude Comp.1

sp2.2 <- glm(spe ~ Comp.2 + max_sumTU_ms_corr + shading, data, family = "poisson")
summary(sp2.2)
drop1(sp2.2, test = "Chi") #-> exclude comp2

sp2.3 <- glm(spe ~ shading + max_sumTU_ms_corr, data, family = "poisson")
summary(sp2.3)
drop1(sp2.3, test = "Chi") #-> exclude shading

sp2.4 <- glm(spe ~ max_sumTU_ms_corr, data, family = "poisson")
summary(sp2.4)
1 - (8.6051 / 14.4198) #-> D² = 0.40

range(data$max_sumTU_ms_corr, na.rm = TRUE)
df <- data.frame(max_sumTU_ms_corr = seq(-1.62549514, -0.01420084, length.out = 10))
p_df <- predict(sp2.4,
  newdata = df,
  type = "response",
  se.fit = TRUE
)
tr_df <- transform(df,
  fit = p_df$fit,
  lwr = p_df$fit - (1.96 * p_df$se.fit),
  upr = p_df$fit + (1.96 * p_df$se.fit)
)

spe <- ggplot(tr_df, aes(y = fit, x = max_sumTU_ms_corr)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = data, aes(y = spe, x = max_sumTU_ms_corr), size = 2) +
  labs(
    x = expression("toxicity [max log sum TU]"),
    y = expression("species richness")
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"
  ) +
  xlim(NA, 0.1)
spe


#### 4.2 test for abundance #####
in1.1 <- glm(ind ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "quasipoisson")
summary(in1.1)
drop1(in1.1, test = "F") #-> exclude Comp.2

in1.2 <- glm(ind ~ Comp.1 + max_sumTU_ms_corr + shading, data, family = "quasipoisson")
summary(in1.2)
drop1(in1.2, test = "F") # exclude shading

in1.3 <- glm(ind ~ Comp.1 + max_sumTU_ms_corr, data, family = "quasipoisson")
summary(in1.3)
drop1(in1.3, test = "F") # exclude Comp.1

in1.4 <- glm(ind ~ max_sumTU_ms_corr, data, family = "quasipoisson")
summary(in1.4)
1 - (113.05 / 186.375) # D² = 0.39

range(data$max_sumTU_ms_corr, na.rm = TRUE)
df <- data.frame(max_sumTU_ms_corr = seq(-1.62549514, -0.01420084, length.out = 10))
p_df <- predict(in1.4,
  newdata = df,
  type = "response",
  se.fit = TRUE
)
tr_df <- transform(df,
  fit = p_df$fit,
  lwr = p_df$fit - (1.96 * p_df$se.fit),
  upr = p_df$fit + (1.96 * p_df$se.fit)
)

ind <- ggplot(tr_df, aes(y = fit, x = max_sumTU_ms_corr)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = data, aes(y = ind, x = max_sumTU_ms_corr), size = 2) +
  labs(
    x = expression("toxicity [max log sum TU]"),
    y = expression("number of individuals")
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"
  ) +
  ylim(NA, 200) +
  xlim(NA, 0.1)
ind

#### 4.3 test for size #####
siz1.1 <- glm(CWM.mKLw ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(siz1.1)
drop1(siz1.1, test = "F") #-> exclude Comp.2

siz1.2 <- glm(CWM.mKLw ~ Comp.1 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(siz1.2)
drop1(siz1.2, test = "F") #-> exclude Comp.1

siz1.3 <- glm(CWM.mKLw ~ max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(siz1.3)
drop1(siz1.3, test = "F") #-> exclude tox

siz1.4 <- glm(CWM.mKLw ~ shading, data, family = "gaussian")
summary(siz1.4)
1 - (6.0079 / 9.4986) # D²=0.3674963

range(data$shading, na.rm = TRUE)
range(data$Comp.2, na.rm = TRUE)
range(data$tox_median, na.rm = TRUE)
df <- data.frame((shading <- seq(0, 95, length.out = 10)))
p_df <- predict(siz1.4,
  newdata = df,
  type = "response",
  se.fit = TRUE
)
tr_df <- transform(df,
  fit = p_df$fit,
  lwr = p_df$fit - (1.96 * p_df$se.fit),
  upr = p_df$fit + (1.96 * p_df$se.fit)
)
siz <- ggplot(tr_df, aes(y = fit, x = shading)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = data, aes(y = CWM.mKLw, x = shading), size = 2) +
  labs(
    x = expression("shading [%]"),
    y = expression("CWM for body size")
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"
  ) +
  ylim(NA, 7) +
  xlim(NA, 100)
siz


#### 4.4 test for mobility  #####

bal1.1 <- glm(CWM.Ballon1 ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(bal1.1)
drop1(bal1.1, test = "F") #-> exclude shading

bal1.2 <- glm(CWM.Ballon1 ~ Comp.1 + Comp.2 + max_sumTU_ms_corr, data, family = "gaussian")
summary(bal1.2)
drop1(bal1.2, test = "F") #-> exclude Comp.2

bal1.3 <- glm(CWM.Ballon1 ~ Comp.1 + max_sumTU_ms_corr, data, family = "gaussian")
summary(bal1.3)
drop1(bal1.3, test = "F") #-> exclude Comp.1

bal1.4 <- glm(CWM.Ballon1 ~ max_sumTU_ms_corr, data, family = "gaussian")
summary(bal1.4)
drop1(bal1.3, test = "F")
1 - (0.70088  / 0.92828) # D² = 0.245

## Update 27 March 2019: Mobility shows now a significant relationship with toxicity

#### 4.5 test for shading preference ####

sha1.1 <- glm(CWM.ShaPosS ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(sha1.1)
drop1(sha1.1, test = "F") #-> exclude comp.1

sha1.2 <- glm(CWM.ShaPosS ~ Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(sha1.2)
drop1(sha1.2, test = "F") #-> exclude comp.2

sha1.3 <- glm(CWM.ShaPosS ~ max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(sha1.3)
drop1(sha1.3, test = "F") #-> exclude tox

sha1.4 <- glm(CWM.ShaPosS ~ shading, data, family = "gaussian")
summary(sha1.4)
1 - (0.007834 / 0.014009) # D²=0.4407881

range(data$shading, na.rm = TRUE)
df <- data.frame(shading = seq(0, 95, length.out = 10)) # ,

p_df <- predict(sha1.4,
  newdata = df,
  type = "response",
  se.fit = TRUE
)
tr_df <- transform(df,
  fit = p_df$fit,
  lwr = p_df$fit - (1.96 * p_df$se.fit),
  upr = p_df$fit + (1.96 * p_df$se.fit)
)
sha <- ggplot(tr_df, aes(y = fit, x = shading)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = data, aes(y = CWM.ShaPosS, x = shading), size = 2) +
  labs(
    x = expression("shading [%]"),
    y = expression("CWM for shading")
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"
  ) +
  ylim(NA, 0.4) +
  xlim(NA, 100)
sha

#### 4.6 test for moisture preference #####
moi1.1 <- glm(CWM.MoisPosS ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(moi1.1)
drop1(moi1.1, test = "F") #-> exclude Comp.2

moi1.2 <- glm(CWM.MoisPosS ~ Comp.1 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(moi1.2)
drop1(moi1.2, test = "F") #-> exclude Comp.1

moi1.3 <- glm(CWM.MoisPosS ~ max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(moi1.3)
drop1(moi1.3, test = "F") #-> exclude max_sumTU_ms_corr

moi1.4 <- glm(CWM.MoisPosS ~ shading, data, family = "gaussian")
summary(moi1.4)
1 - (0.017706 / 0.020364) # D²=0.1305245


#### 4.7 test for specilisation ####
nic1.1 <- glm(CWM.Niche ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(nic1.1)
drop1(nic1.1, test = "F") #-> exclude shading

nic1.2 <- glm(CWM.Niche ~ Comp.2 + max_sumTU_ms_corr + Comp.1, data, family = "gaussian")
summary(nic1.2)
drop1(nic1.2, test = "F") #-> exclude tox

nic1.3 <- glm(CWM.Niche ~ Comp.2 + Comp.1, data, family = "gaussian")
summary(nic1.3)
drop1(nic1.3, test = "F") #-> exclude Comp.1

nic1.4 <- glm(CWM.Niche ~ Comp.2, data, family = "gaussian")
summary(nic1.4)
1 - (0.00091742 / 0.00133139) # D²=0.3109307

summary(glm(CWM.Niche ~ CWM.Ballon1, data, family = "gaussian"))
1 - (0.00095269 / 0.00133139) # D²=0.2844396

summary(glm(CWM.Niche ~ no3, data, family = "gaussian"))
summary(glm(CWM.Niche ~ po4, data, family = "gaussian"))
summary(glm(CWM.Niche ~ totalarea, data, family = "gaussian"))

range(data$Comp.2, na.rm = TRUE)
df <- data.frame(Comp.2 = seq(-1.670808, 3.075343, length.out = 10)) # ,

p_df <- predict(nic1.4,
  newdata = df,
  type = "response",
  se.fit = TRUE
)
tr_df <- transform(df,
  fit = p_df$fit,
  lwr = p_df$fit - (1.96 * p_df$se.fit),
  upr = p_df$fit + (1.96 * p_df$se.fit)
)
nic <- ggplot(tr_df, aes(y = fit, x = Comp.2)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = data, aes(y = CWM.Niche, x = Comp.2), size = 2) +
  labs(
    x = expression("second axis of SPCA"),
    y = expression("CWM for niche width")
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"
  )
nic


#### 4.8 test for web  builder  #####

web1.1 <- glm(CWM.web1 ~ Comp.1 + Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(web1.1)
drop1(web1.1, test = "F") #-> exclude Comp.1

web1.2 <- glm(CWM.web1 ~ Comp.2 + max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(web1.2)
drop1(web1.2, test = "F") #-> exclude Comp.2

web1.3 <- glm(CWM.web1 ~ max_sumTU_ms_corr + shading, data, family = "gaussian")
summary(web1.3)
drop1(web1.3, test = "F") #-> exclude tox

web1.4 <- glm(CWM.web1 ~ shading, data, family = "gaussian")
summary(web1.4)
drop1(web1.4, test = "F")
1 - (0.56494 / 0.87031) # D²=0.350875

range(data$shading, na.rm = TRUE)
df <- data.frame(shading = seq(0, 95, length.out = 10)) # ,

p_df <- predict(web1.4,
  newdata = df,
  type = "response",
  se.fit = TRUE
)
tr_df <- transform(df,
  fit = p_df$fit,
  lwr = p_df$fit - (1.96 * p_df$se.fit),
  upr = p_df$fit + (1.96 * p_df$se.fit)
)
web <- ggplot(tr_df, aes(y = fit, x = shading)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = data, aes(y = CWM.web1, x = shading), size = 2) +
  labs(
    x = expression("shading [%]"),
    y = expression("CWM for web builder")
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"
  ) +
  ylim(-1, 0.2) +
  xlim(NA, 100)
web


#### 4.9 combine plots #####
library(cowplot)
res <- plot_grid(spe, ind, siz, sha, nic, web,
  labels = c("A", "B", "C", "D", "E", "F"),
  align = "v", nrow = 3, ncol = 2
)
res

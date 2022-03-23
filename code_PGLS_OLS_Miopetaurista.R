#Miopetaurista paper

#set directory
setwd("")

#load libraries
library(ape) # read.nexus
library(phytools) #open tree
library(nlme) # GLS analysis
library(geiger) #name.check
library(rr2) # R2

#Open data and tree
squirrel.data<-read.csv("Miopetaurista_data_final_Isaac.csv", header=T,row.names = 1)

#Open tree
squirrel.tree<-read.nexus("Calibrated_squirrels.trees")

#Check if the data and tree have the same names
name.check(squirrel.tree, squirrel.data, data.names=NULL)

#Transform data to log10
squirrel.data$Brain_volume_cm3<-log10(squirrel.data$Brain_volume_cm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_cm3"] <- "Brain"

squirrel.data$Body_mass_g<-log10(squirrel.data$Body_mass_g)
names(squirrel.data)[names(squirrel.data) == "Body_mass_g"] <- "Body"

#plot data to check out the relationship
plot(squirrel.data[, c("Body", "Brain")])

#### PGLS
Lambda_Brain_Body<-gls(Brain ~ Body, data=squirrel.data, 
            correlation=corPagel(value=1,phy=squirrel.tree), method="ML")

Lambda_Brain_Body
summary(Lambda_Brain_Body)

sink("Lambda_Brain_Body.txt")
print(summary(Lambda_Brain_Body))
print(R2_pred(Lambda_Brain_Body))
print(R2_lik(Lambda_Brain_Body))
sink()

#### Now do the predicted and residuals
squirrel.data$Predicted_Brain <- predict(Lambda_Brain_Body)
squirrel.data$Residuals_Brain <- residuals(Lambda_Brain_Body)

#export dataframe
write.csv(squirrel.data,'Squirel_data_residuals.csv')

R2.pred(Lambda_Brain_Body) #better to know how much variation is explaine by the model
R2.resid(Lambda_Brain_Body) #similar to OLS calculation of R2

#### Now do the predicted intervals
intervals(Lambda_Brain_Body)

################# Regression for Extant Sciuridae only

squirrel.data<-read.csv("Miopetaurista_data_final_Isaac.csv", header=T)
names(squirrel.data)[names(squirrel.data) == "X"] <- "Species"

#Open tree
squirrel.tree<-read.nexus("Calibrated_squirrels.trees")

#Transform data to log10
squirrel.data$Brain_volume_cm3<-log10(squirrel.data$Brain_volume_cm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_cm3"] <- "Brain"
squirrel.data$Body_mass_g<-log10(squirrel.data$Body_mass_g)
names(squirrel.data)[names(squirrel.data) == "Body_mass_g"] <- "Body"

#subset dataset per Group
Extant<- squirrel.data[ which(squirrel.data$Family=='Extant'), ]

#Find taxa with specific Family
Extant_sp<-subset(squirrel.data, Family=='Extant', select=Species)

#Select taxa per Family
Extant_species<- unlist (Extant_sp)

#subset tree per Family
Extant_Tree<-drop.tip(squirrel.tree,squirrel.tree$tip.label[-match(Extant_species, squirrel.tree$tip.label)])

#Create model PGLS regression line for Sciuridae
Extant_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Extant_Tree), data=Extant)
summary(Extant_Br_B)

sink("PGLS_Brain_Body_Extant_Br_B.txt")
print(summary(Extant_Br_B))
print(R2_pred(Extant_Br_B))
print(R2_lik(Extant_Br_B))
sink()
#### Now do the predicted intervals
intervals(Lambda_Brain_Body)


#################################### Neocortex #####################################

#Transform data to log10
squirrel.data$Brain_surface_cm2<-log10(squirrel.data$Brain_surface_cm2)
names(squirrel.data)[names(squirrel.data) == "Brain_surface_cm2"] <- "Brain"

squirrel.data$Neocortex_surface_cm2<-log10(squirrel.data$Neocortex_surface_cm2)
names(squirrel.data)[names(squirrel.data) == "Neocortex_surface_cm2"] <- "Neocortex"

#plot data to check out the relationship
par(mfrow=c(1,1))
plot(squirrel.data[, c("Brain","Neocortex",)])

#### PGLS
Lambda_Neocortex<-gls(Neocortex ~ Brain, data=squirrel.data, 
            correlation=corPagel(value=1,phy=squirrel.tree), method="ML")

Lambda_Neocortex
summary(Lambda_Neocortex)

#R2 calculation
R2.pred(Lambda_Neocortex)
R2.resid(Lambda_Neocortex)

#### Now do the predicted and residuals
squirrel.data$Predicted_Neocortex <- predict(Lambda_Neocortex)
squirrel.data$Residuals_Neocortex <- residuals(Lambda_Neocortex)

#export dataframe
write.csv(squirrel.data,'Squirel_data__Neocortex_residuals.csv')

#### Now do the predicted intervals
intervals(Lambda_Brain_Body)

################################## Olfactory bulbs #####################################

#Transform data to log10
squirrel.data$Brain_volume_mm3<-log10(squirrel.data$Brain_volume_mm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_mm3"] <- "Brain"

squirrel.data$Body_mass_mg<-log10(squirrel.data$Body_mass_mg)
names(squirrel.data)[names(squirrel.data) == "Body_mass_mg"] <- "Body"

squirrel.data$Olfactory_bulbs_mm3<-log10(squirrel.data$Olfactory_bulbs_mm3)
names(squirrel.data)[names(squirrel.data) == "Olfactory_bulbs_mm3"] <- "OB"

#plot data to check out the relationship
plot(squirrel.data[, c("Brain", "OB")])
plot(squirrel.data[, c("Body", "OB")])

#### PGLS
###### OB vs. Brain
Lambda_OB_Brain<-gls(OB ~ Brain, data=squirrel.data, 
                     correlation=corPagel(value=1,phy=tree), method="ML")

Lambda_OB_Brain
summary(Lambda_OB_Brain)

#R2 calculation
R2.pred(Lambda_OB_Brain)
R2.resid(Lambda_OB_Brain)

#### Now do the predicted and residuals
squirrel.data$Predicted_OB_Brain <- predict(Lambda_OB_Brain)
squirrel.data$Residuals_OB_Brain <- residuals(Lambda_OB_Brain)

#export dataframe
write.csv(squirrel.data,'Squirel_data_OB_Brain_residuals.csv')

####### OB vs. Body
Lambda_OB_Body<-gls(OB ~ Body, data=squirrel.data, 
                    correlation=corPagel(value=1,phy=tree), method="ML")

Lambda_OB_Body
summary(Lambda_OB_Body)

#R2 calculation
R2.pred(Lambda_OB_Body)
R2.resid(Lambda_OB_Body)

#### Now do the predicted and residuals
squirrel.data$Predicted_OB_Body <- predict(Lambda_OB_Body)
squirrel.data$Residuals_OB_Body <- residuals(Lambda_OB_Body)

#export dataframe
write.csv(squirrel.data,'Squirel_data_OB_Body_residuals.csv')

#### Now do the predicted intervals
intervals(Lambda_Brain_Body)

##################################### Paraflocculus #####################################

#Import data to remove NAs from dataset and tree
data<-read.csv("Miopetaurista_data_final_Isaac.csv", header=T)
data.PL<-data[c(1,5,10)] #select column with taxa for the paraflocculus tree
PL.na<-na.omit(data.PL) #eliminate row with NAs
species<-PL.na$Species # used for making the tree without Petinomys in it

#Import tree to take out NAs
tree<-read.nexus("Calibrated_squirrels.trees")
tree_PL<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
plot(tree_PL, cex=.3)

#Make dataset to match with the tree (without Petinomys)
data1<-read.csv("Miopetaurista_data_final_Isaac.csv", header=T,row.names = 1)
squirrel.data1<-na.omit(data1) #eliminate rows with NAs

#Check tree with data
name.check(tree_PL, squirrel.data1, data.names=NULL)

#Save dataset without Petinomys
write.csv(squirrel.data1,'squirrel_data_paraflocculus.csv')

#Save tree without Petinomys
write.nexus(tree_PL, file = "PL_tree.nex")

#### Analyses for paraflocculus

#Import tree
tree<-read.nexus("PL_tree.nex")

#Open dataset without Petinomys
squirrel.data<-read.csv("squirrel_data_paraflocculus.csv", header=T,row.names = 1)

#Check if the data and tree have the same names
name.check(tree, squirrel.data, data.names=NULL)

#Transform data to log10
squirrel.data$Brain_volume_mm3<-log10(squirrel.data$Brain_volume_mm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_mm3"] <- "Brain"

squirrel.data$Body_mass_mg<-log10(squirrel.data$Body_mass_mg)
names(squirrel.data)[names(squirrel.data) == "Body_mass_mg"] <- "Body"

squirrel.data$Petrosal_lobules_mm3<-log10(squirrel.data$Petrosal_lobules_mm3)
names(squirrel.data)[names(squirrel.data) == "Petrosal_lobules_mm3"] <- "PL"

#plot data to check out the relationship
plot(squirrel.data[, c("Brain", "PL")])
plot(squirrel.data[, c("Body", "PL")])

#### PGLS
#PL vs. Brain
Lambda_PL_Brain<-gls(PL ~ Brain, data=squirrel.data, 
                     correlation=corPagel(value=1,phy=tree), method="ML")

Lambda_PL_Brain
summary(Lambda_PL_Brain)

#R2 calculation
R2.pred(Lambda_PL_Brain)
R2.resid(Lambda_PL_Brain)

#### Now do the predicted and residuals
squirrel.data$Predicted_PL_Brain <- predict(Lambda_PL_Brain)
squirrel.data$Residuals_PL_Brain <- residuals(Lambda_PL_Brain)

#export dataframe
write.csv(squirrel.data,'Squirel_data_PL_Brain_residuals.csv')

#### Now do the predicted intervals
intervals(Lambda_Brain_Body)

#PL vs. Body
Lambda_PL_Body<-gls(PL ~ Body, data=squirrel.data, 
                    correlation=corPagel(value=1,phy=tree), method="ML")

Lambda_PL_Body
summary(Lambda_PL_Body)

#R2 calculation
R2.pred(Lambda_PL_Body)
R2.resid(Lambda_PL_Body)

#### Now do the predicted and residuals
squirrel.data$Predicted_PL_Body <- predict(Lambda_PL_Body)
squirrel.data$Residuals_PL_Body <- residuals(Lambda_PL_Body)

#export dataframe
write.csv(squirrel.data,'Squirel_data_PL_Body_residuals.csv')

#### Now do the predicted intervals
intervals(Lambda_Brain_Body)

### END! :)

Plot the regresion line on the graph, ex: abline(Lambda_Brain_Body) 


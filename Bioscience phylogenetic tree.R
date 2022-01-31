# Code and data for Roff (2020) Evolutionary history drives biogeographic patterns of coral reef resilience

rm(list = ls())

library(ape)
library(ggplot2)
library(phytools)
library(plotrix)


# Load Supertree data
# data for phylogenetic tree from Huang & Roy (2005) via Dryad
# https://royalsocietypublishing.org/doi/10.1098/rstb.2014.0010
# http://dx.doi.org/10.5061/dryad.178n3

supertree <- read.nexus(file="Huang&Roy_Supertree.tre") # load tree
supertreetiplabels<-as.data.frame(supertree[[1]]$tip.label) # extract tip labels
names(supertreetiplabels)<-c('speciesname') # name tip labels

supertreeACR <- supertreetiplabels[grep("ACR_Acr", supertreetiplabels$speciesname), ] # extract Acropora

droplist<-supertreetiplabels %>% filter(!speciesname %in% supertreeACR) # create droplist from species naems
row.names(droplist) <- droplist$speciesname # label rownames
supertree <- pblapply(supertree, drop.tip, tip = as.character(c( as.character(rownames(droplist))))) # drop tips and convert multiPhylo object
class(supertree) <- 'multiPhylo' 

# Load functional trait data from McWilliam et al (2018) PNAS:
# https://www.pnas.org/content/115/12/3084.short

traits <- read.csv("biogeographictraits.csv")
traits$growth2 <- substr(traits$growth,1,9) # subset to combine "open" and "closed" branching growth forms

finaltree <- supertree[[1]] # use first tree from Supertree as exemplar
finaltree$tip.label <- substr(finaltree$tip.label,14,100) # tidy tip.labels
finaltree$tip.label <- gsub("_", " ", finaltree$tip.label) # tidy tip.labels

traits2<-traits[!(traits$species=="ACR_Astreopora_expansa"),] # drop outgroup (ACR_Astreopora_expansa)

traits2$growth2 <- gsub("branching", 1, traits2$growth2) # renumber growthforms
traits2$growth2 <- gsub("corymbose", 2, traits2$growth2)
traits2$growth2 <- gsub("digitate", 3, traits2$growth2)
traits2$growth2 <- gsub("encrustin", 4, traits2$growth2)
traits2$growth2 <- gsub("hispidose", 5, traits2$growth2)
traits2$growth2 <- gsub("tables_or", 6, traits2$growth2)

notrait <-  (rnorm(nrow(traits2),0.15,0)) # vector of trait values for plotTree.wBars {phytools}
attributes(notrait)$names <- traits2$species

numerictrait <-  as.integer(traits2$growth2) # convert renumbered growth forms to integer
attributes(numerictrait)$names <- traits2$species

numerictraitcols <- numerictrait # set colours
numerictraitcols[numerictraitcols==1]<-"#DDCC77"
numerictraitcols[numerictraitcols==2]<-"#117733"
numerictraitcols[numerictraitcols==3]<-"#44AA99"
numerictraitcols[numerictraitcols==4]<-"#88CCEE"
numerictraitcols[numerictraitcols==5]<-"#332288"
numerictraitcols[numerictraitcols==6]<-"#CC6677"

numerictraitlegend <- seq(1:6)  # set legend colours
numerictraitlegend[numerictraitlegend==1]<-"#DDCC77"
numerictraitlegend[numerictraitlegend==2]<-"#117733"
numerictraitlegend[numerictraitlegend==3]<-"#44AA99"
numerictraitlegend[numerictraitlegend==4]<-"#88CCEE"
numerictraitlegend[numerictraitlegend==5]<-"#332288"
numerictraitlegend[numerictraitlegend==6]<-"#CC6677"
names(numerictraitlegend) <- as.list(c("branching", "corymbose", "digitate", # set legend names
                                       "encrusting", "hispidose", "table")) 

# plot phylogenetic tree with label colours by growth form
plotTree.wBars(supertree[[1]],notrait, col=numerictraitcols, scale=10,tip.labels=FALSE,fsize=0.2, width=0.5, type="fan",part=0.95,lwd=1)
obj<-axis(1,pos=-2,at=seq(0,35,by=5),cex.axis=0.5,labels=FALSE) # add species age
text(obj,rep(-6,length(obj)),obj,cex=0.6) # add species age labels
text(mean(obj),-10,"time (ma)",cex=0.8) # add species age legend

for(i in 1:(length(obj)-1)){ # add concentric rings for species age
    a1<-atan(-2/obj[i])
    a2<-0.9*2*pi
    draw.arc(0,0,radius=obj[i],a1,a2,lwd=1,
        col=make.transparent("grey",0.5))
}

add.simmap.legend(colors=numerictraitlegend,x=0.9*par()$usr[1], # add legend
    y=0.9*par()$usr[4],prompt=FALSE,fsize=0.9)


# plot tree with species labels
plotTree.wBars(supertree[[1]],notrait, col=numerictraitcols, scale=10,tip.labels=TRUE,fsize=0.2, width=0.5, type="fan",part=0.95,lwd=1)
obj<-axis(1,pos=-2,at=seq(0,35,by=5),cex.axis=0.5,labels=FALSE)
text(obj,rep(-6,length(obj)),obj,cex=0.6)
text(mean(obj),-10,"time (ma)",cex=0.8)



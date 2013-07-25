library(phybase)
gtv  <- read.delim("genetreebstrees.out",sep="\t",header=F)
lizard.taxaname<-species.name(gtv[1,1])
lizard.species <- lizard.taxaname
species.structure<-matrix(0,48,48)
diag(species.structure)<-1

lizard.steac.trees <- rep("",1000)
for (i in 1:1000){
lizard.steac.trees[i] <- steac.sptree(gtv[,i], speciesname=lizard.species, taxaname=lizard.taxaname, species.structure=species.structure, outgroup="I47", method="nj")
}
write.tree.string(lizard.steac.trees,file="lizard.steac.trees.nex")

DNAcon<-read.csv("data/metadataDNA.conc.2021.2022.csv")

DNAcon21<-DNAcon[(DNAcon$Year=="2021"),]
DNAcon22<-DNAcon[(DNAcon$Year=="2022"),]

boxplot(DNAcon$Qubit.DNA..ng.ul~DNAcon$Organism, ylab="DNA (ng/ul)", xlab="Species")
boxplot(DNAcon21$Qubit.DNA..ng.ul~DNAcon21$Organism, ylab="DNA (ng/ul)", xlab="Species")
boxplot(DNAcon22$Qubit.DNA..ng.ul~DNAcon22$Organism, ylab="DNA (ng/ul)", xlab="Species")

calanoid<-DNAcon[(DNAcon$Organism=="Calanoid"),]
boxplot(calanoid$Qubit.DNA..ng.ul~calanoid$Number.of.organisms, 
        ylab="DNA (ng/ul)", xlab="# of organisms", main="Calanoids")

Daphnia<-DNAcon[(DNAcon$Organism=="Daphnia"),]
boxplot(Daphnia$Qubit.DNA..ng.ul~Daphnia$Number.of.organisms, 
        ylab="DNA (ng/ul)", xlab="# of organisms", main= "Daphnia")

Ceriodaphnia<-DNAcon[(DNAcon$Organism=="Ceriodaphnia"),]
boxplot(Ceriodaphnia$Qubit.DNA..ng.ul~Ceriodaphnia$Number.of.organisms, 
        ylab="DNA (ng/ul)", xlab="# of organisms", main= "Ceriodaphnia")


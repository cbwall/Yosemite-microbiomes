getwd()
IDs<-read.csv("Davis.IDs.csv")
head(IDs)

# well combined from rows and columns
IDs$well<-paste(IDs$row,IDs$column, sep= "")

# make UCDavis ID as required
IDs$final.ID<-interaction(IDs$submitter, IDs$UCD.number, IDs$gene, sep="_")

write.csv(IDs, "Davis.IDs.csv")

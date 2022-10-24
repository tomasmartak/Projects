### 1. Libraries
library(reshape2)
library(dplyr)
library(eulerr)
library(ggplot2)
library(RVenn)

### 2. Dataset Loading
MSdata <- as_tibble(read.csv2("raw/MSdata.csv", sep = ",")) #make formatting nice and usable
#class(MSdata)
#head(MSdata)
#colnames(MSdata)
#View(MSdata)

### 3. Filtering (based on the criteria that I received)
# Replace NaNs with 0. Replace anything under 21 (not including) with 0 --> data = irrelevant, as the intensity is too low.
# Also removes all proteins where intensity = NaN/0/<21 across the board.
MSdata[is.na(MSdata)] <- 0           
for (i in seq_along(MSdata[,1:5])) {
  MSdata[[i]][MSdata[[i]] <21] <- 0 #replace everything <21 with 0
}
rm(i) # I don't want i stored. Then again I don't know how else to remove it, and don't care to find out how for now, so here goes.

MSdata <- MSdata[rowSums(MSdata[,1:5])>0,] #removes all proteins where values were 0/NaN
MSdata <- MSdata[(MSdata$unique.pepts)>=2,] #removes all proteins where numbers of unique peptides measured > 2.
for (i in seq_along(MSdata[,1:5])) {
  MSdata[[i]][MSdata[[i]] >=21] <- 1 #replace everything >=21 with 1; intensities are qualitative, not quantitative. We want a simple True/False/NA
}
rm(i) #ditto

### Remove anything after 1st semicolon in MSdata$name.gene/prot --> we want only 1 gene/protein name
MSdata$name.gene <- gsub("\\;.*","",MSdata$name.gene)
MSdata$name.prot <- gsub("\\;.*","",MSdata$name.prot)
#head(MSdata)
#tidy up data into a basic useable dataframe
MS.tidy <- MSdata[,c(8,3,4,5,1,2)] #orders in tidiest format based on protein name
#View(MS.tidy)
#MS.tidy <- dplyr::rename(MS.tidy, protein = name.prot)
MS.tidy <- dplyr::rename(MS.tidy, protein = name.prot)
MS.tidy <- as.data.frame(MS.tidy) #convert to dataframe only, as reshape can't handle tibbles; just a dplyr bug
#head(MS.tidy)
#View(MS.tidy) #opens in new tab for ease of access

### ***VENN DIAGRAM***
#we start with MS.heat --> it's a good format; Protein names are irrelevant, so I removed them for MS.venn
MS.venn = MS.tidy[,2:6] #can't include protein names --> this bit will have to be done manually
#head(MS.venn)
MS.venn <- venn(MS.venn) #here you can use either venn or euler as a command; this will change the shape of the plot;
                         #euler = based on size; venn = based on looking nice; also, can add shape = "ellipse" for more faithful 
                         #depiction (venn misses the X1/X4 only overlap, and solo X5) 
#head(MS.tidy)
#MS.venn shows relevant data
plot(MS.venn,
     quantities = T,
     lty = 1:3,
     labels = list(font = 4))

MS.venn.values <- MS.venn$original.values
MS.venn.values
print(MS.venn)
#getting actual venn values; for this, MS.heat is ideal if  a subset is made based on "value" = 1 and removing the respective column.
MS.mini <- subset(MS.heat, value==1)[,1:2][,c(2,1)]
head(MS.mini)
#View(MS.mini)
#convert MS.mini into a list of lists, with each list named after a sample (X1:5) and including all proteins with intensity >21
l1 <- list(subset(MS.mini, sample=="X1")$protein)
l2 <- list(subset(MS.mini, sample=="X2")$protein)
l3 <- list(subset(MS.mini, sample=="X3")$protein)
l4 <- list(subset(MS.mini, sample=="X4")$protein)
l5 <- list(subset(MS.mini, sample=="X5")$protein)
MS.list <- c(l1,l2,l3,l4,l5)
names(MS.list) <- c("X1", "X2","X3","X4", "X5") #sets names for list; not strictly necessary, but makes indexing more legible
#summary(MS.list)

### convert MS.list into an RVenn object and do the blasted combinations
# Is there a better way to do this? Yes, there is. Can I quickly find it? No. Am I going to use this script again, though? Absolutely not.
# Next time I'm actually going to follow a proper MS pipeline, so I'm not going to spend more time optimalising this disgusting stuff.
# Why did I even do this? This was the output I was told to make. Do I personally think it is very useful? Nope.
MS.RVenn <- Venn(MS.list)
summary(MS.RVenn)
a <- list(MS.list$X1[!(MS.list$X1 %in% RVenn::unite(MS.RVenn, c(2:5)))])
b <- list(MS.list$X2[!(MS.list$X2 %in% RVenn::unite(MS.RVenn, c(1,3:5)))])
c <- list(MS.list$X3[!(MS.list$X3 %in% RVenn::unite(MS.RVenn, c(1,2,4,5)))])
d <- list(MS.list$X4[!(MS.list$X4 %in% RVenn::unite(MS.RVenn, c(1:3,5)))])
e <- list(MS.list$X5[!(MS.list$X5 %in% RVenn::unite(MS.RVenn, c(1:4)))])
f <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,5)), RVenn::unite(MS.RVenn, c(2:4))))
g <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,4)), RVenn::unite(MS.RVenn, c(2,3,5))))
h <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,2)), RVenn::unite(MS.RVenn, c(3:5))))
i <- list(setdiff(RVenn::overlap(MS.RVenn, c(2,5)), RVenn::unite(MS.RVenn, c(1,3,4))))
j <- list(setdiff(RVenn::overlap(MS.RVenn, c(2,3)), RVenn::unite(MS.RVenn, c(1,4,5))))
k <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,3)), RVenn::unite(MS.RVenn, c(2,4,5))))
l <- list(setdiff(RVenn::overlap(MS.RVenn, c(3,4)), RVenn::unite(MS.RVenn, c(1,2,5))))
m <- list(setdiff(RVenn::overlap(MS.RVenn, c(4,5)), RVenn::unite(MS.RVenn, c(1:3))))
n <- list(setdiff(RVenn::overlap(MS.RVenn, c(3,5)), RVenn::unite(MS.RVenn, c(1,2,4))))
o <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,2,4)), RVenn::unite(MS.RVenn, c(3,5))))
p <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,3,4)), RVenn::unite(MS.RVenn, c(2,5))))
q <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,3,4,5)), MS.list$X2))
r <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,4,5)), RVenn::unite(MS.RVenn, c(2,3))))
s <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,2,4,5)), MS.list$X3))
t <- list(setdiff(RVenn::overlap(MS.RVenn, c(1,2,5)), RVenn::unite(MS.RVenn, c(3,4))))
u <- list(setdiff(RVenn::overlap(MS.RVenn, c(1:3,5)), MS.list$X4))
v <- list(setdiff(RVenn::overlap(MS.RVenn, c(1:3)), RVenn::unite(MS.RVenn, c(4,5))))
w <- list(setdiff(RVenn::overlap(MS.RVenn, c(1:4)), MS.list$X5))
x <- list(setdiff(RVenn::overlap(MS.RVenn, c(2:4)), RVenn::unite(MS.RVenn, c(1,5))))
y <- list(setdiff(RVenn::overlap(MS.RVenn, c(2:5)), MS.list$X1))
z <- list(setdiff(RVenn::overlap(MS.RVenn, c(3:5)), RVenn::unite(MS.RVenn, c(1,2))))
alpha <- list(RVenn::overlap(MS.RVenn))
Venn.summary <- as.list(c(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,alpha))
names(Venn.summary) <- c(letters, "alpha") #set names for easier legibility
#View(Venn.summary)

#nicely create a csv. file with the Venn.summary; this can be looked at along with the Venn diagram
for(i in 1:27){
  write.table(Venn.summary[i], 'Annotation-proteins.csv', append= T, sep=',', row.names = TRUE)
}



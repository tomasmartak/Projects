### 0. Information
# This script was created to enable simple web scraping from Gene Ontology datasets.
# The dummy link here links to p53 expression from cardiomyopathy and control human septum biopsy gene arrays.
# The goal of this exercise was to create a simple script backbone to use in RStudio for different genes.
      # Note: now that I know more about image recognition, more is possible. However I haven't had the time or neccessity to do so,
      #       so I just manually lable things for now. 

### 1. Libraries
library(ggplot2)
library(ggpubr)
library(tibble) #absolutely unnecessary, but then you won't have that lovely tibble. Look it in the eyes and tell it you don't want it. C'mon. Thought so.
library(xml2) #table scraper from online html
library(rvest) #table scraper from online html

### 2. Functions
#create recurrence function; I was too lazy to learn to use Roxygen2 (might do later if I feel cute)
#What does the function do? Easy. You paste in the Title entry corresponding of the row that includes the 1st recurrent cancer data.
recurrence <- function(divider){
  i = 1
  while(df$Title[i] != divider){
    df$Title[i] <- "non-recurrent"
    i = i + 1
  }
  for(a in i:nrow(df)){
    if(df$Title[a] != "non-recurrent"){
      df$Title[a] <- "recurrent"
    }
  }
  df$Title <<- df$Title #global assignment operator; very useful if you want the data created to actually leave this "nice" function.
} #see info above; this function is nice.

### 3. Define the URL of the GEO graphs webpage you want to analyse
webpage_url <- "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS2206:IMAGp998J225591" #enter webpage here; this is dummy data for p53 expression
first_recurrence <- "P72 dilated cardiomyopathy heart biopsy"#enter 1st recurrence threshold value for df$Title

### 4. Run this code to get a graph.
webpage <- xml2::read_html(webpage_url) #reads url
df <- rvest::html_table(webpage)[[2]][,1:3] %>% tibble::as_tibble() %>% tidyr::drop_na() # the index here is [[2]]; [,1:3] removes useless af rank
graph_title <- toString(rvest::html_table(webpage)[[1]][2,2]) # saves the graph title
dataset_name <- gsub("/","-",rvest::html_table(webpage)[[1]][1,2]) # saves name of dataset for ease of saving the file
recurrence(first_recurrence) 
ggplot(data = df, aes(x = Title, y = Value))+
  geom_bracket(data = df, xmin = "non-failing", xmax = "dilated cardiomyopathy", label = "", size = 0.75, y.position = 1.05*max(df$Value))+
  geom_boxplot(aes(fill = Title)) +
  geom_jitter(shape=16, width = 0.1) +
  expand_limits(y = 1.07*max(df$Value)) + #just in case the lable is outside the scope
  stat_compare_means(aes(label = paste0(..method.., ", ", "p =", ..p.format..)), 
                     method = "t.test", label.x = 1.35, label.y = 1.07*max(df$Value)) + #put ..p.format.. as label for p value
  labs(title = graph_title, x = "Disease state", y = "log2 ratio")

paste0(dataset_name,": ",graph_title) # done to get quick information on the dataset in question
#aes(label = ..p.signif..) <- for less cluttered lable for stat_compare_means
#vjust = -2 <- for being lazy with stat_compare_means if you don't want to bugger around with x and y coordinates


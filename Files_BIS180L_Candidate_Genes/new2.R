library(tidyverse)
GO <- read_csv("GO_Terms.txt", col_names = F)
Trans <- read_table2("transcript_names.txt", col_names = F)
Trans$X1 <- sub("::.*","",Trans$X1)
Trans$X2 <- sub("UniRef90_","",Trans$X2)

New <- inner_join(Trans,GO, by = c("X2" = "X1"))
New <- New[,-2]

write_csv(New, "Genes_and_GO_Terms.tab", col_names = F)
system("sed 's/NA//g' Genes_and_GO_Terms.tab | 
        sed 's/,*$//' |
        sed 's/,/\t/' > Genes_and_GO_Terms.tab")

sed -e "s/ //g" uniref_GO_terms.tab | #Removing all spaces
sed -e "s/;/,/g" | sed -e "s/^.*	//" | #Change semicolons to commas
sed -e "s/^[a-zA-Z,]*\[/\[/" | #Remove all characters from start of line until "[" is reached
sed -e "s/,[a-zA-Z,-]*\[/\[/g" | # Remove all characters after a comma, until "["
sed -e "s/^[a-zA-Z,-]*\[/\[/" > bad
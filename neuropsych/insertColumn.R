#insert a new column into your dataframe

#inputs: existing DF
#        new Column
#        c = the index of the column


insertColumn <- function(existingDF, newcolumn, c) {
  existingDF[,seq(c+1,ncol(existingDF)+1)] <- existingDF[,seq(c,ncol(existingDF))]
  existingDF[,c] <- newcolumn
  existingDF
}

#taking the same ideas, but applying as a column instead of a row
#Reference -- https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended
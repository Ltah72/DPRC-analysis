#insert a new row into your dataframe

#inputs: existing DF
#        new Row
#        r = the index of the row


insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}


#Reference -- https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended
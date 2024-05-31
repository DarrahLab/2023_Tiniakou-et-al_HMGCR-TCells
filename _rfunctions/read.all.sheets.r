  read.all.sheets <- function(filename, tibble = FALSE) {
  #Author: Jeromy Anglim 
  #Stack overflow https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames
  #Requires readxl library
  #His commentary below. 
  # I added suppressMessages to prevent new names from printing to console.  AAG 16Feb2023 
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <-  suppressMessages(readxl::excel_sheets(filename))
  x <- suppressMessages(lapply(sheets, function(X) readxl::read_excel(filename, sheet = X)))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}
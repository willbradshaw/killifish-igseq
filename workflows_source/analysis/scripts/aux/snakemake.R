#########################################################
## AUXILIARY R FUNCTIONS FOR MANIPULATING SNAKEMAKE VARIABLES ##
#########################################################

#-----------------------------------------------------------------------------
# I/O
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Manipulating range TSV tables by label
#-----------------------------------------------------------------------------

filter_by_label <- function(range_table, pattern, complement=FALSE){
  matches <- grepl(pattern, range_table[["label"]])
  if (complement) matches <- !matches
  return(filter(range_table, matches))
}

substitute_labels <- function(range_table, sub_patterns, sub_replacements){
  multisub(range_table, "label", label_patterns, label_replacements)
}


#-----------------------------------------------------------------------------
# Generalised substitution methods
#-----------------------------------------------------------------------------

sub_table <- function(table, column, pattern, replacement){
  # Perform a pattern substitution on a column of a tibble or dataframe
  table[[column]] <- gsub(pattern, replacement, table[[column]])
  return(table)
}

multisub <- function(table, column, patterns, replacements){
  # Perform multiple substitutions on a table column based on lists of patterns
  for (n in seq(length(patterns))){
    table <- sub_table(table, column, patterns[[n]], replacements[[n]])
  }
  return(table)
}
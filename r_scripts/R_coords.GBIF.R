coords.GBIF <- function(spvector) 
{
  print(" ... acquiring GPS coordinates from GBIF, please wait")
  keys <- pbsapply(spvector, function(x) name_suggest(x, rank="species")$key[1])
  keys <- keys[!sapply(keys, is.null)] 

  ls  <- list()
  for (key in keys) 
  {
    tryCatch(
    {
      tmp  <- occ_search(taxonKey=key, hasCoordinate=TRUE, limit=1000, 
                        fields=c("name", "key", "decimalLatitude", "decimalLongitude", 
                                 "country", "elevation"), return='data')

      name <- tmp$name[1]
      ls[[name]] <- tmp
      print(paste0(name, " ... done! - rows: ", length(tmp$name)))
    }, 
    error=function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  }

  return(rbindlist(lapply(ls, as.data.frame.list), fill=TRUE, use.names = TRUE))
}
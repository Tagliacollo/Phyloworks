library(pbapply)
library(rfishbase)
library(data.table)

chrs.FishBase <- function(spvector)
{
  print("... acquiring species data from FishBase, please wait")
  keys <- pbsapply(spvector, function(x) validate_names(x))
  keys <- keys[!sapply(keys, is.null)] 

  ls <- list()
  for (key in keys)
  {
 tryCatch(
    {
      name  <- key
      main  <- species(spvector, fields = c("sciname", "DemersPelag", "Vulnerability", 
                                            "Length", "PD50", "Comments"))
      ecol  <- ecology(spvector, fields = c("sciname", "FoodTroph", "FoodSeTroph", 
                                            "FoodRemark", "AddRems", "Herbivory2"))
      stock <- stocks(spvector, fields  = c("sciname", "EnvTemp", "Remark", 
                                            "ResilienceRemark"))
      tmp   <- Reduce(function(...) merge(..., all=TRUE), list(main, ecol, stock))      

      ls[[name]] <- tmp
      print(paste0(name, " ... done"))
    }, 
    error=function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  }

  return(rbindlist(lapply(ls, as.data.frame.list), fill=TRUE, use.names = TRUE))
}  

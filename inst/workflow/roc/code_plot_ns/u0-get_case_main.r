### Convert i.case to case.main for plotting title.

get.case.main <- function(i.case, model){
  ### cases with phi.
  if(i.case == paste(model, "_ad_wphi_pm", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: posterior mean", sep = "")
  } else if(i.case == paste(model, "_ad_wphi_scuo", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: SCUO", sep = "")
  } else if(i.case == paste(model, "_ad_wphi_true", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: True", sep = "")
  } else if(i.case == paste(model, "_ad_wphi_bInit", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: bInit", sep = "")

  ### cases without phi.
  } else if(i.case == paste(model, "_ad_wophi_pm", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: posterior mean", sep = "")
  } else if(i.case == paste(model, "_ad_wophi_scuo", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: SCUO", sep = "")
  } else if(i.case == paste(model, "_ad_wophi_true", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: True", sep = "")
  } else if(i.case == paste(model, "_ad_wophi_bInit", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: bInit", sep = "")

  ### cases without phi following with phi.
  } else if(i.case == paste(model, "_ad_wphi_wophi_pm", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: wphi,posterior mean", sep = "")
  } else if(i.case == paste(model, "_ad_wphi_wophi_scuo", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: wphi,SCUO", sep = "")

  } else{
    i.case.main <- "case.name Not Found"
  }
  i.case.main
}

new.page <- function(workflow.name, i.case, model, default.layout = TRUE){
  if(default.layout){
    nf <- layout(matrix(c(1, 1, 2, 3, 4, 5, 6, 7),
                        nrow = 4, ncol = 2, byrow = TRUE),
                 c(1, 1), c(1, 8, 8, 8), respect = FALSE)
  }

  # Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5,
       paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  invisible()
} # End of new.page().

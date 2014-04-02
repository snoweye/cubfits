# Propose new B conditioning on all orhter parameters.
#
# These functions are for all amino acids.

get.my.drawBConditionalAll <- function(type){
  if(!any(type[1] %in% .CF.CT$init.fit)){
    stop("fit is not found.")
  }
  ret <- eval(parse(text = paste("my.drawBConditionalAll.",
                                 type[1], sep = "")))
  assign("my.drawBConditionalAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.drawBConditionalAll().


# Draw new B.
my.drawBConditionalAll.current <- function(bCurr, phi.Curr, y, n, reu13.df.obs){
  ### Note that phi.new = phi.Curr is the E[Phi] rather than phi.Obs.
  bFit <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Curr, y, n,
                                        phi.new = phi.Curr, coefstart = bCurr)

  ### Based on the above new fits of parameters to draw new beta (M, S_1, S_2).
  ret <- lapply(1:length(reu13.df.obs),
           function(i.aa){ # i'th amino acid.
             tmp <- .cubfitsEnv$my.drawBConditionalFit(
                      bFit[[i.aa]], bCurr[[i.aa]], phi.Curr, y[[i.aa]], n[[i.aa]],
                      reu13.df.aa = reu13.df.obs[[i.aa]])
             # update S/M's acceptance.
             my.update.acceptance("B", tmp$accept, i.aa)
             tmp
           })
  ret
} # End of my.drawBConditionalAll.current().

my.drawBConditionalAll.random <- function(bCurr, phi.Curr, y, n, reu13.df.obs){
  ### Note that phi.new = phi.Curr is the E[Phi] rather than phi.Obs.
  bFit <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Curr, y, n,
                                        phi.new = phi.Curr)

  ### Based on the above new fits of parameters to draw new beta (M, S_1, S_2).
  ret <- lapply(1:length(reu13.df.obs),
           function(i.aa){ # i'th amino acid.
             tmp <- .cubfitsEnv$my.drawBConditionalFit(
                      bFit[[i.aa]], bCurr[[i.aa]], phi.Curr, y[[i.aa]], n[[i.aa]],
                      reu13.df.aa = reu13.df.obs[[i.aa]])
             # update S/M's acceptance.
             my.update.acceptance("B", tmp$accept, i.aa)
             tmp
           })
  ret
} # End of my.drawBConditionalAll.random().

### Get the specific function according to the options.
get.my.update.DrawScale <- function(adaptive){
  if(!any(adaptive[1] %in% .CF.CT$adaptive)){
    stop("adaptive is not found.")
  }
  ret <- eval(parse(text = paste("my.update.DrawScale.",
                                 adaptive[1], sep = "")))
  assign("my.update.DrawScale", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.update.DrawScale().


### No adaptive.
my.update.DrawScale.none <- function(var.name, phi.DrawScale,
    update.curr.renew = TRUE){
  phi.DrawScale
} # End of my.update.DrawScale.none().


### Update scaling factors for every gene.
my.update.DrawScale.simple <- function(var.name, phi.DrawScale,
    update.curr.renew = TRUE){
  ### Update new scaling factors.
  ret <- my.DrawScale.scaling(var.name, .cubfitsEnv$curr.renew)

  ### Update curr.renew to global.
  .cubfitsEnv$DrawScale[[var.name]][[.cubfitsEnv$curr.renew + 1]] <- ret

  if(update.curr.renew){
    .cubfitsEnv$curr.renew <- .cubfitsEnv$curr.renew + 1
  }

  ret
} # End of my.update.DrawScale().

my.DrawScale.scaling <- function(var.name, curr.window){
  if(curr.window > 2){
    prev.scale <- .cubfitsEnv$DrawScale[[var.name]][[curr.window - 1]]
    prev.accept <- .cubfitsEnv$adaptive[[var.name]][[curr.window - 1]] /
                   .CF.AC$renew.iter
  }
  curr.scale <- .cubfitsEnv$DrawScale[[var.name]][[curr.window]]
  curr.accept <- .cubfitsEnv$adaptive[[var.name]][[curr.window]] /
                 .CF.AC$renew.iter
  ret <- curr.scale

  ### Smaller than the target.lower.
  id <- which(curr.accept <= .CF.AC$target.accept.lower)
  if(length(id) > 0){
    tmp <- .CF.AC$scale.decrease
    if(curr.window > 2){
      tmp <- lapply(id, function(i){
                 if((curr.accept[i] < prev.accept[i] &&
                     curr.scale[i] < prev.scale[i]) ||
                    (curr.accept[i] > prev.accept[i] &&
                     curr.scale[i] > prev.scale[i])){
                    .CF.AC$scale.increase
                 } else{
                    .CF.AC$scale.decrease
                 }
               })
      tmp <- do.call("c", tmp)
    }
    ret[id] <- curr.scale[id] * tmp
  }

  ### Larger than the target.upper.
  id <- which(curr.accept >= .CF.AC$target.accept.upper)
  if(length(id) > 0){
    tmp <- .CF.AC$scale.increase
    if(curr.window > 2){
      tmp <- lapply(id, function(i){
                 if((curr.accept[i] < prev.accept[i] &&
                     curr.scale[i] < prev.scale[i]) ||
                    (curr.accept[i] > prev.accept[i] &&
                     curr.scale[i] > prev.scale[i])){
                    .CF.AC$scale.decrease
                 } else{
                    .CF.AC$scale.increase
                 }
               })
      tmp <- do.call("c", tmp)
    }
    ret[id] <- curr.scale[id] * tmp
  }

  ### Replace too small and too large numbers.
  ret[ret > .CF.AC$sigma2.upper] <- .CF.AC$sigma2.upper
  ret[ret < .CF.AC$sigma2.lower] <- .CF.AC$sigma2.lower

  ret
} # End of my.DrawScale.scaling().

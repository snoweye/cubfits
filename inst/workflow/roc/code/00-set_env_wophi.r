### This is a basic configuration for workflow "without phi".

### Specify the case name.
workflow.name <- "wophi"

### Specify model.
model <- "roc"
file.data.fasta <- "genome.fasta"
file.data.tsv <- ""

### For default path.
prefix <- list()
prefix$root <- "./"
prefix$param <- paste(prefix$root, "param/", sep = "")

prefix$all.out <- "./all.out/"
prefix$data <- paste(prefix$all.out, "data/", sep = "")
prefix$subset <- paste(prefix$all.out, "subset/", sep = "")
prefix$output <- paste(prefix$all.out, "output/", sep = "")
prefix$table <- paste(prefix$all.out, "table/", sep = "")

### For ploting.
prefix$plot <- paste(prefix$all.out, "plot/", sep = "")
prefix$plot.diag <- paste(prefix$all.out, "plot/diag/", sep = "")
prefix$plot.match <- paste(prefix$all.out, "plot/match/", sep = "")
prefix$plot.single <- paste(prefix$all.out, "plot/single/", sep = "")
prefix$plot.trace <- paste(prefix$all.out, "plot/trace/", sep = "")
prefix$plot.AA <- paste(prefix$all.out, "plot/AA/", sep = "")

### For ploting without scaling.
prefix$plot.nps <- paste(prefix$all.out, "plot_nps/", sep = "")
prefix$plot.nps.diag <- paste(prefix$all.out, "plot_nps/diag/", sep = "")
prefix$plot.nps.match <- paste(prefix$all.out, "plot_nps/match/", sep = "")
prefix$plot.nps.single <- paste(prefix$all.out, "plot_nps/single/", sep = "")
prefix$plot.nps.trace <- paste(prefix$all.out, "plot_nps/trace/", sep = "")
prefix$plot.nps.AA <- paste(prefix$all.out, "plot_nps/AA/", sep = "")

### For code.
prefix$code <- paste(cubfits::get.workflow(model = model), "/",
                     "code/", sep = "")
prefix$code.plot <- paste(cubfits::get.workflow(model = model), "/",
                          "code_plot/", sep = "")
prefix$code.plot.nps <- paste(cubfits::get.workflow(model = model), "/",
                             "code_plot_nps/", sep = "")

### Specify data files.
file.data <- list()
file.data$fasta <- paste(prefix$param, file.data.fasta, sep = "")
if(exists("file.data.tsv") && file.data.tsv != ""){
  file.data$tsv <- paste(prefix$param, file.data.tsv, sep = "")
} else{
  file.data$tsv <- ""
}

### All case names.
case.names <- c("ad_wophi_scuo")
case.names <- paste(model, "_", case.names, sep = "")

### Basic information.
run.info <- list()
run.info$nIter <- 4000
run.info$burnin <- 1000
run.info$phi.DrawScale <- 1

### For configuration.
run.info$dump <- FALSE
run.info$prefix.dump <- paste(prefix$output, "tmp/dump_", sep = "")
run.info$parallel <- "task.pull"
run.info$adaptive <- "simple"

### For MCMC.
range <- list()
range$subset <- 3001:5000
range$thinning <- 10

### For simulation only.
simulation <- list()
simulation$EPhi <- FALSE
simulation$Eb <- FALSE
simulation$seed <- 1234
simulation$sdlog <- 1.5

### For plotting.
ci.prob <- c(0.025, 0.975)

### For modeling of logmixture.
p.nclass <- 2


### CAUTION: for extra changes globally.
# suppressMessages(library(cubfits, quietly = TRUE))
# .CF.CT$init.fit <- "RW_Norm"
# .CF.CT$type.p <- "lognormal_RW"
# .CF.CT$model.Phi <- "logmixture"

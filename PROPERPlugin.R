library(PROPER)
library(ggplot2)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


input <- function(inputfile) {
        pfix = prefix()
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  print(parameters)
   # Need to get the three files
   #csvfile <<- paste(pfix, parameters["csvfile", 2], sep="/")

   #myData <<- read.csv(csvfile)
   #mdrrClass <<- readLines(paste(pfix, parameters["classes", 2], sep="/"))

}

run <- function() {

}

output <- function(outputfile) {

sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = as.integer(parameters["ngenes", 2]), p.DE=as.numeric(parameters["de", 2]),
      lOD="cheung", lBaselineExpr="cheung")

## ----echo=TRUE,eval=FALSE--------------------------------------------------
    sim.opts.Bottomly = RNAseq.SimOptions.2grp(as.integer(parameters["ngenes", 2]), p.DE=as.numeric(parameters["de", 2]),
      lOD="bottomly", lBaselineExpr="bottomly")

## ----echo=TRUE,eval=FALSE,result=FALSE-------------------------------------
    simres = runSims(Nreps = c(3, 5, 7, 10), sim.opts=sim.opts.Cheung,
      DEmethod="edgeR", nsims=20)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
    powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
      stratify.by="expr", delta=0.5)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
    summaryPower(powers)

## ----eval=FALSE,echo=TRUE--------------------------------------------------
    plotPower(powers)

## ----eval=FALSE,echo=TRUE--------------------------------------------------
    plotPowerTD(powers)

## ----eval=FALSE,echo=TRUE--------------------------------------------------
    plotFDcost(powers)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
    plotAll(powers)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
    power.seqDepth(simres, powers)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
    powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
      strata = c(0, 10, 2^(1:7)*10, Inf), filter.by="expr",
      strata.filtered=1, stratify.by="expr", delta=0.5)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
      powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
        stratify.by="dispersion", target.by="effectsize", delta=1)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
      powers = comparePower(simres, alpha.type="pval", alpha.nominal=0.001,
        stratify.by="dispersion", target.by="effectsize", delta=1)

      ggsave(outputfile)

## ----echo=TRUE, result=TRUE------------------------------------------------
}

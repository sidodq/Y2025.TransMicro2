#' @title Perform CCA abundance table
#' @param memberships membership of rows (samples) returned by sbm::estimateBipartiteSBM
#' @description
#'Perform CCA using functions from Ravigné et al. 2022 (based on Massol et al 2021)
#'
#'
#' @return
#' XX
#'
#' @examples
#' wholenetwork_cca_results_NULL_ITS <- analysis.function.alt2(incid=as.data.frame(FB_ITS), memb_obs=memb_spl_obs_ITS, var1=cca_data_ITS$FruitType, var2=cca_data_ITS$LifeStage, configs = NULL)
#'
#' @export

analysis.function.alt2 <- function(incid,memb_obs, var1, var2, configs = NULL){#industrial processing of CCA, two permutations
  n <- dim(incid)[1]
  modules <- lme4::dummy(as.factor(memb_obs), levelsToKeep = levels(as.factor(memb_obs)))

  one.1 <- vegan::cca(modules ~ var1, na.action=na.omit)
  one.2 <- vegan::cca(modules ~ var2, na.action=na.omit)
  two.1 <- vegan::cca(modules ~ var1 + var2, na.action=na.omit)
  one.1.alone <- vegan::cca(modules ~ var1 + Condition(var2), na.action=na.omit)
  one.2.alone <- vegan::cca(modules ~ Condition(var1) + var2, na.action=na.omit)

  xxx <- list(one.1,one.2,two.1,one.1.alone,one.2.alone)

  chi2 <- sapply(1:5, function(x) chi2.computation(xxx[[x]]))
  Forms <- sapply(1:5, function(x) extractTerms(xxx[[x]]))
  F <- sapply(1:5, function(x) F.computation(xxx[[x]],n))
  P <- sapply(1:5, function(x) anova(xxx[[x]],permutations = permute::how(nperm=NPERM-1))$P[1])

  if(is.null(configs)) {
    results<-data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F,"P"=P)
  } else {
    depth <- dim(configs)[2]
    all.F <- sapply(1:depth,function(x)
      only.F.alt2(as.vector(configs[,x]),var1,var2))
    #only.F.alt2(as.vector(configs[,x])[[1]],var1,var2))
    all.ecdf <- apply(all.F,1,ecdf)
    P2 <- sapply(1:5,function(x) 1-((all.ecdf[[x]])(F[x])))
    results <- data.frame("Formulas"=Forms,"chi2"=chi2,"F"=F,"P"=P,"P2"=P2)
  }
  results
}

only.F.alt2 <- function(memberships, var1, var2) {
  n <- length(memberships)

  modules <-
    lme4::dummy(as.factor(memberships), levelsToKeep = levels(as.factor(memberships)))

  one.1 <- vegan::cca(modules ~ var1, na.action = na.omit)
  one.2 <- vegan::cca(modules ~ var2, na.action = na.omit)
  two.1 <- vegan::cca(modules ~ var1 + var2, na.action = na.omit)
  one.1.alone <-
    vegan::cca(modules ~ var1 + Condition(var2), na.action = na.omit)
  one.2.alone <-
    vegan::cca(modules ~ Condition(var1) + var2, na.action = na.omit)

  xxx <-
    list(
      one.1,
      one.2,
      two.1,
      one.1.alone,
      one.2.alone)

  sapply(1:5, function(x)
    F.computation(xxx[[x]], n))
}

chi2.computation<-function(ccaobj){
  ccaobj$CCA$tot.chi
}

extractTerms<-function(model){
  paste0(formula(model$terms))[3]
}

F.computation<-function(ccaobj,size){
  if(is.null(ccaobj)){
    0
  }
  else{
    ccaobj$CCA$tot.chi*(size-ccaobj$CA$rank-nulltozero(ccaobj$pCCA$rank))/(ccaobj$CA$tot.chi*ccaobj$CCA$qrank)
  }
}

nulltozero<-function(x){
  if(is.null(x)){
    x<-0
  }
  x
}

#Post Simulation Analysis

loadlib <- function(){
  library(stringr) #rational string functions
  library(plyr)  #rational wrapper for apply like functions
  library(rlist)
  library(utils)
  library(moments)
  library(boot)
  library(plotrix)
}

loadlib()

fltn_ev_bnd_lst <- function(eblnm,cntrst) # flatten an ev_aM_clbrtd list object
{
  #
  ebl=get(eblnm)
  
  if (str_detect(string = eblnm,pattern = "big")) ss <- 200
  if (str_detect(string = eblnm,pattern = "small")) ss <- 50
  if (str_detect(string = eblnm,pattern = "tiny")) ss <- 20
  if (str_detect(string = eblnm,pattern = "huge")) ss <- 2000
  rws <- length(ebl)
  num_bnds <- dim(ebl[[1]]$ev_summry$ev_bnds)[1]
  trials=data.frame(trial=1:rws)
  ev0 <- rep(NA,rws)
  mn_ev <- rep(NA,rws)
  md_ev <- rep(NA,rws)
  bnds <- matrix(NA,nrow = rws,ncol = num_bnds)
  colnames(bnds) <- rownames(ebl[[1]]$ev_summry$ev_bnds)
  for (i in 1:rws){
    ev0[i] <- ebl[[i]]$ev_summry$ev0[1,cntrst]
    mn_ev[i] <- mean(ebl[[i]]$ev_summry$mn_ev[ ,cntrst])
    md_ev[i] <- median(ebl[[i]]$ev_summry$mn_ev[ ,cntrst])
    bnds[i, ] <- ebl[[i]]$ev_summry$ev_bnds[,cntrst]
  }
  dtfrm <- cbind.data.frame(trials,ev0,mn_ev,md_ev,bnds)
  out <- list(ebl=ebl,cntrst=cntrst,dtfrm=dtfrm,ss=ss)
  return(out)
} # end of fltn_ev_bnd_lst


fnd_estmtr_qual <- function(eblnm,Eev,ER,medev)
{
  NA20=function(x){ifelse(is.na(x),0,x)} #replaces NA with 0
  ebl <- get(eblnm)
  rws <- length(ebl)
  cntrst_nms <- intersect(colnames(ebl[[1]]$ev_summry$ev0),names(Eev))
  # Quality stats only for calculated for model comparitons with estimated expected values
  num_cntrsts <- length(cntrst_nms)
  ev_rmse.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(ev_rmse.mat) <- c("ev0","mn_ev","med_ev")
  rownames(ev_rmse.mat) <- cntrst_nms
  ev_rmsmde.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(ev_rmsmde.mat) <- c("ev0","mn_ev","med_ev")
  rownames(ev_rmsmde.mat) <- cntrst_nms
  ev_made.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(ev_made.mat) <- c("ev0","mn_ev","med_ev")
  rownames(ev_made.mat) <- cntrst_nms
  ev_madmde.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(ev_madmde.mat) <- c("ev0","mn_ev","med_ev")
  rownames(ev_madmde.mat) <- cntrst_nms
  
  RL_rmse.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(RL_rmse.mat) <- c("RL1","mn_RL","med_RL")
  rownames(RL_rmse.mat) <- cntrst_nms
  RL_made.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(RL_made.mat) <- c("RL1","mn_RL","med_RL")
  rownames(RL_made.mat) <- cntrst_nms
  
  ev_bias.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(ev_bias.mat) <- c("ev0","mn_ev","med_ev")
  rownames(ev_bias.mat) <- cntrst_nms
  RL_bias.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(RL_bias.mat) <- c("RL1","mn_RL","med_RL")
  rownames(RL_bias.mat) <- cntrst_nms
  
  ev_medbias.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(ev_medbias.mat) <- c("ev0","mn_ev","med_ev")
  rownames(ev_medbias.mat) <- cntrst_nms
  RL_medbias.mat <- matrix(NA,nrow=num_cntrsts,ncol=3)
  colnames(RL_medbias.mat) <- c("RL1","mn_RL","med_RL")
  rownames(RL_medbias.mat) <- cntrst_nms
  
  
  for (cntrst in cntrst_nms){
    sqrd_err <- matrix(NA,ncol = 3,nrow = rws)
    colnames(sqrd_err) <- c("D2_ev0","D2_mnev","D2_medev")
    RLsq_err <- matrix(NA,ncol = 3,nrow = rws)
    colnames(RLsq_err) <- c("D2_RL1","D2_mnRL","D2_medRL")
    
    ev_err <- matrix(NA,ncol = 3,nrow = rws)
    colnames(ev_err) <- c("D_ev0","D_mnev","D_medev")
    RL_err <- matrix(NA,ncol = 3,nrow = rws)
    colnames(RL_err) <- c("D_RL1","D_mnRL","D_medRL")
    
    
    sqrd_mde <- matrix(NA,ncol = 3,nrow = rws)
    colnames(sqrd_mde) <- c("D2_ev0","D2_mnev","D2_medev")
    mad_err <- matrix(NA,ncol = 3,nrow = rws)
    colnames(mad_err) <- c("ad_ev0","ad_mnev","ad_medev")
    RL_mad_err <- matrix(NA,ncol = 3,nrow = rws)
    colnames(RL_mad_err) <- c("ad_RL1","ad_mnRL","ad_medRL")
    mad_mde <- matrix(NA,ncol = 3,nrow = rws)
    colnames(mad_mde) <- c("ad_ev0","ad_mnev","ad_medev")
    for (trl in 1:rws) {
      sqrd_err[trl,"D2_ev0"] <- (ebl[[trl]]$ev_summry$ev0[1,cntrst]-Eev[cntrst])^2
      sqrd_err[trl,"D2_mnev"] <- (mean(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-Eev[cntrst])^2
      sqrd_err[trl,"D2_medev"] <- (median(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-Eev[cntrst])^2
      
      RLsq_err[trl,"D2_RL1"] <- (ebl[[trl]]$ev_summry$RL[1,cntrst]-ER[cntrst])^2
      RLsq_err[trl,"D2_mnRL"] <- (mean(ebl[[trl]]$ev_summry$mn_RL[cntrst])-ER[cntrst])^2
      RLsq_err[trl,"D2_medRL"] <- (median(ebl[[trl]]$ev_summry$md_RL[cntrst])-ER[cntrst])^2

      ev_err[trl,"D_ev0"] <- (ebl[[trl]]$ev_summry$ev0[1,cntrst]-Eev[cntrst])
      ev_err[trl,"D_mnev"] <- (mean(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-Eev[cntrst])
      ev_err[trl,"D_medev"] <- (median(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-Eev[cntrst])
      
      RL_err[trl,"D_RL1"] <- (ebl[[trl]]$ev_summry$RL[1,cntrst]-ER[cntrst])
      RL_err[trl,"D_mnRL"] <- (mean(ebl[[trl]]$ev_summry$mn_RL[cntrst])-ER[cntrst])
      RL_err[trl,"D_medRL"] <- (median(ebl[[trl]]$ev_summry$md_RL[cntrst])-ER[cntrst])

      sqrd_mde[trl,"D2_ev0"] <- (ebl[[trl]]$ev_summry$ev0[1,cntrst]-medev[cntrst])^2
      sqrd_mde[trl,"D2_mnev"] <- (mean(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-medev[cntrst])^2
      sqrd_mde[trl,"D2_medev"] <- (median(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-medev[cntrst])^2
      
      mad_err[trl,"ad_ev0"] <- abs(ebl[[trl]]$ev_summry$ev0[1,cntrst]-Eev[cntrst])
      mad_err[trl,"ad_mnev"] <- abs(mean(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-Eev[cntrst])
      mad_err[trl,"ad_medev"] <- abs(median(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-Eev[cntrst])
      
      RL_mad_err[trl,"ad_RL1"] <- abs(ebl[[trl]]$ev_summry$RL[1,cntrst]-ER[cntrst])
      RL_mad_err[trl,"ad_mnRL"] <-abs(mean(ebl[[trl]]$ev_summry$mn_RL[cntrst])-ER[cntrst])
      RL_mad_err[trl,"ad_medRL"] <- abs(median(ebl[[trl]]$ev_summry$md_RL[cntrst])-ER[cntrst])

      mad_mde[trl,"ad_ev0"] <- abs(ebl[[trl]]$ev_summry$ev0[1,cntrst]-medev[cntrst])
      mad_mde[trl,"ad_mnev"] <- abs(mean(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-medev[cntrst])
      mad_mde[trl,"ad_medev"] <- abs(median(ebl[[trl]]$ev_summry$mn_ev[,cntrst])-medev[cntrst])
    }
    
    
    ev_bias.mat[cntrst,"ev0"] <- (mean(ev_err[,"D_ev0"]))
    ev_bias.mat[cntrst,"mn_ev"] <- (mean(ev_err[,"D_mnev"]))
    ev_bias.mat[cntrst,"med_ev"] <- (mean(ev_err[,"D_medev"]))
    
    RL_bias.mat[cntrst,"RL1"] <- (mean(RL_err[,"D_RL1"]))
    RL_bias.mat[cntrst,"mn_RL"] <- (mean(RL_err[,"D_mnRL"]))
    RL_bias.mat[cntrst,"med_RL"] <- (mean(RL_err[,"D_medRL"]))

    ev_medbias.mat[cntrst,"ev0"] <- (median(ev_err[,"D_ev0"]))
    ev_medbias.mat[cntrst,"mn_ev"] <- (median(ev_err[,"D_mnev"]))
    ev_medbias.mat[cntrst,"med_ev"] <- (median(ev_err[,"D_medev"]))
    
    RL_medbias.mat[cntrst,"RL1"] <- (median(RL_err[,"D_RL1"]))
    RL_medbias.mat[cntrst,"mn_RL"] <- (median(RL_err[,"D_mnRL"]))
    RL_medbias.mat[cntrst,"med_RL"] <- (median(RL_err[,"D_medRL"]))

    ev_rmse.mat[cntrst,"ev0"] <- sqrt(mean(sqrd_err[,"D2_ev0"]))
    ev_rmse.mat[cntrst,"mn_ev"] <- sqrt(mean(sqrd_err[,"D2_mnev"]))
    ev_rmse.mat[cntrst,"med_ev"] <- sqrt(mean(sqrd_err[,"D2_medev"]))
    
    RL_rmse.mat[cntrst,"RL1"] <- sqrt(mean(RLsq_err[,"D2_RL1"]))
    RL_rmse.mat[cntrst,"mn_RL"] <- sqrt(mean(RLsq_err[,"D2_mnRL"]))
    RL_rmse.mat[cntrst,"med_RL"] <- sqrt(mean(RLsq_err[,"D2_medRL"]))

    
    ev_rmsmde.mat[cntrst,"ev0"] <- sqrt(mean(sqrd_mde[,"D2_ev0"]))
    ev_rmsmde.mat[cntrst,"mn_ev"] <- sqrt(mean(sqrd_mde[,"D2_mnev"]))
    ev_rmsmde.mat[cntrst,"med_ev"] <- sqrt(mean(sqrd_mde[,"D2_medev"]))
    
    ev_made.mat[cntrst,"ev0"] <- 1.4826*(mean(mad_err[,"ad_ev0"]))
    ev_made.mat[cntrst,"mn_ev"] <- 1.4826*(mean(mad_err[,"ad_mnev"]))
    ev_made.mat[cntrst,"med_ev"] <- 1.4826*(mean(mad_err[,"ad_medev"]))
    
    RL_made.mat[cntrst,"RL1"] <- 1.4826*(mean(RL_mad_err[,"ad_RL1"]))
    RL_made.mat[cntrst,"mn_RL"] <- 1.4826*(mean(RL_mad_err[,"ad_mnRL"]))
    RL_made.mat[cntrst,"med_RL"] <- 1.4826*(mean(RL_mad_err[,"ad_medRL"]))

    ev_madmde.mat[cntrst,"ev0"] <- 1.4826*(mean(mad_mde[,"ad_ev0"]))
    ev_madmde.mat[cntrst,"mn_ev"] <- 1.4826*(mean(mad_mde[,"ad_mnev"]))
    ev_madmde.mat[cntrst,"med_ev"] <- 1.4826*(mean(mad_mde[,"ad_medev"]))
    
  }
  ev_rmse <- aaply(.data = ev_rmse.mat,.margins = 2,.fun = mean)
  ev_rmsmde <- aaply(.data = ev_rmsmde.mat,.margins = 2,.fun = mean)
  ev_made <- aaply(.data = ev_made.mat,.margins = 2,.fun = mean)
  ev_madmde <- aaply(.data = ev_madmde.mat,.margins = 2,.fun = mean)
  
  RL_rmse <- aaply(.data = RL_rmse.mat,.margins = 2,.fun = mean)
  RL_made <- aaply(.data = RL_made.mat,.margins = 2,.fun = mean)
  
  ev_bias <- aaply(.data = ev_bias.mat,.margins = 2,.fun = mean)
  RL_bias <- aaply(.data = RL_bias.mat,.margins = 2,.fun = mean)
  ev_medbias <- aaply(.data = ev_medbias.mat,.margins = 2,.fun = median)
  RL_medbias <- aaply(.data = RL_medbias.mat,.margins = 2,.fun = median)
  
  
  EQ_ev <- rbind(ev_rmse,ev_rmsmde,ev_made,ev_madmde,ev_bias, ev_medbias)
  colnames(EQ_ev) <- c("ev0", "mn_ev", "med_ev")
  EQ_RL <- rbind(RL_rmse, RL_made, RL_bias, RL_medbias)
  colnames(EQ_RL) <- c("RL1", "mn_RL", "med_RL")
  
  out <- list(EQ_ev=EQ_ev, EQ_RL=EQ_RL, ev_rmse.mat=ev_rmse.mat, 
              ev_rmsmde.mat=ev_rmsmde.mat, ev_made.mat= ev_made.mat, 
              ev_madmde.mat=ev_madmde.mat, RL_rmse.mat=RL_rmse.mat, 
              RL_made.mat=RL_made.mat, ev_bias.mat=ev_bias.mat, 
              ev_medbias.mat=ev_medbias.mat, RL_bias.mat=RL_bias.mat,
              RL_medbias.mat=RL_medbias.mat)
  return(out)
}

plot_bnds <- function(flt_bnds,clbrtd=T,Eev, ylm.0=T, mx_ylm=200, mrgnl_lvl=4, strng_lvl=7)
{
  #mx_ylm keeps plot real in rare pathological cases
  cntrst <- flt_bnds$cntrst
  ss_txt <- str_c("Sample size is ",flt_bnds$ss,".",sep = "")
  refmd_txt <- str_c("Reference model is ", alt_nm[str_sub(string=cntrst,start=1,end=3)],".",sep="")
  altmd_txt <- str_c("Alternative model is ", alt_nm[str_sub(string=cntrst,start=5,end=7)],".",sep="")
  if (clbrtd) {
    ll <- 7
    ul <- 8
  } else {
    ll <- 5
    ul <- 6
  }
  
  lwbnd.nm <- names(flt_bnds$dtfrm)[ll]
  lwi <- str_locate(string = lwbnd.nm, pattern = "L")-1
  hibnd.nm <- names(flt_bnds$dtfrm)[ul]
  hii <- str_locate(string = hibnd.nm, pattern = "U")-1
  lwnum <- as.numeric(str_sub(string = lwbnd.nm,start = 1,end = lwi))*100
  lwstr <- as.character(lwnum[1])
  hinum <- as.numeric(str_sub(string = hibnd.nm,start = 1,end = hii))*100
  histr <- as.character(hinum[1])
  rng_txt <- str_c("(",lwstr,"%, ", histr,"%)",sep = "")
  
  if (clbrtd){
    bnd_txt <- " calibrated percentile bounds."
  } else {
    bnd_txt <- " percentile bounds."
  }
  limit_txt <- str_c("Limits are: ",rng_txt,bnd_txt,sep="")
  mn_ttl <- str_c(ss_txt,limit_txt,"\n",refmd_txt,altmd_txt,sep=" ")
  trial <- flt_bnds$dtfrm$trial
  ev0 <- flt_bnds$dtfrm$ev0
  ulim <- flt_bnds$dtfrm[,ul]
  llim <- flt_bnds$dtfrm[,ll]
  xlim <- c(1,dim(flt_bnds$dtfrm)[1])
  if (ylm.0) {ylm <- range(flt_bnds$dtfrm[,ll],flt_bnds$dtfrm[,ul],0)
  } else ylm <- range(flt_bnds$dtfrm[,ll],flt_bnds$dtfrm[,ul])
  ylm[1] <- max(ylm[1],-mx_ylm)
  ylm[1] <- min(ylm[1],-strng_lvl)
  ylm[2] <- min(ylm[2],mx_ylm)
  ylm[2] <- max(ylm[2],strng_lvl)
  
  ylim <- ylm
  par(mgp=c(2,1,0),cex.lab=1.5)
  plot(x = flt_bnds$dtfrm[,"trial"],y = flt_bnds$dtfrm[,"ev0"],xlim=xlim,ylim=ylim,pch=20,
       xlab="Trial", ylab=expression(Delta*SIC[ra]),main=mn_ttl)
  axis(side = 4, at=c(-strng_lvl,-mrgnl_lvl,0,mrgnl_lvl,strng_lvl),labels = T,
       cex.axis=1,las=1,hadj=0.5,tick=F)
  abline(h=Eev[cntrst])
  abline(h=0,lty=2)
  abline(h=c(-mrgnl_lvl,mrgnl_lvl),lty=3)
  abline(h=c(-strng_lvl,strng_lvl),lty=4)
  plotrix::dispersion(x = trial,y=ev0,ulim = ulim,llim = llim,intervals = F,lty=1,arrow.cap = 0.005, type = "a")
  #  dispersion(x = trial,y=ev0,ulim = ulim,llim = llim,intervals = F,lty=1,arrow.cap = 0.005)
}

plot_bnds_lst <- function(ebl.nm,clbrtd=T,Eev,toPDF=F,PDF.nm="bnd_plt.pdf",
                          ylm.0=T, mx_ylm=100, mrgnl_lvl=4)
{
  ebl <- get(ebl.nm)
  cntrst.vec <- names(ebl[[1]]$ev_summry$mn_RL)
  if (toPDF) pdf(file = PDF.nm)
  for (cntrst in cntrst.vec) {
    flt_bnds <- fltn_ev_bnd_lst(ebl = ebl.nm,cntrst = cntrst)
    plot_bnds(flt_bnds = flt_bnds,clbrtd ,Eev=Eev,ylm.0 = ylm.0,
              mx_ylm = mx_ylm, mrgnl_lvl = mrgnl_lvl)
  }
  if (toPDF) dev.off()
}


FltnBt <- function(dblbtobj) # Flattens a double boot strap into a single one
{
  out <- dblbtobj[[1]]
  for (i in 2:length(dblbtobj)){
    out$NPec <- rbind(out$NPec,dblbtobj[[i]]$NPec) 
    out$Pec <- rbind(out$Pec,dblbtobj[[i]]$Pec) 
  }
  return(out)
} # end of FltnBt

Ev4md <- function(ecmat,refnm)# makes tabel of evidence from matrix of criteria values
{#refnm is the name of the model presumed best or a vector of names to be used as reference models. 
  #2019-01-08 08:29:48 MST bug fix Ev4md now correctly handles single row matricies
  #2019-01-19 11:08:11 MST bug fix Ev4md now correctly handles single column matricies
  out <- NULL
  if (is.vector(ecmat)) {ecmat <- t(ecmat)} # turns  a row ecmat back into a matrix 
  nmods <- dim(ecmat)[2]
  for (i in 1:length(refnm)) {
    evmat <- (ecmat - ecmat[,refnm[i]])
    nms <- setdiff(colnames(ecmat),refnm[i])
    evmat <- evmat[,nms,drop=FALSE]
    cnms <- str_c(refnm[i],setdiff(colnames(ecmat),refnm[i]),sep=".")
    colnames(evmat) <- cnms
    out <- cbind(out,evmat)
  }
  attr(out,"evType") <- attr(ecmat,"evType")
  return(out)
}# makes table of evidence from matrix of criteria values

Ev4 <- function(ecmat,cmp_lst) # Ev4md wrapper--allows finer comparison control 
{
  {# cmp_lst: is a list of lists. Each element has 2 vectors of 1 or more model names.
    #        :alt_nms and ref_nms.    
  }# end of header documentation
  {# main body
    out <- NULL
    lst_len <- length(cmp_lst)
    for (i in 1:lst_len){
      all_nms.i <- unlist(cmp_lst[[i]])
      ecmat.i <- ecmat[,all_nms.i]
      out <- cbind(out,Ev4md(ecmat=ecmat.i,refnm = cmp_lst[[i]]$ref_nms))
    }
    return(out)
  }# end of main body
  
}# end Ev4


aM2RL <- function(ev_smmry) # transforms old data set from aM to RL 
{
   out <- ev_smmry 
   if(is.element(el = "mn_aM",set = names(ev_smmry))){
     mn_RL <- 1-ev_smmry$mn_aM
     md_RL <- 1-ev_smmry$md_aM
     RL_LB <- 1-ev_smmry$aM_UB
     rownames(RL_LB) <- str_c(100-as.numeric(str_sub(rownames(RL_LB),1,str_locate(rownames(RL_LB),"%")[,1]-1)),"%_LB",sep="")
     RL <- 1 - ev_smmry$aM
   }
   out=list.remove(out,c("mn_aM","md_aM","aM_UB","aM","bca_aM_UB"))
   out$mn_RL <- mn_RL
   out$md_RL <- md_RL
   out$RL_LB <- RL_LB
   out$RL <- RL
   return(out)
} # end of aM2RL

aM2RL.lst <- function(ev_c.lst)
{
  rplc <- function(X){
    out <- X
    out$ev_summry <- aM2RL(X$ev_summry)
    return(out)
  }
  out <- rlist::list.apply(ev_c.lst,rplc)
  return(out)
}

is.in_rng <- function(x,rng) {
  out <- ((x >= min(rng)) & (x <= max(rng)))
return(out)
}
  
fnd_cvr <- function(ev_bnd.lst,Eev.vec, ER.vec, medev.vec,
                    lbpc="0.05Lprcnt", ubpc ="0.95Uprcnt",
                    lbcl="0.05Lclbrtd", ubcl ="0.95Uclbrtd",
                    mrgnl_lvl=4, closure_type="cc") # coverage stats
{
  {# initializations
    num_trls <- length(ev_bnd.lst)
    num_cntrsts <- dim(ev_bnd.lst[[1]]$ev_summry$ev0)[2]
    cntrst_nms <- colnames(ev_bnd.lst[[1]]$ev_summry$ev0)
    num_RL_LB <- dim(ev_bnd.lst[[1]]$ev_summry$RL_LB)[1]
    RL_LB_nms <- rownames(ev_bnd.lst[[1]]$ev_summry$RL_LB)
    cvr_pcntl.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    cvr_clbrt.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    colnames(cvr_pcntl.mat) <- cntrst_nms
    colnames(cvr_clbrt.mat) <- cntrst_nms
    mislo_pcntl.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    mislo_clbrt.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    colnames(mislo_pcntl.mat) <- cntrst_nms
    colnames(mislo_clbrt.mat) <- cntrst_nms
    mishi_pcntl.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    mishi_clbrt.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    colnames(mishi_pcntl.mat) <- cntrst_nms
    colnames(mishi_clbrt.mat) <- cntrst_nms
    intlen_pcntl.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    intlen_clbrt.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    colnames(intlen_pcntl.mat) <- cntrst_nms
    colnames(intlen_clbrt.mat) <- cntrst_nms
    secure_pcntl.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    misldng_pcntl.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    secure_clbrt.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    misldng_clbrt.mat <- matrix(NA,nrow=num_trls,ncol=num_cntrsts)
    colnames(secure_pcntl.mat) <- cntrst_nms
    colnames(misldng_pcntl.mat) <- cntrst_nms
    colnames(secure_clbrt.mat) <- cntrst_nms
    colnames(misldng_clbrt.mat) <- cntrst_nms
    RL_LB_cvr.mat <- matrix(NA,nrow=(num_RL_LB + 1),ncol=num_cntrsts)
    colnames(RL_LB_cvr.mat) <- cntrst_nms
    rownames(RL_LB_cvr.mat) <- c("ER", RL_LB_nms)
    RL_LB_cvr.mat["ER",] <- ER.vec[cntrst_nms] # first row of cover mat has true values for contrasts considered
  }# end initializations
  { # p_cvr calculations
    for (trl in 1:num_trls) {
      for (cntrst in cntrst_nms){
        Eev_cntrst <- Eev.vec[cntrst]
        
        prcntrng <- c(ev_bnd.lst[[trl]]$ev_summry$ev_bnds[lbpc,cntrst],
                      ev_bnd.lst[[trl]]$ev_summry$ev_bnds[ubpc,cntrst])
        cvr_pcntl.mat[trl,cntrst] <- is.in_rng(x=Eev_cntrst,rng=prcntrng)
        
        clbrtdrng <- c(ev_bnd.lst[[trl]]$ev_summry$ev_bnds[lbcl,cntrst],
                       ev_bnd.lst[[trl]]$ev_summry$ev_bnds[ubcl,cntrst])
        cvr_clbrt.mat[trl,cntrst] <- is.in_rng(x=Eev_cntrst,rng=clbrtdrng)
        
        mislo_pcntl.mat[trl,cntrst] <- (Eev_cntrst > ev_bnd.lst[[trl]]$ev_summry$ev_bnds[ubpc,cntrst])
        mishi_pcntl.mat[trl,cntrst] <- (Eev_cntrst < ev_bnd.lst[[trl]]$ev_summry$ev_bnds[lbpc,cntrst])
        mislo_clbrt.mat[trl,cntrst] <- (Eev_cntrst > ev_bnd.lst[[trl]]$ev_summry$ev_bnds[ubcl,cntrst])
        mishi_clbrt.mat[trl,cntrst] <- (Eev_cntrst < ev_bnd.lst[[trl]]$ev_summry$ev_bnds[lbcl,cntrst])
        
        intlen_pcntl.mat[trl,cntrst] <- max(prcntrng)-min(prcntrng)
        intlen_clbrt.mat[trl,cntrst] <- max(clbrtdrng)-min(clbrtdrng)
        
        secure_pcntl.mat[trl,cntrst] <- (ev_bnd.lst[[trl]]$ev_summry$ev_bnds[lbpc,cntrst] > mrgnl_lvl)  
        secure_clbrt.mat[trl,cntrst] <- (ev_bnd.lst[[trl]]$ev_summry$ev_bnds[lbcl,cntrst] > mrgnl_lvl)  
        misldng_pcntl.mat[trl,cntrst] <- (ev_bnd.lst[[trl]]$ev_summry$ev_bnds[ubpc,cntrst] < -mrgnl_lvl)  
        misldng_clbrt.mat[trl,cntrst] <- (ev_bnd.lst[[trl]]$ev_summry$ev_bnds[ubcl,cntrst] < -mrgnl_lvl) 
        
      }
    }
    p_cvr_pcntl <- aaply(.data = cvr_pcntl.mat,.margins = 2,.fun = mean, na.rm=T)
    p_cvr_clbrt <- aaply(.data = cvr_clbrt.mat,.margins = 2,.fun = mean, na.rm=T)
    p_mislo_pcntl <- aaply(.data=mislo_pcntl.mat,.margins = 2,.fun = mean, na.rm=T)
    p_mishi_pcntl <- aaply(.data=mishi_pcntl.mat,.margins = 2,.fun = mean, na.rm=T)
    p_mislo_clbrt <- aaply(.data=mislo_clbrt.mat,.margins = 2,.fun = mean, na.rm=T)
    p_mishi_clbrt <- aaply(.data=mishi_clbrt.mat,.margins = 2,.fun = mean, na.rm=T)
    p_secure_pcntl <- aaply(.data=secure_pcntl.mat,.margins = 2,.fun = mean, na.rm=T)
    p_secure_clbrt <- aaply(.data=secure_clbrt.mat,.margins = 2,.fun = mean, na.rm=T)
    p_misldng_pcntl <- aaply(.data=misldng_pcntl.mat,.margins = 2,.fun = mean, na.rm=T)
    p_misldng_clbrt <- aaply(.data=misldng_clbrt.mat,.margins = 2,.fun = mean, na.rm=T)
    mn_intlen_pcntl <- aaply(.data=intlen_pcntl.mat,.margins = 2,.fun = mean, na.rm=T)
    mn_intlen_clbrt <- aaply(.data=intlen_clbrt.mat,.margins = 2,.fun = mean, na.rm=T)
    
    names(p_cvr_pcntl) <- names(p_cvr_clbrt) <- cntrst_nms
    names(p_mislo_pcntl) <- names(p_mislo_clbrt) <- cntrst_nms
    names(p_mishi_pcntl) <- names(p_mishi_clbrt) <- cntrst_nms
    names(p_secure_pcntl) <- names(p_secure_clbrt) <- cntrst_nms
    names(p_misldng_pcntl) <- names(p_misldng_clbrt) <- cntrst_nms
    names(mn_intlen_pcntl) <- names(mn_intlen_clbrt) <- cntrst_nms
    
    cvr_stats <- rbind(p_cvr_pcntl,p_cvr_clbrt,p_mislo_pcntl,p_mislo_clbrt,
                       p_mishi_pcntl,p_mishi_clbrt,p_secure_pcntl,p_secure_clbrt,
                       p_misldng_pcntl,p_misldng_clbrt)
    
    mn_intlen <- rbind(mn_intlen_pcntl,mn_intlen_clbrt)
  } # end p_cvr calculations
  { # RL calculations
    for (rllb in RL_LB_nms) {
      for (cntrst in cntrst_nms){
        ER_cntrst <- ER.vec[cntrst]
        RL_cvr_cntrstXtrl <- rep(NA, times = num_trls)
        for (trl in 1:num_trls){
          rng <- c(ev_bnd.lst[[trl]]$ev_summry$RL_LB[rllb,cntrst],1)
          RL_cvr_cntrstXtrl[trl] <- is.in_rng(x=ER_cntrst,rng=rng)
        }
        RL_LB_cvr.mat[rllb,cntrst] <- mean(RL_cvr_cntrstXtrl)
        
      } # end for cntrst
      
    } # end for rllb
  } # end RL calculations
  out <- list(cvr_stats=cvr_stats, mn_intlen=mn_intlen, RL_LB_cvr.mat=RL_LB_cvr.mat)
  return(out)
} # end fnd_cvr

conv_cntrst.nm <- function(cntrst.nm){
  cn <- cntrst.nm
  out <- str_c(alt_nm[str_sub(string = cn,start = 1,end = 3)]," ",
               alt_nm[str_sub(string = cn,start = 5,end = 8)]) 
  return(out)
}


#####

m13nms=c( "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m25")
crct_nms <- c("m13.m10", "m13.m11", "m13.m15", "m13.m17", "m13.m18", "m13.m19", "m13.m20") 
mis_nms <- c("m10.m17", "m18.m20", "m15.m25", "m10.m15", "m19.m20", "m17.m25")

pcntl_cvr.nms <- c("p_cvr_pcntl", "p_mislo_pcntl", "p_mishi_pcntl",
                   "p_secure_pcntl" , "p_misldng_pcntl" )
clbrt_cvr.nms <- c("p_cvr_clbrt", "p_mislo_clbrt", "p_mishi_clbrt",
                   "p_secure_clbrt" , "p_misldng_clbrt" )


alt_nm <- rep(NA,12)
names(alt_nm) <- m13nms
alt_nm["m10"] <- "UVWXYQ" 
alt_nm["m11"] <- "_VWXYQ" 
alt_nm["m12"] <- "U_WXYQ" 
alt_nm["m13"] <- "__WXYQ"
alt_nm["m14"] <- "UVW_YQ"
alt_nm["m15"] <- "UV_XYQ" 
alt_nm["m16"] <- "UVWX__"  
alt_nm["m17"] <- "UVWXY_" 
alt_nm["m18"] <- "___XYQ" 
alt_nm["m19"] <- "__WXY_" 
alt_nm["m20"] <- "___XY_" 
alt_nm["m25"] <- "UV_XY_" 



cmp_lst <- list()
cmp_lst[[1]] <- list(alt_nms=m13nms,ref_nms="m13")
cmp_lst[[2]] <- list(alt_nms="m17",ref_nms="m10")
cmp_lst[[3]] <- list(alt_nms="m20",ref_nms="m18")
cmp_lst[[4]] <- list(alt_nms="m25",ref_nms="m15")
cmp_lst[[5]] <- list(alt_nms="m15",ref_nms="m10")
cmp_lst[[6]] <- list(alt_nms="m20",ref_nms="m19")
cmp_lst[[7]] <- list(alt_nms="m25",ref_nms="m17")


##huge data begin

ev_cl_hugeD.lst <- readRDS("ev_cl_hugeD.lst.rds")
ev_cl_hugeD.lst <-aM2RL.lst(ev_cl_hugeD.lst)
nms <- colnames(ev_cl_hugeD.lst[[1]]$ev_summry$ev_bnds)

hugeDm10.dblbt <- readRDS(file = "hugeDm10.dblbt.rds")
hugeDm10.fltbt <- FltnBt(hugeDm10.dblbt)
hugeDm10.flt.Pevmat <- Ev4(hugeDm10.fltbt$Pec,cmp_lst = cmp_lst)
rownames(hugeDm10.flt.Pevmat) <- NULL
hugeDm10_Eev <- aaply(.data = hugeDm10.flt.Pevmat,.margins = 2,.fun = mean)
hugeDm10_medev <- aaply(.data = hugeDm10.flt.Pevmat,.margins = 2,.fun = median)
hugeDm10_ER <- aaply(.data = hugeDm10.flt.Pevmat,.margins = 2,
                      .fun = function(X){mean(X > 0)})

hugeD_ev.EQ <- fnd_estmtr_qual(eblnm = "ev_cl_hugeD.lst",
                               Eev = hugeDm10_Eev[nms],ER = hugeDm10_ER[nms], 
                               medev =hugeDm10_medev[nms])


plot_bnds_lst(ebl.nm = "ev_cl_hugeD.lst",toPDF = T,
              PDF.nm = "ev_cl_hugeD_bnds_plt.pdf" ,clbrtd = T,
              Eev = hugeDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)

plot_bnds_lst(ebl.nm = "ev_cl_hugeD.lst",toPDF = T,
              PDF.nm = "ev_pl_hugeD_bnds_plt.pdf" ,clbrtd = F,
              Eev = hugeDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)



hugeD_cl.cvr <- fnd_cvr(ev_bnd.lst = ev_cl_hugeD.lst,Eev.vec = hugeDm10_Eev,
                       ER.vec = hugeDm10_ER, medev.vec = hugeDm10_medev )

hugeD.crct.pcntl.tbl=hugeD_cl.cvr$cvr_stats[pcntl_cvr.nms,crct_nms]
colnames(hugeD.crct.pcntl.tbl) <-conv_cntrst.nm(colnames(hugeD.crct.pcntl.tbl[,crct_nms]))
hugeD.crct.pcntl.tbl=format(hugeD.crct.pcntl.tbl,digits = 3)
write.csv(hugeD.crct.pcntl.tbl,file = "hugeD.crct.pcntl.tbl.txt",quote = F)

hugeD.mis.pcntl.tbl=hugeD_cl.cvr$cvr_stats[pcntl_cvr.nms,mis_nms]
colnames(hugeD.mis.pcntl.tbl) <-conv_cntrst.nm(colnames(hugeD.mis.pcntl.tbl))
hugeD.mis.pcntl.tbl=format(hugeD.mis.pcntl.tbl,digits = 3)
write.csv(hugeD.mis.pcntl.tbl,file = "hugeD.mis.pcntl.tbl.txt",quote = F)

hugeD.crct.clbrt.tbl=hugeD_cl.cvr$cvr_stats[clbrt_cvr.nms,crct_nms]
colnames(hugeD.crct.clbrt.tbl) <-conv_cntrst.nm(colnames(hugeD.crct.clbrt.tbl[,crct_nms]))
hugeD.crct.clbrt.tbl=format(hugeD.crct.clbrt.tbl,digits = 3)
write.csv(hugeD.crct.clbrt.tbl,file = "hugeD.crct.clbrt.tbl.txt",quote = F)

hugeD.mis.clbrt.tbl=hugeD_cl.cvr$cvr_stats[clbrt_cvr.nms,mis_nms]
colnames(hugeD.mis.clbrt.tbl) <-conv_cntrst.nm(colnames(hugeD.mis.clbrt.tbl))
hugeD.mis.clbrt.tbl=format(hugeD.mis.clbrt.tbl,digits = 3)
write.csv(hugeD.mis.clbrt.tbl,file = "hugeD.mis.clbrt.tbl.txt",quote = F)

hugeD_ev_rmse.tbl <- hugeD_ev.EQ$ev_rmse.mat
rownames(hugeD_ev_rmse.tbl) <-conv_cntrst.nm(rownames(hugeD_ev_rmse.tbl))
hugeD_ev_rmse.tbl=format(hugeD_ev_rmse.tbl,digits = 3)
write.csv(hugeD_ev_rmse.tbl,file = "hugeD_ev_rmse.tbl.txt",quote = F)

hugeD_RL_rmse.tbl <- hugeD_ev.EQ$RL_rmse.mat
rownames(hugeD_RL_rmse.tbl) <-rownames(hugeD_ev_rmse.tbl)
hugeD_RL_rmse.tbl=format(hugeD_RL_rmse.tbl,digits = 3)
write.csv(hugeD_RL_rmse.tbl,file = "hugeD_RL_rmse.tbl.txt",quote = F)


hugeD.EQ_RL.tbl <- format(hugeD_ev.EQ$EQ_RL,digits=3)
write.csv(hugeD.EQ_RL.tbl,file = "hugeD.EQ_RL.tbl.txt",quote = F)

hugeD.EQ_ev.tbl <- format(hugeD_ev.EQ$EQ_ev, digits=3)
write.csv(hugeD.EQ_ev.tbl,file = "hugeD.EQ_ev.tbl.txt",quote = F)

###huge data end


##big data begin
# some adjustment in this group because not all models calculated in all files
bignms=m13nms[1:11]
big.cmp_lst <- list()
big.cmp_lst[[1]] <- list(alt_nms=bignms,ref_nms="m13")
big.cmp_lst[[2]] <- list(alt_nms="m17",ref_nms="m10")
big.cmp_lst[[3]] <- list(alt_nms="m20",ref_nms="m18")
big.cmp_lst[[4]] <- list(alt_nms="m15",ref_nms="m10")
big.cmp_lst[[5]] <- list(alt_nms="m20",ref_nms="m19")


ev_cl_bigD.lst <- readRDS("ev_cl_bigD.lst.rds")
ev_cl_bigD.lst <-aM2RL.lst(ev_cl_bigD.lst)
#nms <- colnames(ev_cl_bigD.lst[[1]]$ev_summry$ev_bnds)
nms <- names(bigDm10_Eev)
bigDm10.dblbt <- readRDS(file = "bigDm10.dblbt.rds")
bigDm10.fltbt <- FltnBt(bigDm10.dblbt)
bigDm10.flt.Pevmat <- Ev4(bigDm10.fltbt$Pec,cmp_lst = big.cmp_lst)
rownames(bigDm10.flt.Pevmat) <- NULL
bigDm10_Eev <- aaply(.data = bigDm10.flt.Pevmat,.margins = 2,.fun = mean)
bigDm10_medev <- aaply(.data = bigDm10.flt.Pevmat,.margins = 2,.fun = median)
bigDm10_ER <- aaply(.data = bigDm10.flt.Pevmat,.margins = 2,
                     .fun = function(X){mean(X > 0)})

bigD_ev.EQ <- fnd_estmtr_qual(eblnm = "ev_cl_bigD.lst",
                               Eev = bigDm10_Eev[nms],ER = bigDm10_ER[nms], 
                               medev =bigDm10_medev[nms])


plot_bnds_lst(ebl.nm = "ev_cl_bigD.lst",toPDF = T,
              PDF.nm = "ev_cl_bigD_bnds_plt.pdf" ,clbrtd = T,
              Eev = bigDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)

plot_bnds_lst(ebl.nm = "ev_cl_bigD.lst",toPDF = T,
              PDF.nm = "ev_pl_bigD_bnds_plt.pdf" ,clbrtd = F,
              Eev = bigDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)



bigD_cl.cvr <- fnd_cvr(ev_bnd.lst = ev_cl_bigD.lst,Eev.vec = bigDm10_Eev,
                        ER.vec = bigDm10_ER, medev.vec = bigDm10_medev )

bigD.crct.pcntl.tbl=bigD_cl.cvr$cvr_stats[pcntl_cvr.nms,crct_nms]
colnames(bigD.crct.pcntl.tbl) <-conv_cntrst.nm(colnames(bigD.crct.pcntl.tbl[,crct_nms]))

bigD.crct.pcntl.tbl=format(bigD.crct.pcntl.tbl,digits = 3)
write.csv(bigD.crct.pcntl.tbl,file = "bigD.crct.pcntl.tbl.txt",quote = F)

bigD.mis.pcntl.tbl=bigD_cl.cvr$cvr_stats[pcntl_cvr.nms,mis_nms]
colnames(bigD.mis.pcntl.tbl) <-conv_cntrst.nm(colnames(bigD.mis.pcntl.tbl))
bigD.mis.pcntl.tbl=format(bigD.mis.pcntl.tbl,digits = 3)
write.csv(bigD.mis.pcntl.tbl,file = "bigD.mis.pcntl.tbl.txt",quote = F)

bigD.crct.clbrt.tbl=bigD_cl.cvr$cvr_stats[clbrt_cvr.nms,crct_nms]
colnames(bigD.crct.clbrt.tbl) <-conv_cntrst.nm(colnames(bigD.crct.clbrt.tbl[,crct_nms]))
bigD.crct.clbrt.tbl=format(bigD.crct.clbrt.tbl,digits = 3)
write.csv(bigD.crct.clbrt.tbl,file = "bigD.crct.clbrt.tbl.txt",quote = F)

bigD.mis.clbrt.tbl=bigD_cl.cvr$cvr_stats[clbrt_cvr.nms,mis_nms]
colnames(bigD.mis.clbrt.tbl) <-conv_cntrst.nm(colnames(bigD.mis.clbrt.tbl))
bigD.mis.clbrt.tbl=format(bigD.mis.clbrt.tbl,digits = 3)
write.csv(bigD.mis.clbrt.tbl,file = "bigD.mis.clbrt.tbl.txt",quote = F)

bigD_RL_rmse.tbl <- bigD_ev.EQ$RL_rmse.mat
rownames(bigD_RL_rmse.tbl) <-rownames(bigD_ev_rmse.tbl)
bigD_RL_rmse.tbl=format(bigD_RL_rmse.tbl,digits = 3)
write.csv(bigD_RL_rmse.tbl,file = "bigD_RL_rmse.tbl.txt",quote = F)


bigD_ev_rmse.tbl <- bigD_ev.EQ$ev_rmse.mat
rownames(bigD_ev_rmse.tbl) <-conv_cntrst.nm(rownames(bigD_ev_rmse.tbl))
bigD_ev_rmse.tbl=format(bigD_ev_rmse.tbl,digits = 3)
write.csv(bigD_ev_rmse.tbl,file = "bigD_ev_rmse.tbl.txt",quote = F)


bigD.EQ_RL.tbl <- format(bigD_ev.EQ$EQ_RL,digits=3)
write.csv(bigD.EQ_RL.tbl,file = "bigD.EQ_RL.tbl.txt",quote = F)

bigD.EQ_ev.tbl <- format(bigD_ev.EQ$EQ_ev, digits=3)
write.csv(bigD.EQ_ev.tbl,file = "bigD.EQ_ev.tbl.txt",quote = F)

###big data end

##small data begin

ev_cl_smallD.lst <- readRDS("ev_cl_smallD.lst.rds")
ev_cl_smallD.lst <-aM2RL.lst(ev_cl_smallD.lst)
nms <- colnames(ev_cl_smallD.lst[[1]]$ev_summry$ev_bnds)

smallDm10.dblbt <- readRDS(file = "smallDm10.dblbt.rds")
smallDm10.fltbt <- FltnBt(smallDm10.dblbt)
smallDm10.flt.Pevmat <- Ev4(smallDm10.fltbt$Pec,cmp_lst = cmp_lst)
rownames(smallDm10.flt.Pevmat) <- NULL
smallDm10_Eev <- aaply(.data = smallDm10.flt.Pevmat,.margins = 2,.fun = mean)
smallDm10_medev <- aaply(.data = smallDm10.flt.Pevmat,.margins = 2,.fun = median)
smallDm10_ER <- aaply(.data = smallDm10.flt.Pevmat,.margins = 2,
                     .fun = function(X){mean(X > 0)})

smallD_ev.EQ <- fnd_estmtr_qual(eblnm = "ev_cl_smallD.lst",
                               Eev = smallDm10_Eev[nms],ER = smallDm10_ER[nms], 
                               medev =smallDm10_medev[nms])


plot_bnds_lst(ebl.nm = "ev_cl_smallD.lst",toPDF = T,
              PDF.nm = "ev_cl_smallD_bnds_plt.pdf" ,clbrtd = T,
              Eev = smallDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)

plot_bnds_lst(ebl.nm = "ev_cl_smallD.lst",toPDF = T,
              PDF.nm = "ev_pl_smallD_bnds_plt.pdf" ,clbrtd = F,
              Eev = smallDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)



smallD_cl.cvr <- fnd_cvr(ev_bnd.lst = ev_cl_smallD.lst,Eev.vec = smallDm10_Eev,
                        ER.vec = smallDm10_ER, medev.vec = smallDm10_medev )

smallD.crct.pcntl.tbl=smallD_cl.cvr$cvr_stats[pcntl_cvr.nms,crct_nms]
colnames(smallD.crct.pcntl.tbl) <-conv_cntrst.nm(colnames(smallD.crct.pcntl.tbl[,crct_nms]))

smallD.crct.pcntl.tbl=format(smallD.crct.pcntl.tbl,digits = 3)
write.csv(smallD.crct.pcntl.tbl,file = "smallD.crct.pcntl.tbl.txt",quote = F)

smallD.mis.pcntl.tbl=smallD_cl.cvr$cvr_stats[pcntl_cvr.nms,mis_nms]
colnames(smallD.mis.pcntl.tbl) <-conv_cntrst.nm(colnames(smallD.mis.pcntl.tbl))
smallD.mis.pcntl.tbl=format(smallD.mis.pcntl.tbl,digits = 3)
write.csv(smallD.mis.pcntl.tbl,file = "smallD.mis.pcntl.tbl.txt",quote = F)

smallD.crct.clbrt.tbl=smallD_cl.cvr$cvr_stats[clbrt_cvr.nms,crct_nms]
colnames(smallD.crct.clbrt.tbl) <-conv_cntrst.nm(colnames(smallD.crct.clbrt.tbl[,crct_nms]))
smallD.crct.clbrt.tbl=format(smallD.crct.clbrt.tbl,digits = 3)
write.csv(smallD.crct.clbrt.tbl,file = "smallD.crct.clbrt.tbl.txt",quote = F)

smallD.mis.clbrt.tbl=smallD_cl.cvr$cvr_stats[clbrt_cvr.nms,mis_nms]
colnames(smallD.mis.clbrt.tbl) <-conv_cntrst.nm(colnames(smallD.mis.clbrt.tbl))
smallD.mis.clbrt.tbl=format(smallD.mis.clbrt.tbl,digits = 3)
write.csv(smallD.mis.clbrt.tbl,file = "smallD.mis.clbrt.tbl.txt",quote = F)

smallD_RL_rmse.tbl <- smallD_ev.EQ$RL_rmse.mat
rownames(smallD_RL_rmse.tbl) <-rownames(smallD_ev_rmse.tbl)
smallD_RL_rmse.tbl=format(smallD_RL_rmse.tbl,digits = 3)
write.csv(smallD_RL_rmse.tbl,file = "smallD_RL_rmse.tbl.txt",quote = F)

smallD_ev_rmse.tbl <- smallD_ev.EQ$ev_rmse.mat
rownames(smallD_ev_rmse.tbl) <-conv_cntrst.nm(rownames(smallD_ev_rmse.tbl))
smallD_ev_rmse.tbl=format(smallD_ev_rmse.tbl,digits = 3)
write.csv(smallD_ev_rmse.tbl,file = "smallD_ev_rmse.tbl.txt",quote = F)


smallD.EQ_RL.tbl <- format(smallD_ev.EQ$EQ_RL,digits=3)
write.csv(smallD.EQ_RL.tbl,file = "smallD.EQ_RL.tbl.txt",quote = F)

smallD.EQ_ev.tbl <- format(smallD_ev.EQ$EQ_ev, digits=3)
write.csv(smallD.EQ_ev.tbl,file = "smallD.EQ_ev.tbl.txt",quote = F)

###small data end

##tiny data begin

ev_cl_tinyD.lst <- readRDS("ev_cl_tinyD.lst.rds")
ev_cl_tinyD.lst <- aM2RL.lst(ev_cl_tinyD.lst)
nms <- colnames(ev_cl_tinyD.lst[[1]]$ev_summry$ev_bnds)

tinyDm10.dblbt <- readRDS(file = "tinyDm10.dblbt.rds")
tinyDm10.fltbt <- FltnBt(tinyDm10.dblbt)
tinyDm10.flt.Pevmat <- Ev4(tinyDm10.fltbt$Pec,cmp_lst = cmp_lst)
rownames(tinyDm10.flt.Pevmat) <- NULL
tinyDm10_Eev <- aaply(.data = tinyDm10.flt.Pevmat,.margins = 2,.fun = mean)
tinyDm10_medev <- aaply(.data = tinyDm10.flt.Pevmat,.margins = 2,.fun = median)
tinyDm10_ER <- aaply(.data = tinyDm10.flt.Pevmat,.margins = 2,
                     .fun = function(X){mean(X > 0)})

tinyD_ev.EQ <- fnd_estmtr_qual(eblnm = "ev_cl_tinyD.lst",
                               Eev = tinyDm10_Eev[nms],ER = tinyDm10_ER[nms], 
                               medev =tinyDm10_medev[nms])


plot_bnds_lst(ebl.nm = "ev_cl_tinyD.lst",toPDF = T,
              PDF.nm = "ev_cl_tinyD_bnds_plt.pdf" ,clbrtd = T,
              Eev = tinyDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)

plot_bnds_lst(ebl.nm = "ev_cl_tinyD.lst",toPDF = T,
              PDF.nm = "ev_pl_tinyD_bnds_plt.pdf" ,clbrtd = F,
              Eev = tinyDm10_Eev, ylm.0 = T,mrgnl_lvl = 4)



tinyD_cl.cvr <- fnd_cvr(ev_bnd.lst = ev_cl_tinyD.lst,Eev.vec = tinyDm10_Eev,
                        ER.vec = tinyDm10_ER, medev.vec = tinyDm10_medev )

tinyD.crct.pcntl.tbl=tinyD_cl.cvr$cvr_stats[pcntl_cvr.nms,crct_nms]
colnames(tinyD.crct.pcntl.tbl) <-conv_cntrst.nm(colnames(tinyD.crct.pcntl.tbl[,crct_nms]))

tinyD.crct.pcntl.tbl=format(tinyD.crct.pcntl.tbl,digits = 3)
write.csv(tinyD.crct.pcntl.tbl,file = "tinyD.crct.pcntl.tbl.txt",quote = F)

tinyD.mis.pcntl.tbl=tinyD_cl.cvr$cvr_stats[pcntl_cvr.nms,mis_nms]
colnames(tinyD.mis.pcntl.tbl) <-conv_cntrst.nm(colnames(tinyD.mis.pcntl.tbl))
tinyD.mis.pcntl.tbl=format(tinyD.mis.pcntl.tbl,digits = 3)
write.csv(tinyD.mis.pcntl.tbl,file = "tinyD.mis.pcntl.tbl.txt",quote = F)

tinyD.crct.clbrt.tbl=tinyD_cl.cvr$cvr_stats[clbrt_cvr.nms,crct_nms]
colnames(tinyD.crct.clbrt.tbl) <-conv_cntrst.nm(colnames(tinyD.crct.clbrt.tbl[,crct_nms]))
tinyD.crct.clbrt.tbl=format(tinyD.crct.clbrt.tbl,digits = 3)
write.csv(tinyD.crct.clbrt.tbl,file = "tinyD.crct.clbrt.tbl.txt",quote = F)

tinyD.mis.clbrt.tbl=tinyD_cl.cvr$cvr_stats[clbrt_cvr.nms,mis_nms]
colnames(tinyD.mis.clbrt.tbl) <-conv_cntrst.nm(colnames(tinyD.mis.clbrt.tbl))
tinyD.mis.clbrt.tbl=format(tinyD.mis.clbrt.tbl,digits = 3)
write.csv(tinyD.mis.clbrt.tbl,file = "tinyD.mis.clbrt.tbl.txt",quote = F)

tinyD_RL_rmse.tbl <- tinyD_ev.EQ$RL_rmse.mat
rownames(tinyD_RL_rmse.tbl) <-rownames(tinyD_ev_rmse.tbl)
tinyD_RL_rmse.tbl=format(tinyD_RL_rmse.tbl,digits = 3)
write.csv(tinyD_RL_rmse.tbl,file = "tinyD_RL_rmse.tbl.txt",quote = F)

tinyD_ev_rmse.tbl <- tinyD_ev.EQ$ev_rmse.mat
rownames(tinyD_ev_rmse.tbl) <-conv_cntrst.nm(rownames(tinyD_ev_rmse.tbl))
tinyD_ev_rmse.tbl=format(tinyD_ev_rmse.tbl,digits = 3)
write.csv(tinyD_ev_rmse.tbl,file = "tinyD_ev_rmse.tbl.txt",quote = F)


tinyD.EQ_RL.tbl <- format(tinyD_ev.EQ$EQ_RL,digits=3)
write.csv(tinyD.EQ_RL.tbl,file = "tinyD.EQ_RL.tbl.txt",quote = F)

tinyD.EQ_ev.tbl <- format(tinyD_ev.EQ$EQ_ev, digits=3)
write.csv(tinyD.EQ_ev.tbl,file = "tinyD.EQ_ev.tbl.txt",quote = F)

###tiny data end















#pdf(file = "cv_convergence.pdf")
jpeg(file = "cv_convergence.jpg")
par(mfrow=c(4,1), mar=c(2,1,1,1))
Plot.dblbt(M.lst = tinyD_m13_m19.M.lst,fltbt.M = tinyD_m13_m19.fltbt.M,evnm = "m13.m19",mnTitle = "",xlm.0=T)
title(main = "Sample Size = 20      ",line = -1, adj=1)
Plot.dblbt(M.lst = bigD_m13_m19.M.lst,fltbt.M = bigD_m13_m19.fltbt.M,evnm = "m13.m19",mnTitle = "",xlm.0=T)
title(main = "Sample Size = 200     ",line = -1, adj=1)
Plot.dblbt(M.lst = hugeD_m13_m19.M.lst,fltbt.M = hugeD_m13_m19.fltbt.M,evnm = "m13.m19",mnTitle = "",xlm.0=T)
title(main = "Sample Size = 2,000   ",line = -1, adj=1)
Plot.dblbt(M.lst = jumboD_m13_m19.M.lst,fltbt.M = jumboD_m13_m19.fltbt.M,evnm = "m13.m19",mnTitle = "",xlm.0=T)
title(main = "Sample Size = 20,000 ",line = -1, adj=1)
par(mfrow=c(1,1))
dev.off()




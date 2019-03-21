#Singe file to generate simulation results

loadlib<-function(){
  library(stringr) #rational string functions
  library(foreach)  # iteration and looping compatible with parallel processing
  library(mvtnorm)  # r multivariate normal random numbers    
  library(tictoc) # timing functions
  library(plyr)  #rational wrapper for apply like functions
  library(lava)  # Latent variable analysis
  library(rlist)
  library(utils)
  library(moments)
  library(doFuture)
}

loadlib()
##### Functions



CutBySize <- function (m, block.size, nb = ceiling(m/block.size)) {
  #Function to break a long index in blocks as evenly as possible for parallel computing
  #https://privefl.github.io/blog/a-guide-to-parallelism-in-r/ written by Florian Prive'
  #m=max index, can cut by block.size or number of blocks (nb)
  #All actual blocks sizes will be within 1 of each other
  if (nb > m)
    nb <- m
  int <- m/nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  cbind(lower, upper, size)
}#break a long index in blocks


nobs.lvmfit <- function(lvm.ob) #returns the number of observation in lvm.fit object
{
  out <- lvm.ob$data$n
  return(out)
} #required to make IC functions work



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


Km1m2 <- function(m1,p1,X1,m2,p2,X2,rspns=Y,n=10000)
{
  Y1 <- sim(x = m1,n = n,p = p1,X = X1)$Y 
  Y2 <- sim(x = m2,n = n,p = p2,X = X2)$Y
  K12 <- mkKPQ(obP = Y1,obQ = Y2)$KPQ
  out <- list(K12=K12)
  return(out)
}


clbrt_btstrp <- function #calibrate existing double bootstrap
(t0,t1,t2, cnf_pnts=c(0.05,0.95), bnd_type=c("L","U"), cnsrvtv=T,note="") # arguments
{ # function body
  { # clbrt_btstrp header documentation
    # t0: scalar theta_hat from original data 
    # t1: vector of level 1 theta_hat
    # t2: matrix level 2 theta_hat. Columns of t2 should equal length of t1. 
    #   : A t2 column should be estimates from bootstraps of a single level 1 bootstrap.
    #   : It may be good if dim(t2)[1] = length(t1) but it is not necessary,
    # cnf_pnts: vector of length 1 or 2.  these are the confidence bounds
    # bnd_type: vector of bound types "L"=lower, "U"=upper. length must equal length cnf_pts 
    #         : If cnf_pnt*dim(t2)[1] is not a "whole number" indices adjusted to make desired bnd conservative
    #         : L will use floor and U will use ceiling
    # Version 2: include conservative/liberal switch
  } # end of header documentation
  { # functions
    sgn_chngd <- function(x1,x2){sign(x1) != sign(x2)  }
    
    what_cvr <- function(ibnd,t0=t0,t2srt=t2srt){
      # proportion t1 <= t2 quantile
      IB <- TRUE
      dm <- dim(t2srt)[1]
      if (ibnd <= 1 ) {inbd <-1
      IB <- FALSE
      } else if (ibnd  >= dm) {
        ibnd <- dm
        IB <- FALSE
      }
      bnds <- t2srt[ibnd,]                       
      {cvr <- mean((t0 <= bnds ))} 
      attr(cvr,"internal_bound") <- IB
      return(cvr)
    } #end what_cvr
    

    limit_ci <- function(ci,dm) #keeps ci within [1:dm]
    {
      {# limits index and sets "internal_bound" attr
      } # end header documentatio
      if ((ci <= 1) ) {
        ci <- 1
        attr(ci,"internal_bound") <- FALSE
      }
      if ((ci >= dm) ) {
        ci <- dm
        attr(ci,"internal_bound") <- FALSE
      }
      return(ci)
    }
    
    what_i4p <- function(p,bnd_type,dm){
      {# If p=0 or 1 the attribute "internal_bound" is set to FALSE
      } # what_i4p header documentation
      idx <- switch(EXPR = bnd_type,
                    L = floor(p*dm),
                    U = ceiling(p*dm))
      attr(idx,"internal_bound") <- T
      if (idx <= 1) {idx <- 1
      attr(idx,"internal_bound") <- F}
      if (idx >= (dm)) {idx <- dm
      attr(idx,"internal_bound") <- F}
      if (!is.null(attr(p,"internal_bound"))) {
        if (attr(p,"internal_bound") == FALSE) {attr(idx,"internal_bound") <- FALSE}
      }
      return(idx)
    } # end what_i4p
    
    what_p4i <- function(i,dm)
    {
      p <- (i/dm) 
      attr(p,"internal_bound") <-TRUE
      if (p<=(1/dm)) 
      {p <- 1/dm
      attr(p,"internal_bound") <- FALSE
      }
      if (p >= 1)
      {P <- 1
      attr(p, "internal_bound") <- FALSE
      }
      return(p)
    } # end what_p4i
    
    what_incrmnt <-function (av_cvr,pbnd) # find change of index
    {
      # All increments calculated as lower bounds
      # if coverage (left tail) is low, look a higher index (increment = 1)
      # if coverage correct, no change is needed (increment = 0)
      # if coverage is high increment = -1
      incrmnt <- sign( pbnd - av_cvr )
      return(incrmnt)    
    } # end what_incrmnt
    
    what_fnl_incrmnt <- function #adjust final index to maintain a conservative bnd
    (incrmnt=incrmnt, bnd_type=bnd_type, cnsrvtv=cnsrvtv) # arguments
    {
      if(incrmnt==0) {fi <- 0}
      if (cnsrvtv) {
        if(incrmnt==-1) {
          if(bnd_type=="L") {fi <- 0
          } else { fi <- 1}
        }
        if(incrmnt==1) {
          if(bnd_type=="L") {fi <- -1
          } else { fi <- 0}
        }
      } else {  # cnsrvtv == F
        if(incrmnt==-1) {
          if(bnd_type=="L") {fi <- 1
          } else { fi <- 0}
        }
        if(incrmnt==1) {
          if(bnd_type=="L") {fi <- 0
          } else { fi <- -1}
        }
      }
      return(fi)
    } # end  what_fnl_incrmnt
    

    clbrt_p <- function(nom_p, bnd_type, t0, t2srt,cnsrvtv=T)# find calibrated p 
    {
      {# nom_p: nominal confidence point  
        # t1: vector of level 1 theta each estimate from a level 1 bootstrap; 
        # t2srt: matrix of level 2 bootstrap theta, each col from a level1 bootstrap.
        #      : cols of t2srt correspond to elements of t1  
        # bnd_type: either "L" or "U". adjusts final indexes so conservative
        # The root finding is done iteratively.  Integer steps away from the initial 
        # guess of the percentile CI index until a sign change in the objective function.
        # More efficiency might be achieved with larger initial step sizes.
      }# end documentation
      t2rws <- dim(t2srt)[2]
      t2cls <- dim(t2srt)[1]
      isunsrtd.vec <- aaply(.data=t2srt,.margins = 2,.fun = is.unsorted)
      if (any(isunsrtd.vec)) # check is t2srt has any unsorted columns
      {
        for (j in 1:t2cls) {if (isunsrtd.vec[j]) t2srt[,j] <- sort(t2srt[,j])}
      }
      
      nom_i <- what_i4p(p = nom_p,dm = t2rws,bnd_type) #index corresponding to nominal P
      ci <- nom_i # initialize calibrated index
      cvr_nom_p <- what_cvr(ibnd = nom_i, t0 = t0, t2srt = t2srt) 
      cvr <- cvr_nom_p
      old_incrmnt <- what_incrmnt(av_cvr = cvr_nom_p, pbnd = nom_p)
      if (cvr_nom_p == nom_p) {good_as_it_gets <- T
      } else {good_as_it_gets <- F}
      while ( !good_as_it_gets) # if coverage not exact will try at least once
      {
        incrmnt <- what_incrmnt(av_cvr = cvr,pbnd = nom_p)
        ci <- ci + incrmnt
        ci <- limit_ci(ci,dm=t2rws) # keep ci, between 1 and t2rws
        cvr <- what_cvr(ibnd = ci, t0=t0, t2srt=t2srt)
        if ((cvr == nom_p) | (attr(cvr,"internal_bound") == FALSE)) {good_as_it_gets <- T}
        if (sgn_chngd(old_incrmnt,incrmnt)) {
          good_as_it_gets <- T
          fnl_incrmnt <- what_fnl_incrmnt(incrmnt=incrmnt, bnd_type=bnd_type,cnsrvtv = cnsrvtv )
          ci <- ci + fnl_incrmnt # final increment keeps desired interval conservative
          ci <- limit_ci(ci,dm=t2rws)
        }
        old_incrmnt <- incrmnt
      } # end while
      clbrtd_p <- what_p4i(i = ci,dm = t2rws)
      return(clbrtd_p)
    }#end clbrt_p()
  } # end functions
  {# initializations
    version=1
    num_bnds <- length(cnf_pnts)
    t1len=length(t1)
    t2rws=dim(t2)[1]
    {run <- vector(mode = "list", length = 7)
      names(run) <- c("note","Version","start","stop","elapsed","top_B","bttm_B")
      run$note <- note
      run$start <- Sys.time()
      run$version <- version
      run$top_B <- t1len
      run$bttm_B <- t2rws
    } # run info
    t1srt <- sort(x = t1,decreasing = F) #used for percentile intevals. Note not t1!
    t2srt <- t(aaply(.data = t2,.margins = 2,.fun = sort)) #sorts within each column
    #despite the absence of decreasing=F in the sort statment t2srt is properly sorted
    rownames(t2srt) <- NULL  # each column is sorted independently any rownames are now meaningless
    bnd.df <- as.data.frame(matrix(NA, nrow=num_bnds, ncol=6))
    bnd.df[,5] <- rep(TRUE,num_bnds) # will be set to false if bound not internal
    bnd.df[,6] <- rep(TRUE,num_bnds) # will be set to false if bound not internal
    colnames(bnd.df) <- c("bnd_type","cnf_pt","%bnd","clbrtd_bnd","%bnd_IB","clbrtd_bnd_IB")
  } # end initializations
  { # main loop
    all_bounds_internal <- TRUE
    for (i in 1:num_bnds) {
      bnd.df[i,1] <- bnd_type[i]
      bnd.df[i,2] <- cnf_pnts[i]
      ibnd <- what_i4p(p = cnf_pnts[i],bnd_type = bnd_type[i],dm = t1len)
      bnd.df[i,3] <- t1srt[ibnd]
      if (attr(ibnd,"internal_bound") == FALSE) {
        bnd.df[i,5]  <- FALSE
        all_bounds_internal <- FALSE
      }
      cp <- clbrt_p(nom_p = cnf_pnts[i],bnd_type = bnd_type[i],t0=t0,t2srt = t2srt)
      ibnd <- what_i4p(p = cp,bnd_type = bnd_type[i],dm=t1len)
      bnd.df[i,4] <- t1srt[ibnd]
      if (attr(ibnd,"internal_bound") == FALSE) {
        bnd.df[i,6] <- FALSE
        all_bounds_internal <- FALSE
      }
    }
  } # end of main loop
  {
    run$stop <- Sys.time()
    run$elapsed <- run$stop - run$start
    out <- list(bnd.df=bnd.df, all_bounds_interal=all_bounds_internal, run=run)
  } # end of clean up
  return(out)
} # end clbrt_btstrp 



est_M_clbrtd <- function(data0,fixdcov=NULL,mdnms,cmp_lst,evType="BIC",top_B=100,
                         bttm_B=100, cnsrvtv=T,sv_ev_bt.mat=F,sv_Bidx=F,
                         ev_cnf_bnds=c(0.05, 0.5, 0.5, 0.95),
                         ev_cnf_bnd_type=c("L","L","U","U"),
                         Rl_cnf_bnds=c(0.02, 0.05, 0.5, 0.95, 0.98 ), 
                         prlll=T,note="",talkback=T)
{
  # estimate M' and confidence limits
  # output 1) list of fitted models 2) summary table 3) list of run info 4)optional Matrix of ev comparisons.
  # 5) optional matrix of indicies for top level bootstrap
  # Based on est_M version 2.  Foreach combine type changed to list to allow including full 2nd level bootstrap evmat.
  # The access to bottom bootstrap at higher level allows bootstrap calibration
  # 2019-02-03 23:07:03 EST Version 2: cnsrvtv switch toggles between conservative and liberal bounds
  #2019-02-10 20:52:52 MST Version 3: output corrected to include md_aM and aM
  # 2019-03-06 16:43:00 MST Version 4: switching output from aM to Rl (Rl == 1-aM)
  {# initializing
    Version <- 4
    #startup bookkeeping
    {run <- vector(mode = "list", length = 9)
      names(run) <- c("note","Version","start","stop","elapsed",
                      "top_B","bttm_B","smpl_sz","cnsrvtv")
      run$note <- note
      run$start <- Sys.time()
      run$Version <- Version
      run$top_B <- top_B
      run$bttm_B <- bttm_B
      run$smpl_sz <- dim(data0)[1]
      run$cnsrvtv <- cnsrvtv
    }
    
    ec <- get(evType) # defining evidence criterion function
    mdlst <- plyr::alply(.data=mdnms,.margins=1,.fun=get) # creates a list of models from a list of names
    names(mdlst) <- mdnms #names list elements
    nummds <- length(mdlst)# number of models in list
    fitlst <- vector(mode = "list",length = nummds)
    names(fitlst) <- mdnms
    for (i in 1:nummds) {
      fitlst[[i]] <- lava::estimate(x = mdlst[[i]],data = data0)
    }
    ec_vec <- plyr::laply(.data = fitlst,.fun = ec)
    names(ec_vec) <- mdnms
    ec_mat <- matrix(ec_vec,nrow=1)
    colnames(ec_mat) <- mdnms
    ev0 <- (Ev4(ecmat = ec_mat,cmp_lst = cmp_lst))
    
    if (prlll) {
      numwrkrs <- getDoParWorkers() # number of available cores for ||processing
    } else {
      numwrkrs <- 1
    }
    ss <- dim(data0)[1]
    Bidx_top <- balsmpl(B=top_B,sz = ss, incl_bt0 = T) #indecies for creating balanced bootstrap
    Blki_top <- vector("list",length = numwrkrs) #list to hold blocks of Bidx
    cuts <- CutBySize(top_B,nb=numwrkrs) # indexes to break Bids into blocks
    for (i in 1:numwrkrs) {
      Blki_top[[i]] <- matrix(Bidx_top[cuts[i,1]:cuts[i,2],],ncol=ss,byrow=T) #matrix create correct form when R=1
      rownames(Blki_top[[i]]) <- rownames(Bidx_top)[cuts[i,1]:cuts[i,2]]  #keeps proper rownames even if R=1
    }
  }# end initializing
  if (prlll) # parallel processing
  { ev_bt.lst <- foreach(i = 1:numwrkrs, .inorder = F, .combine=c, #because blk results are list of list, this just makes longer list of list
                         .packages = c("foreach","lava","stringr","plyr","mvtnorm","stats4"),
                         .export = c("fes","fev","BIC","KICS", "KKIC","balsmpl","HIC","nobs.lvmfit",
                                     "Ev4md","Ev4","CutBySize","m10","m11","m12","m13",
                                     "m14","m15","m16","m17","m18","m19",
                                     "m20","m25","data0", "NPPboot","fltn_ary.lst"),
                         .verbose = talkback) %dopar% {
                           #                          sink("log.txt", append=TRUE)
#                           cat("Starting iteration",i,"\n")
                           nrwBlk <- dim(Blki_top[[i]])[1] #number of rows in block
                           Blk_rslts <- foreach (j = 1:nrwBlk,  #default .combine is list so .combine omitted
                                                 .errorhandling = 'remove') %do% {
                                                   idx=Blki_top[[i]][j,]
                                                   NPd <- data0[idx,]
                                                   # NPecvec <- foreach (k = 1:nummds, .combine = "c") %do% {
                                                   #   ec(lava::estimate(mdlst[[k]],data=NPd)) #create rows of ec
                                                   # }
                                                   NPec <- NPPboot(data = NPd,mdnms = mdnms,R=bttm_B,
                                                                   evType = evType,NPbt = T,
                                                                   Pbt = F,bt0 = T,prlll = F)$NPec
                                                   evmat <- Ev4(ecmat = NPec,cmp_lst = cmp_lst)
                                                   ev1 <- evmat[1,]
                                                   names(ev1) <- colnames(evmat)
                                                   mn_ev <- plyr::aaply(.data = evmat,.margins = 2,.fun = mean)
                                                   names(mn_ev) <- colnames(evmat)
                                                   aM <- aaply(.data = evmat,.margins = 2,.fun = function(X){mean(X <= 0)})
                                                   names(aM) <- colnames(evmat)
                                                   top_bt <- rownames(Blki_top[[i]])[j]
                                                   out <- list(mn_ev=mn_ev,  aM=aM, ev1=ev1, evmat=evmat,top_bt=top_bt)
                                                   return(out)
                                                 }
                           return(Blk_rslts) # combine blocks into overall list
                           
                         }
  
  } else 
  { # sequential processing
    ev_bt.lst <- foreach(i = 1:numwrkrs, .inorder = F, .combine=c, #because blk results are list of list, this just makes longer list of list
                         .packages = c("foreach","lava","stringr","plyr","mvtnorm","stats4"),
                         .export = c("fes","fev","BIC","KICS", "KKIC","balsmpl","HIC","nobs.lvmfit",
                                     "Ev4md","CutBySize","m10","m11","m12","m13",
                                     "m14","m15","m16","m17","m18","m19",
                                     "m20","data0", "NPPboot","fltn_ary.lst"),
                         .verbose = talkback) %do% {
                           nrwBlk <- dim(Blki_top[[i]])[1] #number of rows in block
                           Blk_rslts <- foreach (j = 1:nrwBlk,  #default .combine is list so .combine omitted
                                                 .errorhandling = 'remove') %do% {
                                                   idx=Blki_top[[i]][j,]
                                                   NPd <- data0[idx,]
                                                   NPecvec <- foreach (k = 1:nummds, .combine = "c") %do% {
                                                     ec(lava::estimate(mdlst[[k]],data=NPd)) #create rows of ec
                                                   }
                                                   NPec <- NPPboot(data = NPd,mdnms = mdnms,R=bttm_B,
                                                                   evType = evType,NPbt = T,
                                                                   Pbt = F,bt0 = T,prlll = F)$NPec
                                                   evmat <- Ev4(ecmat = NPec,cmp_lst = cmp_lst)
                                                   ev1 <- evmat[1,]
                                                   names(ev1) <- colnames(evmat)
                                                   mn_ev <- plyr::aaply(.data = evmat,.margins = 2,.fun = mean)
                                                   names(mn_ev) <- colnames(evmat)
                                                   aM <- aaply(.data = evmat,.margins = 2,.fun = function(X){mean(X <= 0)})
                                                   names(aM) <- colnames(evmat)
                                                   top_bt <- rownames(Blki_top[[i]])[j]
                                                   out <- list(mn_ev=mn_ev,  aM=aM, ev1=ev1, evmat=evmat,top_bt=top_bt)
                                                   return(out)
                                                 }
                           return(Blk_rslts) # combine blocks into overall list
                         }
    
  } # end of if else
  saveRDS(ev_bt.lst,file = "ev_bt.lst.rds") #use for debuging if crash in clean up
  evmat <- NULL
  for (i in 1:length(ev_bt.lst)){
    evmat <- rbind(evmat,ev_bt.lst[[i]]$evmat[1,])
  }
  #ev_bt.mat <- fltn_ary.lst(ev_bt.lst)
  ncmp <- dim(evmat)[2]
  mn_ev <- NULL
  for (i in 1:top_B){
    mn_ev <- rbind(mn_ev,ev_bt.lst[[i]]$mn_ev)
  }
  grnd_mn_ev<- plyr::aaply(.data = mn_ev[,1:ncmp],.margins = 2,.fun = mean)
  Rl <- NULL
  for (i in 1:top_B){
    Rl <- rbind(Rl,1 - ev_bt.lst[[i]]$aM)
  }
  mn_Rl<- plyr::aaply(.data = Rl[,1:ncmp],.margins = 2,.fun = mean)
  md_Rl<- plyr::aaply(.data = Rl[,1:ncmp],.margins = 2,.fun = median)
  Rl_bnds=t(plyr::aaply(.data = Rl,.margins = 2,.fun = quantile,probs=Rl_cnf_bnds))
  rownames(Rl_bnds) = str_c(rownames(Rl_bnds),"_bnd",sep="")
  # bca_Rl_bnds <- matrix(rep(NA,ncmp*length(Rl_cnf_bnds)),ncol=ncmp)
  # colnames(bca_Rl_bnds) <- colnames(Rl_bnds)
  # rownames(bca_Rl_bnds) <- str_c(rownames(Rl_bnds),".BCa",sep="")
  # for (i in 1:length(Rl_cnf_bnds)) {
  #   cnf <- 1-((1 - Rl_cnf_bnds[i])*2)
  #   for (j in 1:ncmp) {
  #     bca_Rl_bnds[i,j] <- bootBCa(estimate = Rl[1,j],estimates = Rl[,j],type = "bca",conf.int = cnf,n = top_B,seed = .Random.seed)[2]
  #   }
  # }
  ev_bnds <-matrix(NA,nrow = 2*length(ev_cnf_bnds),ncol=ncmp) #matrix to hold conf limits on evidence
  colnames(ev_bnds) <- colnames(Rl)
  rownames(ev_bnds) <- str_c(ev_cnf_bnds,ev_cnf_bnd_type,rep(c("prcnt","clbrtd"),each=length(ev_cnf_bnds)))
  intrnl_bnds <-matrix(NA,nrow = 2*length(ev_cnf_bnds),ncol=ncmp) #matrix to hold conf limits on evidence
  colnames(intrnl_bnds) <- colnames(Rl)
  rownames(intrnl_bnds) <- str_c(ev_cnf_bnds,ev_cnf_bnd_type,rep(c("prcnt","clbrtd"),each=length(ev_cnf_bnds)))
  for (i in 1:ncmp){
    t1 <- evmat[,i]
    t2 <- NULL 
    for (j in 1:top_B){
      t2 <- cbind(t2,ev_bt.lst[[j]]$evmat[,i])
    }
    bnd.df <- clbrt_btstrp(t0 = ev0[i],t1 = t1,t2 = t2, cnf_pnts = ev_cnf_bnds, bnd_type = ev_cnf_bnd_type,cnsrvtv = cnsrvtv)$bnd.df
    bnd.vec <-c(bnd.df[,"%bnd"], bnd.df[,"clbrtd_bnd"])
    IB.vec <-c(bnd.df[,"%bnd_IB"], bnd.df[,"clbrtd_bnd_IB"])
    ev_bnds[,i] <- bnd.vec
    intrnl_bnds[,i] <- IB.vec
  }
  
  smmry <- list(ev0=ev0, mn_ev=mn_ev, ev_bnds=ev_bnds, intrnl_bnds=intrnl_bnds, 
                mn_Rl=mn_Rl, md_Rl=md_Rl, Rl_bnds=Rl_bnds, Rl=Rl)#, bca_Rl_bnds=bca_Rl_bnds)
  if (!sv_ev_bt.mat) ev_bt.mat <- NULL
  if (!sv_Bidx)  Bidx_top <- NULL
  run$stop <- Sys.time()
  run$elapsed <- run$stop - run$start
  out <- list(fttd_mds=fitlst,ev_summry=smmry,run=run,ev_bt.mat=ev_bt.mat,Bidx_top=Bidx_top)
  return(out)
}# end est_M_clbrtd 



NPPboot <- function(data=NULL,gennm=NULL,genpar=NULL, mdnms, simnm=NULL, simpar=NULL, 
                    simsigma=0, fixdcov=NULL, ss=-1, R, 
                    evType="BIC",note="",
                    NPbt=T, Pbt=T, bt0=T, prlll=T)
{
  #NPPboot generates model comparisons w/ nonparametric & parametric bootstraps.
  #data=dataset, if data is NULL, initial data are generated.
  #gennm and genpar are the name and parameters used to generate initial data.
  #gennm and genpar can be passed in to identify a data set.
  #mdnms=vector of names of models to be compared. 
  #simmd=model to simulate from.
  #simpar=parameter vector for simulations. 
  #simsigma is the vcov for incorparating estimation error in parameters.
  #by default simsigma is a zero matrix resulting in no variation.
  #fixdcov is a matrix of covariates to be fixed in data. 
  #fixdcov names should match variables in data.
  #ss= sample size for both NP and P bootstraps.
  #If ss<1 ss is dim(data)[1]. If no data given, ss must be specified. 
  #R=number of bootstraps. note=descriptive note.
  #evType is the name, as a string, of an evidence criterion (differences are evidence)
  #The boolean parameters NPbt and Pbt allow you to switch off either the NonParametric bootstrap
  #or the Parametric bootstrap. The boolean argument bt0 controls whether original estimates returned as part of bootstrap
  # The boolean argument prlll controls whether the calculations run serially or parallelly.
  #Version 1:F 2018-11-26 14:09. Including run time and version number
  #Version 2: 2018-12-01 13:37:52 MST. Switiching to input of model names rather than models
  #Version 3: 2018-12-03 15:38:32 MST. simulation of initial data now an option
  #Version 4: 2018-12-04 15:12:21 MST. Reogranized output
  #Version 5: 2018-12-21 08:00:23 MST. Added support for KKIC as a evType 
  #Version 6: 2018-12-24 21:37:17 MST. Fixed bug where covariate missing in Pboot
  #Version 7: 2018-12-30 19:26:58 MST. Putting parallel proc switch
  #Version 8: 2019-01-03 10:08:25 MST. R=1 now doesn't cause crash
  #Version 9: 2019-01-06 20:19:32 MST. revised balsmpl includes original estimate in first posisiton  
  #Version 10: 2019-01-08 09:59:39 MST. moved some stuff only required for sim inside Pbt=T {}
  #                                     Other changes to make run when PbT = false        
  # Initializing
  {
    Version <- 10
    
    
    
    #startup bookkeeping  
    {run <- vector(mode = "list", length = 5)
      names(run) <- c("note","Version","start","stop","elapsed")
      run$note <- note
      run$start <- Sys.time()
      run$Version <- Version
    }
    
    if (is.null(data)) # Simulate data if none given
    { 
      if (is.null(gennm)) stop("genmd is required if data is generated")
      genmd <- get(gennm)
      if (is.null(genpar)) stop("genpar is required if data is generated")
      if (ss<1) stop("ss is required if data is generated")
      data=sim(x = genmd,n =ss,p=genpar,X=fixdcov)
    }
    
    if (ss < 1) ss <- dim(data)[1]
    ec <- get(evType) # defining evidence criterion function
    mdlst <- plyr::alply(.data=mdnms,.margins=1,.fun=get) # creates a list of models from a list of names
    names(mdlst) <- mdnms #names list elements
    nummds=length(mdlst)# number of models in list
    if (prlll) {
      numwrkrs <- getDoParWorkers() # number of available cores for ||processing
    } else {
      numwrkrs <- 1
    }
    Bidx <- balsmpl(B=R,sz = ss, incl_bt0 = T) #indecies for creating balanced bootstrap
    Blki <- vector("list",length = numwrkrs) #list to hold blocks of Bidx
    cuts <- CutBySize(R,nb=numwrkrs) # indexes to break Bids into blocks
    for (i in 1:numwrkrs) {
      Blki[[i]] <- matrix(Bidx[cuts[i,1]:cuts[i,2],],ncol=ss,byrow=T) #matrix create correct form when R=1
    }
  }# end of initializing
  
  if (NPbt) # Nonparametric bootstraping if requested
  {
    if (prlll)
    {NPec <- foreach(i = 1:numwrkrs, .combine=rbind, .inorder = F,
                     .packages = c("foreach","lava","stringr"),
                     .export = c("fes","fev","KICS", "KKIC","balsmpl","HIC","nobs.lvmfit"),
                     .verbose = F) %dopar% {
                       nrwBlk <- dim(Blki[[i]])[1] #number of rows in block
                       tmp <- foreach (j = 1:nrwBlk, .combine=rbind,
                                       .errorhandling = 'remove') %do% {
                                         idx=Blki[[i]][j,]
                                         NPd <- data[idx,]
                                         NPecvec <-foreach (k = 1:nummds, .combine = "c") %do% {
                                           ec(lava::estimate(mdlst[[k]],data=NPd)) #create rows of ec
                                         }
                                         NPecvec #combine rows into blocks
                                       }
                       tmp # combine blocks into matrix of ec
                     }
    } else {
      NPec <- foreach(i = 1:numwrkrs, .combine=rbind, .inorder = F,
                      .packages = c("foreach","lava","stringr"),
                      .export = c("fes","fev","KICS", "KKIC","balsmpl","HIC","nobs.lvmfit"),
                      .verbose = F) %do% {
                        nrwBlk <- dim(Blki[[i]])[1] #number of rows in block
                        if (i==0) { tmp <- NULL
                        } else {
                          tmp <- foreach (j = 1:nrwBlk, .combine=rbind,
                                          .errorhandling = 'remove') %do% {
                                            idx=Blki[[i]][j,]
                                            NPd <- data[idx,]
                                            NPecvec <- foreach (k = 1:nummds, .combine = "c") %do% {
                                              ec(lava::estimate(mdlst[[k]],data=NPd)) #create rows of ec
                                            }
                                            NPecvec #combine rows into blocks
                                          }
                        }
                        tmp # combine blocks into matrix of ec
                      }
    }
  } else NPec <- NULL # if NP bootstrapping not requested
  
  if (Pbt) # parametric bootstrapping if requested
  {
    simmd <- get(simnm) # pulling in the model associated with the name
    
    if (is.null(simpar))#need to fit simmd if no simpar
    {
      esimmd <- estimate(x = simmd,data = data)
      simpar <- coef(esimmd)
      simnm <- str_c("e",simnm,sep="") #name change to indicate internal estimation
      if (simsigma==0) {
        simsigma <- diag(rep(0,length(simpar)))
      }else {
        simsigma <- vcov(esimmd)
        simnm <- str_c(simnm,"v",sep="") #name change indicates using vcov in pbt
      }
    }
    if (simsigma==0) simsigma <- diag(rep(0,length(simpar))) # zeroing simsigma 
    
    if (prlll) {
      Pec <- foreach(i = 1:numwrkrs, .combine=rbind, .inorder = F,
                     .packages = c("foreach","lava","mvtnorm","stringr"),
                     .export = c("fes","fev","KICS", "KKIC","balsmpl","HIC","nobs.lvmfit")
      ) %dopar% {
        nrwBlk <- dim(Blki[[i]])[1] #number of rows in block
        tmp <-foreach (j = 1:nrwBlk, .combine = rbind,
                       .errorhandling = 'remove') %do% {
                         idx=Blki[[i]][j,]
                         if (is.null(fixdcov)) {
                           X <- NULL
                         } else {
                           X <- fixdcov[idx,] #selects same fixed cov as NP bootstrap
                         }
                         simp <- as.vector(mvtnorm::rmvnorm(n=1,mean = simpar,sigma = simsigma))
                         #note if simsigma==0 simp is unchanged simpar
                         names(simp) <- names(simpar)
                         Pd <- lava::sim(simmd,n=ss,p=simp,X=X) # generating data w/ Parametric bootstrap
                         misscov <- setdiff(colnames(fixdcov),colnames(Pd))
                         Pd <- cbind(Pd,X[,misscov]) #adds in covariate not used in simmd
                         Pecvec <- foreach (k = 1:nummds, .combine = "c" ) %do% {
                           ec(lava::estimate(mdlst[[k]],data=Pd)) #create rows of ec
                         }
                         Pecvec #combine rows into blocks
                       }
        tmp # combine blocks into matrix of ec
      }
      
    } else {
      Pec <- foreach(i = 1:numwrkrs, .combine=rbind, .inorder = F,
                     .packages = c("foreach","lava","mvtnorm","stringr"),
                     .export = c("fes","fev","KICS", "KKIC","balsmpl","HIC","nobs.lvmfit")
      ) %do% {
        nrwBlk <- dim(Blki[[i]])[1] #number of rows in block
        tmp <-foreach (j = 1:nrwBlk, .combine = rbind,
                       .errorhandling = 'remove') %do% {
                         idx=Blki[[i]][j,]
                         if (is.null(fixdcov)) {
                           X <- NULL
                         } else {
                           X <- fixdcov[idx,] #selects same fixed cov as NP bootstrap
                         }
                         simp <- as.vector(mvtnorm::rmvnorm(n=1,mean = simpar,sigma = simsigma))
                         #note if simsigma==0 simp is unchanged simpar
                         names(simp) <- names(simpar)
                         Pd <- lava::sim(simmd,n=ss,p=simp,X=X) # generating data w/ Parametric bootstrap
                         misscov <- setdiff(colnames(fixdcov),colnames(Pd))
                         Pd <- cbind(Pd,X[,misscov]) #adds in covariate not used in simmd
                         Pecvec <- foreach (k = 1:nummds, .combine = "c" ) %do% {
                           ec(lava::estimate(mdlst[[k]],data=Pd)) #create rows of ec
                         }
                         Pecvec #combine rows into blocks
                       }
        tmp # combine blocks into matrix of ec
      }
    }
  } else Pec <- NULL # if parametric bootstrap not requested
  
  # Post run clean up
  {
    if (!is.null(NPec)) #setting names and attributes 
    {
      if (is.vector(NPec)) NPec <- t(NPec) #make sure NPec is matrix
      colnames(NPec) <- names(mdlst)
      {attr(NPec,"evType") <- evType}
    }
    if (!is.null(Pec))  #setting names and attributes 
    {
      if (is.vector(Pec)) Pec <- t(Pec) #make sure NPec is matrix
      colnames(Pec) <- names(mdlst)
      attr(Pec,"evType") <- evType
    }
    
    run$stop <- Sys.time()
    run$elapsed <- run$stop - run$start
    
    if (is.null(gennm)) {genmd <- NULL}
    if (is.null(simnm)) {simmd <- NULL}
    
    out=list(NPec=NPec,Pec=Pec,data=data,
             gennm=gennm,genmd=genmd,genpar=genpar,
             mdnms=mdnms,mdlst=mdlst,
             simnm=simnm,simmd=simmd,simpar=simpar,
             simsigma=simsigma,fixdcov=fixdcov,ss=ss,R=R,evType=evType,
             run=run,NPbt=NPbt,Pbt=Pbt,bt0=bt0,prlll=prlll)
    
  } #end of clean up
  
  return(out) #return results
} # end of NPPboot




fltn_ary.lst <- function(ary.lst) # Flattens a list of arrays into a single array
{# retains rownames, uses colnames from first array
  out <- ary.lst[[1]]
  if (length(ary.lst)>1){
    for (i in 2:length(ary.lst)){
      out <- rbind(out,ary.lst[[i]])
    }
  }
  return(out)
} # end of fltn_ary.lst






conv_cntrst.nm <- function(cntrst.nm){
  cn <- cntrst.nm
  out <- str_c(alt_nm[str_sub(string = cn,start = 1,end = 3)]," ",
               alt_nm[str_sub(string = cn,start = 5,end = 8)]) 
  return(out)
}




###Utility functions
balsmpl=function(B,sz,mltp=1,incl_bt0=T,rwname="bt"){
  #Returns a BXsz matrix of indices for a balanced bootstrap.  
  #All indicies occur = number of times in matrix  
  #Balanced bootstraps are more efficient than regular bootstraps
  #Use each row of the matrix to select a bootstrapped sample
  #Each row (i.e. bootstrap) is given an identifier in the rowname e.g. "bt324"
  #If incl_b0=T then first row is 1:sz (i.e original estimate)
  idxes=sample(rep(rep(1:sz,each=mltp),times=B))
  out=matrix(idxes,nrow=B)
  rownames(out) <- stringr::str_c(rwname,1:B,sep = "")
  if (incl_bt0) {
    out <- rbind(1:sz,out)
    rownames(out) <- stringr::str_c(rwname,0:B,sep = "")
  }
  return=out
} #makes indecies for a balanced bootstrap


#######################################################################################################
##models

m13nms=c( "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m25")

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



m10 <- lvm()#(Y~X+Y~W+Y~U+Y~V+Y~Z+Y~ZZ) no latent y, no latent predictors, predictors not random
{
  regression(m10) <- list(
    Y~X,
    Y~W,
    Y~U,
    Y~V,
    Y~Z,
    Y~ZZ)
}

m11 <- lvm() #m10-(Y~U)
{
  regression(m11) <- list(
    Y~X,
    Y~W,
    Y~V,
    Y~Z,
    Y~ZZ)
}

m12 <- lvm()#m10 -(Y~V)
{
  regression(m12) <- list(
    Y~X,
    Y~W,
    Y~U,
    Y~Z,
    Y~ZZ)
}



m13 <- lvm() #m10 - (Y~U + Y~V)
{
  regression(m13) <- list(
    Y~X,
    Y~W,
    Y~Z,
    Y~ZZ)
}

m14 <- lvm()#m10 - (Y~X)
{
  regression(m14) <- list(
    Y~W,
    Y~U,
    Y~V,
    Y~Z,
    Y~ZZ)
}


m15 <- lvm()#m10 - (Y~W)
{
  regression(m15) <- list(
    Y~X,
    Y~U,
    Y~V,
    Y~Z,
    Y~ZZ
  )
}

m16 <- lvm()#m10 - (Y~Z + Y~ZZ)
{
  regression(m16) <- list(
    Y~X,
    Y~W,
    Y~U,
    Y~V
  )
}

m17 <- lvm()#m10 - (Y~ZZ)
{
  regression(m17) <- list(
    Y~X,
    Y~W,
    Y~U,
    Y~V,
    Y~Z
  )
}


m18 <- lvm() #m10 - (Y~U + Y~V + Y~W)
{
  regression(m18) <- list(
    Y~X,
    Y~Z,
    Y~ZZ)
}

m19 <- lvm() #m10 - (Y~U + Y~V + Y~ZZ)
{
  regression(m19) <- list(
    Y~X,
    Y~W,
    Y~Z)
}

m20 <- lvm() #m10 - (Y~U + Y~V + Y~W +ZZ)
{
  regression(m20) <- list(
    Y~X,
    Y~Z)
}


m21 <- lvm() #m10-(Y~X + Y~U)
{
  regression(m21) <- list(
    Y~W,
    Y~V,
    Y~Z,
    Y~ZZ)
}

m22 <- lvm()#m10 -(Y~X + Y~V)
{
  regression(m22) <- list(
    Y~W,
    Y~U,
    Y~Z,
    Y~ZZ)
}



m23 <- lvm() #m10 - (Y~X, Y~U + Y~V)
{
  regression(m23) <- list(
    Y~W,
    Y~Z,
    Y~ZZ)
}

m24 <- lvm()#m10 - (Y~W + Y~X)
{
  regression(m24) <- list(
    Y~U,
    Y~V,
    Y~Z,
    Y~ZZ)
}

m25 <- lvm()#m10 - (Y~W + Y~ZZ)
{
  regression(m25) <- list(
    Y~X,
    Y~U,
    Y~V,
    Y~Z)
}


m26 <- lvm()#m10 - (Y~X + Y~Z + Y~ZZ)
{
  regression(m26) <- list(
    Y~W,
    Y~U,
    Y~V
  )
}

m27 <- lvm()#m10 - (Y~X + Y~ZZ)
{
  regression(m27) <- list(
    Y~W,
    Y~U,
    Y~V,
    Y~Z
  )
}

####################################################################################
#Script

x_sd <- 1  # standard deviation of latent predictor variable x
x_mn <- 0  # mean of latent predictor variable
w_sd <- 1  # standard deviation of latent predictor variable x
w_mn <- 0  # mean of latent predictor variable
Z_sd <- 1  # standard deviation of predictor variable Z
Z_mn <- 0  # mean of  predictor variable Z
s_sd <- 1  # standared deviation of spurious latent predictor
s_mn <- 0  # mean of spurious latent predictors
y_sd <- 1 # process variation
B0y <- 1   # intercept for y latent response variable
Bx <- 0.5  # slope for x latent predictor variable  
Bw <- 0.25  # slope for w latent predictor variable  
BZ <- 0.5  # slope for z predictor variable  
Bu <- 0  # slope for u spurious predictor
Bv <- 0  # slope for v spurious predictor
BZZ <- -0.15  # slope for z quadratic predictor variable  
ss <- 20000  # sample size

x=rnorm(n=ss, mean=x_mn, sd=x_sd)  # generating the predictor variable
w=rnorm(n=ss, mean=w_mn, sd=w_sd)  # generating the predictor variable
u=rnorm(n=ss, mean=s_mn, sd=s_sd)  # generating the predictor variable
v=rnorm(n=ss, mean=s_mn, sd=s_sd)  # generating the predictor variable

#Z=rnorm(n=ss, mean=Z_mn, sd=Z_sd)  # generating the predictor variable
Z=runif(n=ss, min = -4,max = 4) # generating non-normal predictor variable

ZZ <- Z^2
y =B0y + Bu*u + Bv*v +Bx*x + Bw*w + BZ*Z + BZZ*ZZ + rnorm(n=ss, mean=0, sd=y_sd)                    # calculating latent response (note no process error)

P_ME <- 0      #Measurement error variance for observed predictors x & w
Y_ME <- 0      #Measurement error variance for observed predictor Y

X <- x + rnorm(n=ss,mean=0,sd=sqrt(P_ME)) # constructing X observed with measurement error
W <- w + rnorm(n=ss,mean=0,sd=sqrt(P_ME)) # constructing W observed with measurement error
U <- u + rnorm(n=ss,mean=0,sd=sqrt(P_ME)) # constructing U observed with measurement error
V <- v + rnorm(n=ss,mean=0,sd=sqrt(P_ME)) # constructing V observed with measurement error
Y <- y + rnorm(n=ss,mean=0,sd=sqrt(Y_ME)) # constructing Y observed with measurement error

ZZZ <- cbind(Z,ZZ) # matrix of fixed covariates
covmat=cbind(X,W,U,V,Z,ZZ)

tru.g.par <- c(B0y,x_mn,w_mn,s_mn,s_mn,Bx,Bw,Bu,Bv,BZ,BZZ,y_sd^2,x_sd^2,w_sd^2,s_sd^2,s_sd^2)
names(tru.g.par) <- c("y","x","w","u","v","y~x","y~w","y~u","y~v","y~Z","y~ZZ","y~~y",
                      "x~~x","w~~w","u~~u","v~~v")
d <- data.frame(Y=Y,X=X,W=W,U=U,V=V,Z=Z,ZZ=ZZ)

g.m10.par <- c(B0y,Bx,Bw,Bu,Bv,BZ,BZZ,y_sd^2)
names(g.m10.par) <- c( "Y","Y~X","Y~W","Y~U","Y~V","Y~Z","Y~ZZ","Y~~Y")


#registerDoSEQ()
# stopCluster(cl)
# cl <- makeCluster(19)
# registerDoParallel((cl))
# unregister()
# #registerDoRNG() 
getDoParWorkers()


registerDoFuture()
plan(multiprocess)
fx20 <- covmat[1:20,]
fx50 <- covmat[1:50,]
fx100 <- covmat[1:100,]
fx200 <- covmat[1:200,]
fx400 <- covmat[1:400,]
fx2000 <- covmat[1:2000,]
fx20000 <- covmat[1:20000,]

cmp_lst <- list()
cmp_lst[[1]] <- list(alt_nms=m13nms,ref_nms="m13")
cmp_lst[[2]] <- list(alt_nms="m17",ref_nms="m10")
cmp_lst[[3]] <- list(alt_nms="m20",ref_nms="m18")
cmp_lst[[4]] <- list(alt_nms="m25",ref_nms="m15")
cmp_lst[[5]] <- list(alt_nms="m15",ref_nms="m10")
cmp_lst[[6]] <- list(alt_nms="m20",ref_nms="m19")
cmp_lst[[7]] <- list(alt_nms="m25",ref_nms="m17")





cl <- makeCluster(400, type = "MPI")  # this is supposed to define the cluster size
registerDoFuture()  
plan(cluster, workers = cl) 

B=200 # computation time increases linearly with B. B=200 desired  -- at least B=100
cB=200 # computation time increases as cB^2.  cB as big as possible, keep computation within 1 day cB=400 very nice

# version 4 redo

Sys.time()
tic("ev_cl_tinyD")
base_file_name <- "ev_cl_tinyD4"
tmp_fnm_base <- str_c("tmp_",base_file_name,sep="")
lst_nm <- str_c(base_file_name,"lst",sep=".")
sv_fl_nm <-str_c(lst_nm,"rds",sep=".")
out.lst <- vector(mode="list",length=B)

for (i in 1:B){
  d20 <- sim(x = m10,n = 20,p = g.m10.par,X=fx20)
  out.lst[[i]] <- est_M_clbrtd(data0 = d20,fixdcov = fx20,mdnms = m13nms,
                               cmp_lst = cmp_lst,
                               evType = "BIC",top_B =cB,bttm_B =cB ,
                               cnsrvtv=F, sv_ev_bt.mat = F,sv_Bidx = F,
                               ev_cnf_bnds = c(0.05,0.95),
                               ev_cnf_bnd_type=c("L","U"),
                               Rl_cnf_bnds = c(0.95,0.9,0.75,0.5),
                               prlll = T,talkback = T)
#  tmp_fnm <- str_c(tmp_fnm_base,i,"rds",sep=".")
#  saveRDS(out.lst[[i]], file=tmp_fnm)
}
assign(x = lst_nm,value = out.lst)
saveRDS(out.lst,file = sv_fl_nm)
toc()
Sys.time()

Sys.time()
tic("ev_cl_smallD4")
base_file_name <- "ev_cl_smallD4"
tmp_fnm_base <- str_c("tmp_",base_file_name,sep="")
lst_nm <- str_c(base_file_name,"lst",sep=".")
sv_fl_nm <-str_c(lst_nm,"rds",sep=".")
out.lst <- vector(mode="list",length=B)

for (i in 1:B){
  d50 <- sim(x = m10,n = 50,p = g.m10.par,X=fx50)
  out.lst[[i]] <- est_M_clbrtd(data0 = d50,fixdcov = fx50,mdnms = m13nms,
                               cmp_lst = cmp_lst,
                               evType = "BIC",top_B =cB,bttm_B =cB ,
                               cnsrvtv=F, sv_ev_bt.mat = F,sv_Bidx = F,
                               ev_cnf_bnds = c(0.05,0.95),
                               ev_cnf_bnd_type=c("L","U"),
                               Rl_cnf_bnds = c(0.95,0.9,0.75,0.5),
                               prlll = T,talkback = T)
  # tmp_fnm <- str_c(tmp_fnm_base,i,"rds",sep=".")
  # saveRDS(out.lst[[i]], file=tmp_fnm)
}
assign(x = lst_nm,value = out.lst)
saveRDS(out.lst,file = sv_fl_nm)
toc()
Sys.time()

Sys.time()
tic("ev_cl_medD4")
#B=100
base_file_name <- "ev_cl_med4"
tmp_fnm_base <- str_c("tmp_",base_file_name,sep="")
lst_nm <- str_c(base_file_name,"lst",sep=".")
sv_fl_nm <-str_c(lst_nm,"rds",sep=".")
out.lst <- vector(mode="list",length=B)

for (i in 1:B){
  d100 <- sim(x = m10,n = 100,p = g.m10.par,X=fx100)
  out.lst[[i]] <- est_M_clbrtd(data0 = d100,fixdcov = fx100,mdnms = rdcdnms,
                               cmp_lst = rdcd_lst,
                               evType = "BIC",top_B =cB,bttm_B =cB ,
                               cnsrvtv=F, sv_ev_bt.mat = F,sv_Bidx = F,
                               ev_cnf_bnds = c(0.05,0.95),
                               ev_cnf_bnd_type=c("L","U"),
                               Rl_cnf_bnds = c(0.95,0.9,0.75,0.5),
                               prlll = T,talkback = T,
                               note = "v4 Final run")
  # tmp_fnm <- str_c(tmp_fnm_base,i,"rds",sep=".")
  # saveRDS(out.lst[[i]], file=tmp_fnm)
}
assign(x = lst_nm,value = out.lst)
saveRDS(out.lst,file = sv_fl_nm)
toc()
Sys.time()

Sys.time()
tic("ev_cl_bigD4")
#B=100
base_file_name <- "ev_cl_bigD4"
tmp_fnm_base <- str_c("tmp_",base_file_name,sep="")
lst_nm <- str_c(base_file_name,"lst",sep=".")
sv_fl_nm <-str_c(lst_nm,"rds",sep=".")
out.lst <- vector(mode="list",length=B)

for (i in 1:B){
  d200 <- sim(x = m10,n = 200,p = g.m10.par,X=fx200)
  out.lst[[i]] <- est_M_clbrtd(data0 = d200,fixdcov = fx200,mdnms = rdcdnms,
                               cmp_lst = rdcd_lst,
                               evType = "BIC",top_B =cB,bttm_B =cB ,
                               cnsrvtv=F, sv_ev_bt.mat = F,sv_Bidx = F,
                               ev_cnf_bnds = c(0.05,0.95),
                               ev_cnf_bnd_type=c("L","U"),
                               Rl_cnf_bnds = c(0.95,0.9,0.75,0.5),
                               prlll = T,talkback = T,
                               note = "V4 final run")
  # tmp_fnm <- str_c(tmp_fnm_base,i,"rds",sep=".")
  # saveRDS(out.lst[[i]], file=tmp_fnm)
}
assign(x = lst_nm,value = out.lst)
saveRDS(out.lst,file = sv_fl_nm)
toc()
Sys.time()

Sys.time()
tic("ev_cl_vbigD4")
#B=100
base_file_name <- "ev_cl_vbigD4"
tmp_fnm_base <- str_c("tmp_",base_file_name,sep="")
lst_nm <- str_c(base_file_name,"lst",sep=".")
sv_fl_nm <-str_c(lst_nm,"rds",sep=".")
out.lst <- vector(mode="list",length=B)

for (i in 1:B){
  d400 <- sim(x = m10,n = 400,p = g.m10.par,X=fx400)
  out.lst[[i]] <- est_M_clbrtd(data0 = d400,fixdcov = fx400,mdnms = rdcdnms,
                               cmp_lst = rdcd_lst,
                               evType = "BIC",top_B =cB,bttm_B =cB ,
                               cnsrvtv=F, sv_ev_bt.mat = F,sv_Bidx = F,
                               ev_cnf_bnds = c(0.05,0.95),
                               ev_cnf_bnd_type=c("L","U"),
                               Rl_cnf_bnds = c(0.95,0.9,0.75,0.5),
                               prlll = T,talkback = T,
                               note = "V4 final run")
  # tmp_fnm <- str_c(tmp_fnm_base,i,"rds",sep=".")
  # saveRDS(out.lst[[i]], file=tmp_fnm)
}
assign(x = lst_nm,value = out.lst)
saveRDS(out.lst,file = sv_fl_nm)
toc()
Sys.time()



Sys.time()
tic("ev_cl_hugeD")
base_file_name <- "ev_cl_hugeD4"
tmp_fnm_base <- str_c("tmp_",base_file_name,sep="")
lst_nm <- str_c(base_file_name,"lst",sep=".")
sv_fl_nm <-str_c(lst_nm,"rds",sep=".")
out.lst <- vector(mode="list",length=B)

for (i in 1:B){
  d2000 <- sim(x = m10,n = 2000,p = g.m10.par,X=fx2000)
  out.lst[[i]] <- est_M_clbrtd(data0 = d2000,fixdcov = fx2000,mdnms = m13nms,
                               cmp_lst = cmp_lst,
                               evType = "BIC",top_B =cB,bttm_B =cB ,
                               cnsrvtv=F, sv_ev_bt.mat = F,sv_Bidx = F,
                               ev_cnf_bnds = c(0.20,0.5,0.5,0.8),
                               ev_cnf_bnd_type=c("L","L","U","U"),
                               Rl_cnf_bnds = c(0.02, .5, 0.98),
                               prlll = T,talkback = F,
                               note = "version 4 final run")
  # tmp_fnm <- str_c(tmp_fnm_base,i,"rds",sep=".")
  # saveRDS(out.lst[[i]], file=tmp_fnm)
}
assign(x = lst_nm,value = out.lst)
saveRDS(out.lst,file = sv_fl_nm)
toc()
Sys.time()


#end of singleFile.R


##### Packages we need ##### 
# library(foreign)
# library(rtiff) # voor de plotjes
# library(lavaan)
# library(psych)
# library(Hmisc)


library(foreign);library(gdata);library(RGtk2);library(ppcor);
library(nlme);library(rtiff);library(psych);library(lavaan);library(MASS)
library(shiny)
campbell <- read.spss("C:/Users/admin/Dropbox/APIM/campbell.sav", to.data.frame = T)

##### main function #####
#' Main Function
#' 
#' Description of main function
#' 
#' @param xname TODO: Describe arguments
#' @return Object of class apim.
#' @examples
#' res <- apim_sem(dname = "gender", xname = "neurotocism", yname="distress", 
#'  namex1 = "mneuro", namex2 = "fneuro", namey1 = "mdistr", 
#'  namey2 = "fdistr",labs1 = "husband", labs2 = "wive", 
#'  lab1 = "husbands", lab2 = "wives",data = campbell)
#'  
#' res
#' @export
#' @import lavaan

apim_sem <- function(xname="x", yname="y", namex1, namex2, namey1, namey2,
                        lab1="roles 1", lab2="roles 2", labs1 = "role 1", labs2="rol2", 
                        dname, data, relx = .8, rely = .8, alpha = .05, dvarr=TRUE,  
                        fix1k = TRUE, boottrials=50L, fix2k = TRUE){
 
  md <- createInput(dname = dname, xname = xname, yname=yname, 
                    namex1 = namex1, namex2 = namex2, namey1 = namey1, 
                    namey2 = namey2, labs1 = labs1, labs2 = labs2, 
                    lab1 = lab1, lab2 = lab2, data = campbell)
  
  md <- descdyad(md)
  
  md <- fitBasicAPIM(md)
  
  # md <- k_refit(md,loww=-1.5, upp=0.2, i=1, label=c("husbands","wives"))
  
  md <- FitRestrictedModel(md)
  
  md <- reliab(md)
  
  md <- fitCorrModel(md)
  
  class(md) <- "apim"
  
  return(md)
  
}  

##### Create Input #####

createInput <- function(xname="x", yname="y", namex1, namex2, namey1, namey2,
                        lab1="roles 1", lab2="roles 2", labs1 = "role 1", labs2="rol2", 
                        dname, data, relx = .8, rely = .8, alpha = .05, dvarr=TRUE,  
                        fix1k = TRUE, boottrials=50L, fix2k = TRUE){
  # the dataset:
#     ext <- substr(ifilename, nchar(ifilename)-3+1, nchar(ifilename))
#     if (ext=='sav') MaDa <-  read.spss (ifilename,use.value.labels=FALSE,max.value.labels=Inf,to.data.frame=TRUE)
#     if (ext=='csv') MaDa <- read.csv(file=ifilename,head=TRUE,sep=",")  
    
    MaDa <- data
  # variable names:
    vnames <- list(x=xname, y=yname)
  # change the variable names in the dataset
    names(MaDa)[names(MaDa)==namex1] <- "xv1"
    names(MaDa)[names(MaDa)==namex2] <- "xv2"
    names(MaDa)[names(MaDa)==namey1] <- "yv1"
    names(MaDa)[names(MaDa)==namey2] <- "yv2"
  # desired reliability:
   rel <- list(x=relx, y=rely)
  # labels:
    lab <- list(lab1=lab1, lab2=lab2, labs1=labs1, labs2=labs2)
  # Save everything in a list
    md <- list(varnames=vnames,
             MaDa=MaDa,
             reliab=rel,
             alpha=alpha,
             dvarr= dvarr,
             labels = lab,
             ddname=dname,
             bootstr = boottrials, 
             xv1=namex1,xv2=namex2,yv1=namey1,yv2=namey2,
             # does the user wants to fix k1 and k2 in a new model:
             fixk = fix2k, 
             singleK = fix1k
             )

  return(md)
}







##### Calculate descriptives #####
descdyad  = function (md){
  
  md$totn = nrow(md$MaDa) 
  md$x1min=min(md$MaDa$xv1, na.rm = TRUE) 
  md$x2min=min(md$MaDa$xv2, na.rm = TRUE) 
  md$x1max=max(md$MaDa$xv1, na.rm = TRUE) 
  md$x2max=max(md$MaDa$xv2, na.rm = TRUE) 
  md$x1mean=mean(md$MaDa$xv1, na.rm = TRUE) 
  md$x1sd=sd(md$MaDa$xv1, na.rm = TRUE) 
  md$x1n=describe(md$MaDa$xv1, na.rm = TRUE) [2]
  md$x2min=min(md$MaDa$xv2, na.rm = TRUE) 
  md$x2max=max(md$MaDa$xv2, na.rm = TRUE) 
  md$x2mean=mean(md$MaDa$xv2, na.rm = TRUE) 
  md$x2sd=sd(md$MaDa$xv2, na.rm = TRUE) 
  md$x2n=describe(md$MaDa$xv2, na.rm = TRUE) [2]
  md$y1min=min(md$MaDa$yv1, na.rm = TRUE) 
  md$y1max=max(md$MaDa$yv1, na.rm = TRUE) 
  md$y1mean=mean(md$MaDa$yv1, na.rm = TRUE) 
  md$y1sd=sd(md$MaDa$yv1, na.rm = TRUE) 
  md$y1n=describe(md$MaDa$yv1, na.rm = TRUE) [2]
  md$y2min=min(md$MaDa$yv2, na.rm = TRUE) 
  md$y2max=max(md$MaDa$yv2, na.rm = TRUE) 
  md$y2mean=mean(md$MaDa$yv2, na.rm = TRUE) 
  md$y2sd=sd(md$MaDa$yv2, na.rm = TRUE) 
  md$y2n=describe(md$MaDa$yv2, na.rm = TRUE) [2]
  md$MaDa$sumx = md$MaDa$xv1 + md$MaDa$xv2
  md$MaDa$diffx = (md$MaDa$xv2 - md$MaDa$xv1)^2
  md$msbx=(sd(md$sumx, na.rm = TRUE) )^2
  md$mswx=mean(md$MaDa$diffx, na.rm = TRUE) 
  md$iccx=(md$msbx-md$mswx)/(md$msbx+ md$mswx)
  md$xxn=as.numeric(describe(md$MaDa$sumx, na.rm = TRUE)[2])
  md$Fx=md$msbx/md$mswx
  md$pFx=1-pf(md$Fx,md$xxn-1,md$xxn)
  md$MaDa$sumy = md$MaDa$yv1 + md$MaDa$yv2
  md$MaDa$diffy = (md$MaDa$yv2 - md$MaDa$yv1)^2
  md$z2T = qnorm(1-md$alpha/2)
  md$CIv=(1- md$alpha)*100
  if (md$alpha< .01) md$CIv=rround((1-md$alpha)*100,1)
  if (md$alpha< .001) md$CIv=rround((1-md$alpha)*100,2)  
  md$msby=(sd(md$MaDa$sumy, na.rm = TRUE) )^2
  md$mswy=mean(md$MaDa$diffy, na.rm = TRUE) 
  md$iccy=(md$msby-md$mswy)/(md$msby+ md$mswy)
  md$yyn=as.numeric(describe(md$MaDa$sumy, na.rm = TRUE) [2])
  md$Fy=md$msby/md$mswy
  md$pFy=1-pf(md$Fy,md$yyn-1,md$yyn)
  md$di=0
  if (nrow(table(md$MaDa$xv1) )==2 & nrow(table(md$MaDa$xv2) )==2) mi$di=1
  md$p1=md$lab$labs1
  md$p2=md$lab$labs2
  md$missdat= 0
  md$mtxt = ""
  if (md$totn > md$x1n | md$totn > md$x2n |md$totn > md$y1n |md$totn > md$y2n)
  { md$missdat= 1
    md$mtxt = paste0("The number of cases for Person 1 on ",md$varnames$x," is ",md$x1n,", for Person 2 on ",md$varnames$x," is ",md$x2n,", for Person 1 on ",md$varnames$y," is ",md$y1n,", and for Person 2 on ",md$varnames$y," is ",md$y2n,".") 
  }

  return(md)
}





##### Fit Basic APIM Model #####
fitBasicAPIM <- function(md, treatmissing = "ML", bootstr = md$bootstr){
  
  datafile <- md$MaDa
  
  ## lavaan run: 1st output
  APIM1_out <- '
  yv1  ~ a1*xv1
  yv2  ~ a2*xv2
  yv1  ~ p12*xv2
  yv2  ~ p21*xv1
  xv1 ~ mx1*1
  xv2 ~ mx2*1
  yv1 ~ my1*1
  yv2 ~ my2*1
  xv1 ~~ vx1*xv1
  xv2 ~~ vx2*xv2
  yv1 ~~ vy1*yv1
  yv2 ~~ vy2*yv2
  xv2 ~~ cx*xv1
  yv2 ~~ cy*yv1
    a_diff := a1 - a2
    p_diff := p12 - p21
    k1 := p12/a1
    k2 := p21/a2
    k_diff := k1 - k2
    i_diff := my1 - my2
    a_ave := (a1 + a2)/2
    p_ave := (p12 + p21)/2
    sum1 := (p12 + a1)/2
    sum2 := (p21 + a2)/2
    isum := (my1 +my2)/2
    cont1 := a1 - p12
    cont2 := a2 - p21
      '
  
  # If user wanted to bootstrap: do
  if (md$bootstr != 0){
      cat("Start bootstrapping...","Can take awhile...", sep="\n")
      md$fit1_o <- sem(APIM1_out,fixed.x=FALSE,data=datafile,missing=treatmissing, se="boot", bootstrap= bootstr)
      Est1 <- parameterEstimates(md$fit1_o, standardized = TRUE,level=1-md$alpha, boot.ci.type="perc")
      md$Est1 <- as.data.frame(Est1)
      rownames(md$Est1) <- md$Est1[,4]
    # save the results of the bootstrap CI's in md
      md$kl1<- md$Est1["k1","ci.lower"]
      md$ku1<- md$Est1["k1","ci.upper"]
      md$kl2<- md$Est1["k2","ci.lower"]
      md$ku2<- md$Est1["k2","ci.upper"]
     # if sign difference based on bootstrap CI: save 1, if not, save 0:
      md$kdif.l <- md$Est1["k_diff","ci.lower"]
      md$kdif.u <- md$Est1["k_diff","ci.upper"]
      if (md$kdif.l < 0 & 0 < md$kdif.u){
        md$kdif.yn <- 0
      } else{md$kdif.yn <- 1}
  } else { 
    md$fit1_o <- sem(APIM1_out,fixed.x=FALSE,data=md$MaDa,missing="ML")
    md$Est1 <- parameterEstimates(md$fit1_o, standardized = TRUE)
    md$Est1 <- as.data.frame(md$Est1)
    rownames(md$Est1) <- md$Est1[,4]      
  }
 
# only proceed if the model converged:
  if (md$fit1_o@Fit@iterations > 0L && !md$fit1_o@Fit@converged){
    stop("Model did not converge")
  } else{
    # Save results of APIM1_out in md
    md$k1<- md$Est1["k1","est"]
    md$k2<- md$Est1["k2","est"]
    ppv=fitted(md$fit1_o)$cov
    md$yvar1=ppv["yv1","yv1"]
    md$xvar1=ppv["xv1","xv1"]
    md$yvar2=ppv["yv2", "yv2"]
    md$xvar2=ppv["xv2","xv2"]
    md$ryy=ppv["yv2","yv1"]/ppv["yv1","yv1"]
    md$rxx =ppv["xv2","xv1"]/ppv["xv1","xv1"]
    md$rx1y1=ppv["xv1","yv1"]/sqrt(ppv["yv1","yv1"]*ppv["xv1","xv1"])
    md$rx2y2=ppv["xv2","yv2"]/sqrt(ppv["xv2","xv2"]*ppv["yv2","yv2"])
    md$rx1y2=ppv["xv2","yv1"]/sqrt(ppv["yv1","yv1"]*ppv["xv2","xv2"])
    md$rx2y1=ppv["xv1","yv2"]/sqrt(ppv["yv2","yv2"]*ppv["xv1","xv1"])
    
    # save the means for both dist. members
    ppm=fitted(md$fit1_o)$mean
    md$ymean1=ppm["yv1"]
    md$xmean1=ppm["xv1"]
    md$ymean2=ppm["yv2"]
    md$xmean2=ppm["xv2"]
    
    # save for both dist. members
    md$intr_dp = md$Est1["i_diff","pvalue"]
    md$intr_d.CIl = md$Est1["i_diff","ci.lower"]
    md$intr_d.CIu = md$Est1["i_diff","ci.upper"]
    md$act_dp = md$Est1["a_diff","pvalue"]
    md$act_d.CIl = md$Est1["a_diff","ci.lower"]
    md$act_d.CIu = md$Est1["a_diff","ci.upper"]
    md$part_dp = md$Est1["p_diff","pvalue"]
    md$part_d.CIl = md$Est1["p_diff","ci.lower"]
    md$part_d.CIu = md$Est1["p_diff","ci.upper"]
    md$intr1 = md$Est1["my1","est"]
    md$intr1s = md$Est1["my1","std.all"]
    md$intr1z = md$Est1["my1","z"]
    md$intr1p = md$Est1["my1","pvalue"]
    md$intr1se = md$Est1["my1","se"]
    md$intr1.CIu=md$intr1+md$intr1se*md$z2T
    md$intr1.CIl=md$intr1-md$intr1se*md$z2T
    md$intr2 = md$Est1["my2","est"]
    md$intr2s = md$Est1["my2","std.all"]
    md$intr2z = md$Est1["my2","z"]
    md$intr2p = md$Est1["my2","pvalue"]
    md$intr2se = md$Est1["my2","se"]
    md$intr2.CIu=md$intr2+md$intr2se*md$z2T
    md$intr2.CIl=md$intr2-md$intr2se*md$z2T
    md$intro = md$Est1["i_diff","est"]
    md$introp = md$Est1["i_diff","pvalue"]
    md$introse = md$Est1["i_diff","se"]
    md$intro.CIu=md$intro+md$introse*md$z2T
    md$intro.CIl=md$intro-md$introse*md$z2T
    md$acto = md$Est1["a_ave","est"]
    md$actop = md$Est1["a_ave","pvalue"]
    md$actose = md$Est1["a_ave","se"]
    md$acto.CIu=md$acto+md$actose*md$z2T
    md$acto.CIl=md$acto-md$actose*md$z2T
    md$parto = md$Est1["p_ave","est"]
    md$partop = md$Est1["p_ave","pvalue"]
    md$partose = md$Est1["p_ave","se"]
    md$parto.CIu=md$parto+md$partose*md$z2T
    md$parto.CIl=md$parto-md$partose*md$z2T
  
    # saved the results for both actor (partner) effects
    md$act1 = md$Est1["a1","est"]
    md$act1ss = md$Est1["a1","std.all"]
    md$act1z = md$Est1["a1","z"]
    md$act1p = md$Est1["a1","pvalue"]
    md$act1se = md$Est1["a1","se"]
    md$act1.CIu=md$act1+md$act1se*md$z2T
    md$act1.CIl=md$act1-md$act1se*md$z2T
    md$act2 = md$Est1["a2","est"]
    md$act2ss = md$Est1["a2","std.all"]
    md$act2z = md$Est1["a2","z"]
    md$act2p = md$Est1["a2","pvalue"]
    md$act2se = md$Est1["a2","se"]
    md$act2.CIu=md$act2+md$act2se*md$z2T
    md$act2.CIl=md$act2-md$act2se*md$z2T
    md$part1 = md$Est1["p12","est"]
    md$part1ss = md$Est1["p12","std.all"]
    md$part1z = md$Est1["p12","z"]
    md$part1p = md$Est1["p12","pvalue"]
    md$part1se = md$Est1["p12","se"]
    md$part1.CIu=md$part1+md$part1se*md$z2T
    md$part1.CIl=md$part1-md$part1se*md$z2T    
    md$part2 = md$Est1["p21","est"]
    md$part2ss = md$Est1["p21","std.all"]
    md$part2z = md$Est1["p21","z"]
    md$part2p = md$Est1["p21","pvalue"]
    md$part2se = md$Est1["p21","se"]
    md$part2.CIu=md$part2+md$part2se*md$z2T
    md$part2.CIl=md$part2-md$part2se*md$z2T
    md$iccY=md$Est1["cy","est"]
    md$iccYs=md$Est1["cy","std.all"]
    md$iccYz=md$Est1["cy","z"]
    md$iccYp=md$Est1["cy","pvalue"]
    md$iccY.CIl=md$Est1["cy","ci.lower"]
    md$iccY.CIu=md$Est1["cy","ci.upper"]
    md$iccX=md$Est1 ["cx","est"]
    md$iccXs=md$Est1["cx","std.all"]
    md$iccXz=md$Est1["cx","z"]
    md$iccXp=md$Est1["cx","pvalue"]
    md$iccX.CIl=md$Est1["cx","ci.lower"]
    md$iccX.CIu=md$Est1["cx","ci.upper"]
    
    # saved for both dist. roles
    md$varx1=md$Est1["vx1","est"]
    md$varx2=md$Est1["vx2","est"]
    md$esd1=md$Est1["vy1","est"]
    md$esd2=md$Est1["vy2","est"]
    md$esd1p=onetailp(md$Est1["vy1","pvalue"])
    md$esd2p=onetailp(md$Est1["vy2","pvalue"])
    md$varxm1p=onetailp(md$Est1["vx1","pvalue"])
    md$varxm2p=onetailp(md$Est1["vx2","pvalue"])
    md$esd1.CIl=md$Est1["vy1","ci.lower"]
    md$esd1.CIu=md$Est1["vy1","ci.upper"]
    md$esd2.CIl=md$Est1["vy2","ci.lower"]
    md$esd2.CIu=md$Est1["vy2","ci.upper"]
    md$varxm1.CIl=md$Est1["vx1","ci.lower"]
    md$varxm1.CIu=md$Est1["vx1","ci.upper"]
    md$varxm2.CIl=md$Est1["vx2","ci.lower"]
    md$varxm2.CIl=md$Est1["vx2","ci.upper"]
    md$mmeanx1=md$Est1["mx1","est"]
    md$mmeanx1p=md$Est1["mx1","pvalue"]
    md$mmeanx1.CIl=md$Est1["mx1","ci.lower"]
    md$mmeanx1.CIu=md$Est1["mx1","ci.upper"]
    md$mmeany1=md$Est1["my1","est"]
    md$mmeany1p=md$Est1["my1","pvalue"]
    md$mmeany1.CIl=md$Est1["my1","ci.lower"]
    md$mmeany1.CIu=md$Est1["my1","ci.upper"]
    md$mmeanx2=md$Est1["mx2","est"]
    md$mmeanx2p=md$Est1["mx2","pvalue"]
    md$mmeanx2.CIl=md$Est1["mx2","ci.lower"]
    md$mmeanx2.CIu=md$Est1["mx2","ci.upper"]
    md$mmeany2=md$Est1["my2","est"]
    md$mmeany2p=md$Est1["my2","pvalue"]
    md$mmeany2.CIl=md$Est1["my2","ci.lower"]
    md$mmeany2.CIu=md$Est1["my2","ci.upper"]
    
    # calculate r-squared for both dist. members.
    md$rsq1a= inspect(md$fit1_o, 'r2')[1]
    md$rsq1b= inspect(md$fit1_o, 'r2')[2]
    
    # fit indices
    md$cs1=fitMeasures(md$fit1_o)[2]
    md$df1=fitMeasures(md$fit1_o)[3]
    md$pv1=fitMeasures(md$fit1_o)[4]
    md$tli1=fitMeasures(md$fit1_o)[9]
    md$rmsea1=fitMeasures(md$fit1_o)[23]
    md$rmseal1=fitMeasures(md$fit1_o)[24]
    md$rmseau1=fitMeasures(md$fit1_o)[25]
    }

    return(md)
  }






##### k_refit_function #####
k_refit <- function (md,loww,upp, i, label) {
  
  if ((loww <  -1.0 & upp > 1.0) |  (loww <  0.0 & upp > 1.0)  | (loww <  -1.0 & upp > 0.0)) {
    md$no_k[i] <- paste("The confidence interval of the", label[i],
               "is too wide for determining which single model is most likely. Therefore, the value of k is not fixed to a certain value in a new model to determine its fit.")
    md$k[i] <- NA
  } else if (loww <  1 & upp > 1) { 
    md$k[i] <- 1
  } else if (loww <  0 & upp > 0){
    md$k[i] <- 0
  } else if (loww <  -1 & upp > -1) { 
    md$k[i] <- -1
  } else if (loww > 0 & upp < 1) {
    md$no_k[i] <- paste("It can be concluded that for the ",label[i],
               "the model is in between the actor-only (k = 0) and the couple (k = 1) models. Therefore k cannot be fixed to a specific value in a new model in order to determine its fit.")
    md$k[i] <- NA
  } else if (loww > -1 & upp < 0) {
    md$no_k[i] <- c("It can be concluded that for the ",label[i],
               "the model is in between the contrast (k = -1) and the actor-only (k = 0) models. Therefore k cannot be fixed to a specific value in a new model in order to determine its fit.")
    md$k[i] <- NA
  } else if (upp < -1) {
    md$no_k[i] <- c("It can be concluded that for the ",label[i],
               "the model is more extreme than the contrast model (k = -1). Therefore k cannot be fixed to a specific value in a new model in order to determine its fit.")
    md$k[i] <- NA
  } else if (loww > 1) {
    md$no_k[i] <- c("It can be concluded that for the ",label[i],
               "the model is more extreme than the couple model (k = 1). Therefore k cannot be fixed to a specific value in a new model in order to determine its fit.")
    md$k[i] <- NA
  }
  
  return(md)
}






##### Fit Restricted model #####
FitRestrictedModel <- function(md, fixk = TRUE, treatmissing="ML", bootstrapSingleK = 50L){
  
  singleK = md$singleK
  datafile = md$MaDa
  
   # if user wants to fix k1 and k2
   if (fixk == TRUE){
     # if actor values have substantial values
      if (md$act1p > md$alpha & md$act2p > md$alpha & abs(md$act1ss) < .1 & abs(md$act2ss) < .1){
        stop("The standardized actor effect(s) do not have a substantial value. Estimates of k would be highly unstable. Therefore, all further analyses are aborted.") 
      } else {          
        
        # Run the k_refit-function
        label <- c(md$lab$lab1, md$lab$lab2)
        k <- list()
        no_k <- list()
        for (i in 1:2){
          eval(parse(text=paste("md <- k_refit", "(md, md$kl", i, ", md$ku", 
                                i,", i =", i,",label=label", ")", sep="")))
        }       
        #     rm(k, no_k, label)
        #     md$k <- NULL
        #     md$no_k <- NULL
        #     md$k[2] <- md$k[1]
        
        # Evaluate if both k1 and k2 could be fixed, if so, continue:
        if (!is.na(md$k[1]) & !is.na(md$k[2])){
          md$k1 <- md$k[1]
          md$k2 <- md$k[2]
          # Run the lavaan model with the k's constrained to the previously defined values
          APIM2_2k <- paste0('
           yv1  ~ a1*xv1
           yv2  ~ a2*xv2
           P1 = ~ 0
           P2 = ~ 0
           P1  ~ a1*xv2
           P2  ~ a2*xv1
           yv1  ~ ',md$k[1],'*P1
           yv2  ~ ',md$k[2],'*P2
           xv1 ~ mx1*1
           xv2 ~ mx2*1
           yv1 ~ my1*1
           yv2 ~ my2*1
           xv1 ~~ vx1*xv1
           xv2 ~~ vx2*xv2
           yv1 ~~ vy1*yv1
           yv2 ~~ vy2*yv2
           xv2 ~~ c1*xv1
           yv2 ~~ c2*yv1'
           )
          fit2_2k <- sem(APIM2_2k,fixed.x=FALSE, data=datafile,missing=treatmissing)
          Est2_2k <- parameterEstimates(fit2_2k, standardized = TRUE,level=1-md$alpha)
          
          # Subtract Chi-square test statistic to come to a conclusion
          md$csk2 <- fitMeasures(fit2_2k)["chisq"]
          md$dfk2 <- fitMeasures(fit2_2k)["df"]
          md$pvk2 <- fitMeasures(fit2_2k)["pvalue"]
          
            
        # If one or both values of k cannot be derived from their CI AND no sign difference between k1 and k2
        # AND user specified he/she wants a single k in their model
      } else if (md$singleK == FALSE){
          # if the user does not want to calculate a single k: abort the analyses (since the k's cannot be fixed)
          stop("The analyses is aborted since one or both k's could not be constrained to a single value.")
          
      } else if ((is.na(md$k[1]) | is.na(md$k[2]))  &&  md$kdif.yn == 0 && md$singleK == TRUE){      
          # if the user wants to: specify a model with a single k and bootstrap the CI
          APIM2_1k <- '
          yv1  ~ a1*xv1
          yv2  ~ a2*xv2
          P1 = ~ 0
          P2 = ~ 0
          P1  ~ a1*xv2
          P2  ~ a2*xv1
          yv1  ~ k*P1
          yv2  ~ k*P2
          xv1 ~ mx1*1
          xv2 ~ mx2*1
          yv1 ~ my1*1
          yv2 ~ my2*1
          xv1 ~~ vx1*xv1
          xv2 ~~ vx2*xv2
          yv1 ~~ vy1*yv1
          yv2 ~~ vy2*yv2
          xv2 ~~ c1*xv1
          yv2 ~~ c2*yv1'
          
          cat("Model with a single k is fitted...","Start bootstrapping...","Can take awhile...", sep="\n")
          fit2_kb <- sem(APIM2_1k,fixed.x=FALSE,data=md$MaDa,missing=treatmissing, se = "boot", bootstrap= bootstrapSingleK)
          Est2_kb <- parameterEstimates(fit2_kb, standardized = TRUE,level=1-md$alpha, boot.ci.type="perc")
          
          # save the results in md
          md$kl <- Est2_kb[7,9]
          md$ku <- Est2_kb[7,10]
          md$k <- Est2_kb[7,5] 
        
          # To what value can this k be fixed based on its CI?
          label <- paste(md$lab$lab1,"and", md$lab$lab2, sep=" ")      
          k <- list()
          no_k <- list()
          md <- k_refit(md, md$kl, md$ku, i = 1, label=label)   
          if (!is.na(md$k[1])){
              md$k <- md$k[1]
              APIM3_1k <- paste0('
              yv1  ~ a1*xv1
              yv2  ~ a2*xv2
              P1 = ~ 0
              P2 = ~ 0
              P1  ~ a1*xv2
              P2  ~ a2*xv1
              yv1  ~ ',md$k[1],'*P1
              yv2  ~ ',md$k[1],'*P2
              xv1 ~ mx1*1
              xv2 ~ mx2*1
              yv1 ~ my1*1
              yv2 ~ my2*1
              xv1 ~~ vx1*xv1
              xv2 ~~ vx2*xv2
              yv1 ~~ vy1*yv1
              yv2 ~~ vy2*yv2
              xv2 ~~ c1*xv1
              yv2 ~~ c2*yv1'
              )
          
            fit3_1k <- sem(APIM3_1k,fixed.x=FALSE, data=datafile,missing=treatmissing)
            Est3_1k <- parameterEstimates(fit3_1k, standardized = TRUE,level=1-md$alpha)
            
            # Subtract Chi-square test statistic to come to a conclusion
            md$csk <- fitMeasures(fit3_1k)["chisq"]
            md$dfk <- fitMeasures(fit3_1k)["df"]
            md$pvk <- fitMeasures(fit3_1k)["pvalue"]          
          }}}}
  
 return(md)
}
      


##### Compute Reliability #####
reliab = function (md){ 
  md$relok =1-md$reliab$y*md$reliab$x
  md$reloor=0
  md$texrel=""
  #DAK With Distinguishable dyads, in principle the reliablity for the two members on X could be different; I guess for now we just presume they are the same and say something about that in the output
  md$texrel=paste0("The estimates here presume that the reliability of measurement of ",md$varnames$x," is ",rrround(md$reliab$x,2)," and of ",md$varnames$y," is ",rrround(md$reliab$y,2),".")
  if (md$reliab$x < abs(md$rxx))
  { md$texrel=paste0(md$texrel,"  The proposed reliablity for ",md$varnames$x," must be at least as large as the correlation between the two persons' scores on ",md$varnames$x," which is not the case.")
    md$reloor = 1  } 
  if (md$reliab$y < abs(md$ryy))
  { md$texrel=paste0(md$texrel,"  The proposed reliablity for ",md$varnames$y," must be at least as large as the correlation between the two persons' scores on ",md$varnames$y," which is not the case.")
    md$reloor = 1 } 
  md$reloor = 0
  if (sqrt(md$reliab$x*md$reliab$y) < abs(md$rx1y1) | sqrt(md$reliab$x*md$reliab$y) < abs(md$rx2y2) |  sqrt(md$reliab$x*md$reliab$y) < abs(md$rx1y2)  | sqrt(md$reliab$x*md$reliab$y) < abs(md$rx2y1))
  { md$texrel=paste0(md$texrel,"  The proposed reliablities are too low given the magnitude of the correlations between ",md$varnames$x," and ",md$varnames$y,".")
    md$reloor = 1} 
  if (md$reloor == 1)md$texrel=paste0(md$texrel,"  Correction for measurement error cannot be preformed.  Please pick new and higher values of the reliabilites or correlated errors if you want corrections for measurement error.") 
  # I made these for both dist. members
  md$xesd1=sqrt(md$xvar1*(1-md$reliab$x))
  md$xesd2=sqrt(md$xvar2*(1-md$reliab$x))
  md$yesd1=sqrt(md$yvar1*(1-md$reliab$y))
  md$yesd2=sqrt(md$yvar2*(1-md$reliab$y))
  md <<- as.list(md)
  return(md)
}





##### Run correcting for unreliablity #####
fitCorrModel <- function(md, treatmissing = "ML"){
  
  datafile = md$MaDa
  bootstr = md$bootstr
  
  if (md$relok < 0 | md$reloor != 0){
    stop("The model correcting for unreliability was not fitted")
  } else{
    # lavaan run:
    APIM3 <- paste0('
    xv1 ~~ 0*xv1  
    xv2 ~~ 0*xv2
    yv1 ~~ 0*yv1
    yv2 ~~ 0*yv2
    xx1 ~~ 1*xx1  
    xx2 ~~ 1*xx2
    yy1 ~~ 1*yy1
    yy2 ~~ 1*yy2
    xx1 =~ 0
    xx2 =~ 0
    yy1 =~ 0
    yy2 =~ 0
    xv1  ~ ',md$xesd1,'*xx1
    yv1  ~ ',md$yesd1,'*yy1
    xv2  ~ ',md$xesd2,'*xx2
    yv2  ~ ',md$yesd2,'*yy2
    xfact1  =~ 1*xv1
    yfact1  =~ 1*yv1
    xfact2  =~ 1*xv2
    yfact2  =~ 1*yv2
    yfact1  ~ a1*xfact1 + p12*xfact2 
    yfact2  ~ a2*xfact2 + p21*xfact1
    xv1 ~ mx1*1
    xv2 ~ mx2*1
    yv1 ~ my1*1
    yv2 ~ my2*1
    xfact1  ~~ vx1*xfact1
    xfact2  ~~ vx2*xfact2
    yfact1  ~~ vy1*yfact1
    yfact2  ~~ vy2*yfact2
    xfact1  ~~ xfact2
    yfact1  ~~ yfact2
    ck1 := p12/a1
    ck2 := p21/a2
    ck_diff := ck1 - ck2
    ')
    cat("Start bootstrapping...","Can take awhile...", sep="\n")
    md$fit3 <- sem(APIM3,fixed.x=FALSE,data=datafile,missing=treatmissing, se="boot", bootstrap = md$bootstr)
    Est3 = parameterEstimates(md$fit3, standardized = TRUE,level=1-md$alpha, boot.ci.type="perc")
    
    # Only proceed if the model converged:
    if (md$fit3@Fit@iterations > 0L && !md$fit3@Fit@converged){
      stop("Model did not converge")
    } else{
      
      Est3_sub <- as.data.frame(Est3[c(21:34, 43:45),])
      Est3_sub[c(13:14),4] <- c("icx", "icy")
      rownames(Est3_sub) <- Est3_sub[,4]  
      
      
      # Save estimates  
      # k
      md$ck1 <- Est3_sub["ck1","est"]
      md$ck2 <- Est3_sub["ck2","est"]
      md$ck1.CIl <- Est3_sub["ck1","ci.lower"]
      md$ck2.CIl <- Est3_sub["ck2","ci.lower"]
      md$ck1.CIu <- Est3_sub["ck1","ci.upper"]
      md$ck2.CIu <- Est3_sub["ck2","ci.upper"]
      md$ckdif.CIl <- Est3_sub["ck_diff","ci.lower"]
      md$ckdif.CIu <- Est3_sub["ck_diff","ci.upper"]
      if (md$ckdif.CIl < 0 & 0 < md$ckdif.CIu){
        md$ckdif.yn <- 1
      } else{md$ckdif.yn <- 0}
      
      # actor effects for both dist. members
      md$act13 <- Est3_sub["a1","est"]
      md$act1s3 <- Est3_sub["a1","std.all"]
      md$act1z3 <- Est3_sub["a1","z"]
      md$act1p3 <- Est3_sub["a1","pvalue"]
      md$act1se3 <- Est3_sub["a1","se"]
      md$act1.CIu3<-md$act13+md$act1se3*md$z2T
      md$act1.CIl3<-md$act13-md$act1se3*md$z2T  
      md$act23 <- Est3_sub["a2","est"]
      md$act2s3 <- Est3_sub["a2","std.all"]
      md$act2z3 <- Est3_sub["a2","z"]
      md$act2p3 <- Est3_sub["a2","pvalue"]
      md$act2se3 <- Est3_sub["a2","se"]
      md$act2.CIu3<-md$act23+md$act2se3*md$z2T
      md$act2.CIl3<-md$act23-md$act2se3*md$z2T
      
      # partner effects for both dist. members
      md$part13 <- Est3_sub["p12","est"]
      md$parts13 <- Est3_sub["p12","std.all"]
      md$partz13 <- Est3_sub["p12","z"]
      md$partp13 <- Est3_sub["p12","pvalue"]
      md$part13se <- Est3_sub["p12","se"]
      md$part1.CIu3<-md$part13+md$part13se*md$z2T
      md$part1.CIl3<-md$part13-md$part13se*md$z2T
      md$part23 <- Est3_sub["p21","est"]
      md$parts23 <- Est3_sub["p21","std.all"]
      md$partz23 <- Est3_sub["p21","z"]
      md$partp23 <- Est3_sub["p21","pvalue"]
      md$part23se <- Est3_sub["p21","se"]
      md$part1.CIu3<-md$part23+md$part23se*md$z2T
      md$part1.CIl3<-md$part23-md$part23se*md$z2T
      
      md$iccY3<- Est3_sub["icy","est"]
      md$iccYs3<-Est3_sub["icy","std.all"]
      md$iccYz3<-Est3_sub["icy","z"]
      md$iccYp3<-Est3_sub["icy","pvalue"]
      md$iccY.CIl3<-Est3_sub["icy","ci.lower"]
      md$iccY.CIu3<-Est3_sub["icy","ci.upper"]
      
      md$iccX3<-Est3_sub ["icy","est"]
      md$iccXs3<-Est3_sub["icy","std.all"]
      md$iccXz3<-Est3_sub["icy","z"]
      md$iccXp3<-Est3_sub["icy","pvalue"]
      md$iccX.CIl3<-Est3_sub["icy","ci.lower"]
      md$iccX.CIu3<-Est3_sub["icy","ci.upper"]
      
      # esd for both dist. members
      md$esd13<-Est3_sub["vy1","est"]
      md$esd23<-Est3_sub["vy2","est"]
      md$esd1p3<-onetailp(Est3_sub["vy1","pvalue"])
      md$esd2p3<-onetailp(Est3_sub["vy2","pvalue"])
      md$esd1.CIl3<-Est3_sub["vy1","ci.lower"]
      md$esd1.CIu3<-Est3_sub["vy1","ci.upper"]
      md$esd2.CIl3<-Est3_sub["vy2","ci.lower"]
      md$esd2.CIu3<-Est3_sub["vy2","ci.upper"]
      
      
      # r-squared for both dist. members
      md$rsq3a<- 0
      md$rsq3b<- 0
      if (md$esd13 > 0) md$rsq3a= inspect(md$fit3, 'r2')[9]
      if (md$esd23 > 0) md$rsq3b= inspect(md$fit3, 'r2')[10]
    }
  }
  
  return(md)
  }






##### Text output #####

print.apim <- function(md){
  
out<-capture.output(summary(md$fit1_o))

# write output in a textfile.
cat("\n                                 APIM: Distinguishable Dyads Using lavaan \n ",
    "####################################################################",
    "#              Standard Model (Distinguishable members)            #",
    "#################################################################### \n",
    out, sep = "\n")

# if no sign. difference in actor nor partner effects: advise to use indistinguishable method
if (md$Est1[15,]$pvalue > md$alpha & md$Est1[16,]$pvalue > md$alpha){
  cat("\n ***WARNING*** \n Please use the method for indistinguishable members since a model with equal actor and partner effects is more appropriate for this dataset! \n",sep = "\n")}

cat("\n  \n",
    "#################################################################### \n",
    "#        DataToText: Standard Model (Distinguishable members)      # \n",
    "#################################################################### \n \n",
    "Introduction \n",
    "------------  \n",
    "      The focus of this study is the investigation of the effect of ", md$varnames$x," on ", md$varnames$y," for ",md$lab$lab1," and ",md$lab$lab2,".",
    " Both the effect of own ",md$varnames$x," (actor) and the effect of partner's ", md$varnames$x," (partner) on ", md$varnames$y," are studied. ",
    "The total number of dyads is ", md$totn,        
    sep="")
if (md$missdat == 0) cat(", and there are no missing data. ")
if (md$missdat == 1) cat(", and there are missing data for one or more of the variables:",md$mtxt,".")
cat("The dyad members are treated as if they are distinguishable by the variable ",md$ddname,".",sep = "")
if (.05 != md$alpha){cat("\n For all these analyses, alpha is set at",md$alpha,".",sep = "")}



cat("\n  \n",
    "Analyses \n",
    "-------- \n",
    "      The analyses use structural equation modeling with maximum likelihood estimation using the program lavaan. The tests of coefficients are Z tests. ",
    "Effect sizes for actor and partner effects are partial correlations and d when the predictor is dichotomous. Betas are given twice, one using the overall standard deviation across all persons (o) for standardization and a second using the standard deviation for ",
    md$lab$lab1," and ",md$lab$lab2," separately (s). If betas are to be compared across members, the beta (o) value should be examined. For all these analyses, alpha is set at ",md$alpha,".",
    " The descriptive statistics of the raw variables are contained in Table 1 and the FIML estimated means and standard deviations in Table 2. ", 
    sep ="")


if (md$missdat == 0){cat("Because there are no missing data, the raw and the FIML means should be the same the FIML standard deviations differ from the raw ones by the square root of N - 1 divided by N.", sep ="")}

cat("\n  \n",
    "Results \n",
    "------- \n",
    sep ="")

# Only output the text if model converged:
if (md$fit1_o@Fit@iterations > 0L && !md$fit1_o@Fit@converged){
  cat("\n       ERROR: model did not converge \n \n",sep="")
} else{
  cat("\n       The lavaan model converged after ", md$fit1_o@Fit@iterations, " iterations. \n      A summary of results of the APIM analyses is contained in Table 3 and the overall effects in Table 4.", 
      sep="")
  
  cat(" The variance of the errors for the ", md$lab$lab1, " and ", md$lab$lab2, " are ", rround(md$esd1,3), " and ",rround(md$esd2,3), ", respectively.",
      " The R squared for the ", md$lab$lab1," is ",rrround(md$rsq1a,3),", for the ", md$lab$lab2, " it is ", rrround(md$rsq1b,3),
      ". The partial intraclass correlation for ",md$varnames$y," controlling for actor and partner variables is equal to ",rrround(md$iccYs,3)," and is",knt(md$iccYp,md$alpha),
      " statistically significant (",pval(md$iccYp),", 95% CI [", md$iccY.CIl, ", ", md$iccY.CIu, "])",
      sep="")
  if (md$iccY>.1) cat(" Thus, when one member of the dyad scores high (low) on this construct, the other member also tends to have a high (low) score. \n", sep="")
  if (md$iccY < -.1) cat(" Thus, there is no clear trend in the scores of both members of the dyad. \n",sep="")
  cat("\n        The intercept (the predicted score on ",md$varnames$y," when the actor and partner variables equal zero) for ",md$lab$lab1," is equal to ",rround(md$intr1,3)," and is",knt(md$intr1p,md$alpha),
      " statistically significant (",pval(md$intr1p), ", 95% CI [", round(md$intr1.CIl,2), ", ", round(md$intr1.CIu,2), "]).",  
      " The intercept for ",md$lab$lab2," is equal to ",rround(md$intr2,3)," and is",knt(md$intr2p,md$alpha)," statistically significant (",
      pval(md$intr2p),  ", 95% CI [", round(md$intr2.CIl,2),", ", round(md$intr2.CIu,2), "]).",    
      " The test that the two intercepts are equal is",knt(md$intr_dp,md$alpha)," statistically significant (",
      pval(md$intr_dp),", 95% CI [", round(md$intr_d.CIl,2),", ", round(md$intr_d.CIu,2), "]).", 
      "  The difference in intercepts is equal to ",rround(md$intro,3)," and is",knt(md$introp,md$alpha)," statistically significant (",
      pval(md$introp),", 95% CI [", round(md$intro.CIl,2),", ", round(md$intro.CIu,2), "]).", sep="")
  
  # DAK I think this goes.
  # Only describe the output when actor or partner effects differ signifcantly
  # if (Est1[15,]$pvalue > .05 & Est1[16,]$pvalue > .05){
  #  cat("\n \n The method for indistinguishable members would be more appropriate for this dataset since no difference in actor nor partner effects could be withdrawn.")
  #}
  
  # report about actor and partner effects
  #Converts the actor effect to a partial correlation, an effect size measure
  
  # Check these
  md$act1r = md$act1ss*sqrt((1-md$rxx^2)/(1-md$rx2y1^2))
  md$act2r = md$act2ss*sqrt((1-md$rxx^2)/(1-md$rx1y2^2))
  md$part1r = md$part1ss*sqrt((1-md$ryy^2)/(1-md$rx1y2^2))
  md$part2r = md$part2ss*sqrt((1-md$ryy^2)/(1-md$rx2y1^2))
  
  
  md$xsd = sqrt((md$xvar1 + md$xvar1)/2)
  md$ysd = sqrt((md$yvar1 + md$yvar1)/2)
  
  md$act1so= md$act1*md$xsd/md$ysd
  md$act2so= md$act2*md$xsd/md$ysd
  md$part1so= md$part1*md$xsd/md$ysd
  md$part2so= md$part2*md$xsd/md$ysd
  
  esst1=md$act1r
  esst2=md$act2r
  if (md$di==1) {
    esst1=md$act1/md$esd1
    esst2=md$act2/md$esd2}
  ttu1=paste0(effszm(md$di), rrround(esst1,3),effsz(esst1,md$di),").")
  ttu2=paste0(effszm(md$di), rrround(esst2,3),effsz(esst2,md$di),").")
  
  cat("\n       The actor effect for the ", md$lab$lab1, " is equal to ",rround(md$act1,3)," (",
      pval(md$act1p),", 95% CI [", round(md$act1.CIl,2),", ", round(md$act1.CIu,2), "])",
      ". The the overall standardized effect for the ",md$lab$lab1," is ", rround(md$act1so,3),ttu1,
      "  The actor effect for the ", md$lab$lab2, " is equal to ",rround(md$act2,3)," (",
      pval(md$act2p),", 95% CI [", round(md$act2.CIl,2),", ", round(md$act2.CIu,2),"])",
      " and the overall standardized actor effect for the ",md$lab$lab2," is ", rround(md$act2so,3),ttu2, 
      "  The test that the two actor effects are equal is",knt(md$act_dp,md$alpha)," statistically significant",
      pval(md$act_dp),", 95% CI [", round(md$act_d.CIl,2), ", ", round(md$act_d.CIu,2),"]).",
      " The overall actor effect is equal to ",rround(md$acto,3)," and is",knt(md$actop,md$alpha)," statistically significant (",
      pval(md$actop),", 95% CI [", round(md$acto.CIl,2), ", ", round(md$acto.CIu,2),"]).",
      sep="")
  
  esst1=md$part1r
  esst2=md$part2r
  if (md$di==1) {
    esst1=md$part1/md$esd1
    esst2=md$part2/md$esd2}
  ttu1=paste0(effszm(md$di), rrround(esst1,  3),effsz(esst1,md$di),").")
  ttu2=paste0(effszm(md$di), rrround(esst2,  3),effsz(esst2,md$di),").")
  
  cat("\n  \n")
  
  cat("       The partner effect of the ", md$lab$lab1, " on the ",md$lab$lab2," is equal to ",rround(md$part1,3),", which is",knt(md$part1p,md$alpha)," statistically significant (",
      pval(md$part1p),", 95% CI [", round(md$part1.CIl,2),", ", round(md$part1.CIu,2), "])",
      ", and its overall standardized effect is ",rround(md$part1so,3),ttu1, 
      " The partner effect of the ", md$lab$lab2," on the ", md$lab$lab1," is equal to ",rround(md$part2,3)," and is",knt(md$part2p,md$alpha)," statistically significant (",
      pval(md$part2p),", 95% CI [", round(md$part2.CIl,2),", ", round(md$part2.CIu,2), "])",
      " and its overall standardized partner effect is ",rround(md$part2so,3),ttu2, 
      "  The test that the two partner effects are equal is",knt(md$part_dp,md$alpha)," statistically significant (",
      pval(md$part_dp),", 95% CI [", round(md$part_d.CIl,2),", ", round(md$part_d.CIu,2), "]).",
      "  The overall partner effect is equal to ",rround(md$parto,3)," is",knt(md$partop,md$alpha)," statistically significant (",
      pval(md$partop),", 95% CI [", round(md$parto.CIl,2),", ", round(md$parto.CIu,2), "]).",
      sep="")
  
  if (md$act1p > md$alpha & md$act2p > md$alpha & abs(md$act1ss) < .1 & abs(md$act2ss) < .1){
    stop("The standardized actor effect(s) do not have a substantial value. Therefore, estimates of k would be highly unstable and all further analyses are aborted.") 
  } else{
    cat("\n \n       Because the standardized actor effects of both", md$lab$lab1, "and", md$lab$lab2, 
        "are greater than .1 in absolute value and they are statistically significant, k (i.e. the ratio of the partner effect to the actor effect) can be interpreted in the output. ")
    
    # DAK I think we still report the two k's even if they do not differ.
    #    if (md$kp < .05)
    
    cat("Here, k1 represents the effect of ", md$lab$lab1, " on ", md$lab$lab2, " and k2 the effect of ", md$lab$lab2, " on ", md$lab$lab1,
        #          ".  Because k1 differs significantly from k2, both k's need to be interpreted separately 
        ". The value of k1 equals ", round(md$k1,2), " and k2 equals ", 
        
        round(md$k2,2),".  Their 95% percentile bootstrap confidence intervals range from ", round(md$kl1,2), " to ", round(md$ku1,2), " and ", round(md$kl2,2), 
        " to ", round(md$ku2,2), ", respectively. So, for the ", md$lab$lab1, ", it can be concluded that ", kcon(md$kl1,md$ku1),
        " For the ", md$lab$lab2, ", the confidence interval of k2 shows that ", kcon(md$kl2,md$ku2),
        sep=""
    )

    
    if (!is.na(md$k[1]) & !is.na(md$k[2])){
      cat(" When investigating the relative fit of the model with k1 fixed to", md$k[1], "and k2 fixed to",  md$k[2]
          )   
      if (md$pvk2 < .05){
        cat(", the Chi-Square test statistic with ", md$dfk2, " df is equal to ", round(md$csk2,3), " (p = ", round(md$pvk2,3) , "). This indicates a poor fit for the model with both k's fixed to ",md$k[1],". \n \n", 
            sep = "")
      } else {
        cat(" the Chi-Square test statistic with ", md$dfk2, " df is equal to ", round(md$csk2,3), " (p = ", round(md$pvk2,3) , "). This indicates a good fit for the model with both k's fixed to ",md$k[2],". \n \n", 
            sep = "")
      }
  } else if (md$singleK == FALSE){
    # if the user does not want to calculate a single k: abort the analyses (since the k's cannot be fixed)
    cat("The analyses were aborted since one or both k's could not be constrained to a single value.")
    
  } else if ((is.na(md$k[1]) | is.na(md$k[2]))  &&  md$kdif.yn == 0 && md$singleK == TRUE) {
      cat("\n \n       Since k1 and k2 do not differ significantly (95% CI [", round(md$kdif.l,2),", ", round(md$kdif.u,2), " a new model with a single k for both ", 
          md$lab$lab1, " and ", md$lab$lab2, ". In this new model k is estimated to be ", round(md$k,2),
          ", with its percentile bootstrap confidence interval ranging from ", round(md$kl,2), " to ", round(md$ku,2),
          ". Therefore, it can be concluded that ", kcon(md$kl,md$ku), "\n",    
          sep="")
      if (!is.null(md$k_fix)){
        cat(" When investigating the relative fit of the model with both k's fixed to", md$k_fix)
        if (md$pvk < .05){
          cat("the Chi-Square test statistic with ", md$dfk, " df is equal to ", round(md$csk,3), " (p = ", round(md$pvk,3) , "). This indicates a poor fit for the model with both k's fixed to ",k_refit1(md$kl,md$ku),".", 
              sep = "")
        } else {
          cat(", the Chi-Square test statistic with ", md$dfk, " df is equal to ", round(md$csk,3), " (p = ", round(md$pvk,3) , "). This indicates a good fit for the model with both k's fixed to ",k_refit1(md$kl,md$ku),".", 
              sep = "")
        }
      }
    }
  }
}





# Text with the reliablity of correction run done.
if (md$relok > 0)
{
  out2<-capture.output(summary(md$fit3))
  cat("\n \n \n \n \n####################################################################",
      "#                  Run Correcting for Unreliability                #",
      "#################################################################### \n",
      out2, 
      sep="\n")
}



cat("\n #################################################################### \n",
    "#          DataToText: Run Correcting for Unreliability            # \n",
    "#################################################################### \n \n",
    "    This section describes the lavaan run which corrects for unreliability.")

# Only output the text if the model converge
if (md$fit3@Fit@iterations > 0L && !md$fit3@Fit@converged){
  cat("\n       ERROR: model did not converge \n \n",sep="")
} else{
  cat("\n       The lavaan model converged after ", md$fit3@Fit@iterations, " iterations.",sep="")
  
  md$txt7 = md$texrel
  if (md$reloor == 0) {
    
    cat(md$txt7," The actor effect of the ", md$lab$lab1,  " is equal to ",rround(md$act13,3)," (",
        pval(md$act1p3), ", 95% CI [", round(md$act1.CIl3,2),", ", round(md$act1.CIu3,2), "])",
        ". The overall standardized actor effect is ",rrround(md$act1s3,3),".",
        " The actor effect of the ", md$lab$lab2,  " is equal to ",rround(md$act23,3)," and is",knt(md$act2p3,md$alpha)," statistically significant (",
        pval(md$act2p3), ", 95% CI [", round(md$act2.CIl3,2),", ", round(md$act2.CIu3,2), "])",
        ". The overall standardized actor effect is ",rrround(md$act2s3,3),".",
        "\n The partner effect of the ", md$lab$lab1, " is equal to ",rround(md$part13,3)," and is ",knt(md$partp13,md$alpha)," statistically significant (",
        pval(md$partp13), ", 95% CI [", round(md$part1.CIl3,2),", ", round(md$part1.CIu3,2), "])",
        ". Their standardized partner effects equals ", rrround(md$parts13,3),".",
        " For the ", md$lab$lab2, ", the partner effect is equal to ",rround(md$part23,3)," and is",knt(md$partp23,md$alpha)," statistically significant (",
        pval(md$partp23), ", 95% CI [", round(md$part1.CIl3,2),", ", round(md$part1.CIu3,2),"])",
        ". The standardized effect is ", rrround(md$parts23,3), ".",
        sep="")
    
    
    # Use for-loop to create output for both husbands and wifes
    for (i in 1:2){
      if (eval(parse(text=paste("md$esd",i,"3", sep=""))) < 0) {
        cat("The variance of the errors for the", eval(parse(text=paste("md$lab$lab",i, sep=""))), "is negative which likely indicates that the reliabilities are too high.", 
            "Because of the negative error variance, the partial intraclass correlation for",md$varnames$y,"controlling for actor and partner variables is undefined."
            )
        
      }else if (eval(parse(text=paste("md$esd",i,"3", sep=""))) >= 0){
        cat(" For the ", eval(parse(text=paste("md$lab$lab",i, sep=""))), " the variance of the errors is ", round(eval(parse(text=paste("md$esd",i,"3", sep=""))),3),".",
            sep = ""
        )}
    }
    
    
    cat(" The R squared for the ", md$lab$lab1," is ",rrround(md$rsq3a,3),", for the ", md$lab$lab2, " it is ", rrround(md$rsq3b,3),".",
        sep = "")
    
    if (md$rsq3a > 1 |md$rsq3b > 1) cat(" This correlation is out of range and likely one or both of the reliability values is too large.", 
                                        sep = "") 
    
    
    cat(" The partial intraclass correlation for ",md$varnames$y," controlling for actor and partner variables is equal to ",rrround(md$iccYs3,3),
        " and is",knt(md$iccYp3,md$alpha)," statistically significant (",
        pval(md$iccYp3),", 95% CI [", round(md$iccY.CIl3,2),", ", round(md$iccY.CIu3,2), "]).",
        sep = "") 
    if (md$iccYs3^2 > 1) cat(" (This correlation is out of range and likely one or both of the reliability values is too large.)",
                             sep = "") 
    cat(" The intraclass correlation between the actor and partner variables is equal to ",rrround(md$iccXs3,3)," and is",knt(md$iccXp3,md$alpha),
        " statistically significant (",
        pval(md$iccXp3), " 95% CI [", round(md$iccX.CIl3,2),", ", round(md$iccX.CIu3,2), "]).",
        sep = "") 
    if (md$iccXs3^2 > 1) cat(" (This correlation is out of range and likely one or both of the reliability values is too large.)",
                             sep = "")  
  }
  
  
  
  if (md$act1p3 > md$alpha & md$act2p3 > md$alpha & abs(md$act1s3) < .1 & abs(md$act2s3) < .1){
    stop("The standardized actor effect(s) do not have a substantial value. Therefore, estimates of k would be highly unstable and all further analyses are aborted.") 
  } else{
    cat("\n \n       Because the standardized actor effects of both", md$lab$lab1, "and", md$lab$lab2, 
        "are greater than .1 in absolute value and they are statistically significant, k (i.e. the ratio of the partner effect to the actor effect) can be interpreted in the output. "
        )
    
    # DAK I think we still report the two k's even if they do not differ.
    
    cat("Here, k1 represents the effect of ", md$lab$lab1, " on ", md$lab$lab2, " and k2 the effect of ", md$lab$lab2, " on ", md$lab$lab1,
        ". The value of k1 equals ", round(md$k1,2), " and k2 equals ",       
        round(md$k2,2),".  Their 95% percentile bootstrap confidence intervals range from ", round(md$ck1.CIl,2), " to ", round(md$ck1.CIu,2), " and ", round(md$ck2.CIl,2), 
        " to ", round(md$ck1.CIu,2), ", respectively. So, for the ", md$lab$lab1, ", it can be concluded that ", kcon(md$ck1.CIl,md$ck1.CIu),
        " For the ", md$lab$lab2, ", the confidence interval of k2 shows that ", kcon(md$ck2.CIl,md$ck2.CIu), 
        sep="")  
  }
}
md <<- md
}






###### Tables #######
printTables <- function(md){
  
  # Descriptives of the raw data
    tab1 <- matrix(c(md$varnames$x,  md$lab$lab1,rround(md$x1mean,3),rround(md$x1sd,3),rround(md$x1min,3),rround(md$x1max,3),as.numeric(md$x1n),
                     " ",       md$lab$lab2,rround(md$x2mean,3),rround(md$x2sd,3),rround(md$x2min,3),rround(md$x2max,3),as.numeric(md$x2n),
                     md$varnames$y,  md$lab$lab1,rround(md$y1mean,3),rround(md$y1sd,3),rround(md$y1min,3),rround(md$y1max,3),as.numeric(md$y1n),
                     " ",  md$lab$lab2,rround(md$y2mean,3),rround(md$y2sd,3),rround(md$y2min,3),rround(md$y2max,3),as.numeric(md$y2n))
                   ,nrow=4,ncol=7,byrow=TRUE)
    tab1 <- as.data.frame(tab1)
    colnames(tab1) <- c("Variable","Role","Mean","SD","Minimum","Maximum"," n")
    
    cat("\n \n",
        "Table 1: Descriptive Statistics of the Raw Data \n",
        "----------------------------------------------- \n")
    print(tab1)


  # Descriptive Statistics of the FIML Estimates
    tab2 <- matrix(c(md$varnames$x,  md$lab$lab1,rround(md$xmean1,3),rround(sqrt(md$xvar1),3),
                       "     ",   md$lab$lab2,rround(md$xmean2,3),rround(sqrt(md$xvar2),3),
                       md$varnames$y,  md$lab$lab1,rround(md$ymean1,3),rround(sqrt(md$yvar1),3),
                       "      ",  md$lab$lab2,rround(md$ymean2,3),rround(sqrt(md$yvar2),3))
                   ,ncol=4,byrow=TRUE)
    
    tab2 <- as.data.frame(tab2)
    colnames(tab2) <- c("Variable","Role", "Mean" ,"SD")
    cat("\n \n",
        "Table 2: Descriptive Statistics of the FIML Estimates \n,
        ----------------------------------------------------- \n")
    print(tab2)


  # Table of effect estimates
    ttu=paste0(md$CIv,"% CI")
#     tab3 <- matrix(c(
#       "Intercept", md$lab$lab1,rround(md$intr1,3),rround(md$intr1.CIl,3),"to  ",rround(md$intr1.CIu,3),pval2(md$intr1p),"","","",
#       "Actor","",rround(md$act1,3),rround(md$act1.CIl,3),"to  ",rround(md$act1.CIu,3),pval2(md$act1p),rround(md$act1so,3), rround(md$act1ss,3),rround(md$act1r,3),
#       "Partner","",rround(md$part1,3),rround(md$part1.CIl,3),"to  ",rround(md$part1.CIu,3),pval2(md$part1p),rround(md$part1so,3),rround(md$part1ss,3),rround(md$part1r,3),
#       "k","",rround(md$k1,3),rround(md$kl1,3),"to  ",rround(md$ku1,3),"","","","",
#       "Intercept", md$lab$lab2,rround(md$intr2,3),rround(md$intr2.CIl,3),"to  ",rround(md$intr2.CIu,3),pval2(md$intr2p),"","","",
#       "Actor","",rround(md$act2,3),rround(md$act2.CIl,3),"to  ",rround(md$act2.CIu,3),pval2(md$act2p),rround(md$act2so,3),rround(md$act2ss,3),rround(md$act2r,3),
#       "Partner","",rround(md$part2,3),rround(md$part2.CIl,3),"to  ",rround(md$part2.CIu,3),pval2(md$part2p),rround(md$part2so,3),rround(md$part2ss,3),rround(md$part2r,3),
#       "k","",rround(md$k2,3),rround(md$kl2,3),"to  ",rround(md$ku2,3),"","","","")
#       ,ncol=10,byrow=TRUE)
#     
#     tab3 <- as.data.frame(tab3)
#     colnames(tab3) <- c("Effect", "Role","Estimate" ,"Lower",ttu,"Upper","p value","Beta (o)","Beta (o)","r")
#     cat("\n \n,
#         Table 3: APIM Results for Each Member \n",
#         "------------------------------------- \n")
#     tab3


# Overall APIM Results
  # Not sure about an overall k
  tab4 <- matrix(c(
    "Intercept",rround(md$intro,3),rround(md$intro.CIl,3),"to  ",rround(md$intro.CIu,3),pval2(md$introp),
    "Actor",rround(md$acto,3),rround(md$acto.CIl,3),"to  ",rround(md$acto.CIu,3),pval2(md$actop),
    "Partner",rround(md$parto,3),rround(md$parto.CIl,3),"to  ",rround(md$parto.CIu,3),pval2(md$partop))
    ,ncol=6,byrow=TRUE)
  tab4 <- as.data.frame(tab4)
  colnames(tab4) <- c("Effect", "Estimate" ,"Lower",ttu,"Upper","p value")             
  cat("\n \n",
      "Table 4: Overall APIM Results \n",  
      "----------------------------- \n")
  print(tab4)

  

#   tables <- list("--------------- Table 1: Descriptive Statistics of the Raw Data ---------------" = tab1, 
#               "--------------- Table 2: Descriptive Statistics of the FIML Estimates ---------------" = tab2, 
#               "--------------- Table 3: APIM Results for Each Member ---------------" = tab3,
#               "--------------- Table 4: Overall APIM Results ---------------" = tab4)

  return()
}






###### Citation ######
cite <- function(){
cat("\n \n ****CAUTION****  If you decide to use information contained here in a paper, please make sure that you acknowledge that you have used this program.",
    "Also should you decide to use the exact text included here, you will need to put quotes around that material to avoid plagiarism. Although great effort has been undertaken to ensure the accuracy of results,",
    "no complete guarantee can be made about their accuracy.  It is your responsibility to check the results and text for accuracy.  If you do find an error, please report it to Lara Stas (Lara.Stas@UGent.be) or David A. Kenny (david.kenny@uconn.edu).",
    sep = "")
}






###### Figures ######
figures <- function(temploc){
  ### Figure 1: Standard Model
  figgu=paste0(temploc,'APIM.tif')
  FigU=paste0(temploc,'APIM_SEM_Fig.tif')
  ttu=""
  
  if (md$act1p<md$alpha) ttu="*"
  aaa1=paste0(rround(md$act1,3),ttu," (",rrround(md$act1so,3),")")
  ttu=""
  if (md$act2p<md$alpha) ttu="*"
  aaa2=paste0(rround(md$act2,3),ttu," (",rrround(md$act2so,3),")")
  ttu=""
  
  if (md$part1p<md$alpha) ttu="*"
  ppp1=paste0(rround(md$part1,3),ttu," (",rrround(md$part1so,3),")")
  ttu=""
  if (md$part2p<md$alpha) ttu="*"
  ppp2=paste0(rround(md$part2,3),ttu," (",rrround(md$part2so,3),")")
  ttu=""
  
  if (md$iccXp<md$alpha) ttu="*"
  rxx=paste0(rround(md$iccX,3),ttu," (",rrround(md$iccXs,3),")")
  ttu=""
  if (md$iccYp<md$alpha) ttu="*"
  ree=paste0(rround(md$iccY,3),ttu," (",rrround(md$iccYs,3),")")
  
  img <- readTiff(figgu)
  tiff(filename =FigU, width = 1400, height = 837,   pointsize = 14, bg = "white")
  
  plot(img)
  
  text(224,170,md$varnames$x,cex=1.8)
  text(1015,170,md$varnames$y,cex=1.8)
  text(205,676,md$varnames$x,cex=1.8)
  text(1020,679,md$varnames$y,cex=1.8)
  text(224,140,md$p2,cex=1.8)
  text(1015,140,md$p2,cex=1.8)
  text(205,646,md$p1,cex=1.8)
  text(1020,650,md$p1,cex=1.8)
  # Small change to make figure pretty
  text(124,412,rxx,cex=1.8)
  text(1263,412,ree,cex=1.8)
  text(1005,300,ppp2,cex=1.8)
  text(405,300,ppp1,cex=1.8)
  text(645,180,aaa2,cex=1.8)
  text(645,690,aaa1,cex=1.8)
  text(1290,278,"E2",cex=1.8)
  text(1290,560,"E1",cex=1.8)
  text(1290,278,"E2",cex=1.8)
  text(1290,560,"E1",cex=1.8)
  title(main = "Actor-Partner Interdependence Model (standardized estimates)",cex.main=3.2,line = -3)
  dev.off()



  ### Figure 2: figure with correlations corrected for unreliability.
  if (md$relok > 0 & md$reloor==0)
  {
    figgy=paste0(temploc,'APIML.tif')
    FigL=paste0(temploc,'APIM_SEMLV_Fig.tif')
    ttu=""
    
    if (md$act1p3<md$alpha) ttu="*"
    aaa1=paste0(rround(md$act13,3),ttu," (",rrround(md$act1s3,3),")")
    ttu=""
    if (md$act2p3<md$alpha) ttu="*"
    aaa2=paste0(rround(md$act23,3),ttu," (",rrround(md$act2s3,3),")")
    ttu=""
    
    if (md$partp13<md$alpha) ttu="*"
    ppp1=paste0(rround(md$part13,3),ttu," (",rrround(md$parts13,3),")")
    ttu=""
    if (md$partp23<md$alpha) ttu="*"
    ppp2=paste0(rround(md$part23,3),ttu," (",rrround(md$parts23,3),")")
    ttu=""
    
    if (md$iccXp3<md$alpha) ttu="*"
    rxx=paste0(rround(md$iccX3,3),ttu," (",rrround(md$iccXs3,3),")")
    ttu=""
    if (md$iccYp3<md$alpha) ttu="*"
    ree=paste0(rround(md$iccY3,3),ttu," (",rrround(md$iccYs3,3),")")
    
    
    if (md$esd13 <= 0 | md$esd13 <= 0) ree=paste0(rround(md$iccY3,3),ttu," (undef.)")
    
    img <- readTiff(figgy)
    tiff(filename =FigL, width = 1400, height = 837,   pointsize = 14, bg = "white")
    
    plot(img)
    text(224,155,md$varnames$x,cex=1.8)
    text(1015,165,md$varnames$y,cex=1.8)
    text(217,701,md$varnames$x,cex=1.8)
    text(1020,694,md$varnames$y,cex=1.8)
    text(224,125,md$p2,cex=1.8)
    text(1015,135,md$p2,cex=1.8)
    text(217,671,md$p1,cex=1.8)
    text(1020,664,md$p1,cex=1.8)
    # DAK small change to make figure pretty
    text(124,412,rxx,cex=1.8)
    text(1273,412,ree,cex=1.8)
    text(985,300,ppp2,cex=1.8)
    text(528,300,ppp1,cex=1.8)
    text(645,165,aaa2,cex=1.8)
    text(645,700,aaa1,cex=1.8)
    text(1290,278,"E2",cex=1.8)
    text(1290,560,"E1",cex=1.8)
    text(1290,278,"E2",cex=1.8)
    text(1290,560,"E1",cex=1.8)
    title(main = "Actor-Partner Interdependence Model (standardized estimates)",cex.main=3.2,line = -3)
    dev.off()
  }
}



##### CHECK DEZE #####
# CHECK OF DEZE GEBRUIKT #
# DAK this function is for text output; may not be needed
effectis  = function (eff,effp,effs,alp){
  resu=paste0(" is equal to ",rround(eff,3)," and is",knt(effp,alp)," statistically significant (",pval(effp),") with a standardized effect of ",rrround(effs,3))
}
# CHECK OF KCON NIET BETER IN DE TEKST FUNCTIE ZIT



#' APIM.SEM
#'
#' @name APIM.SEM
#' @docType package
NULL


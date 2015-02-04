##### Helper functions #####



# Helper function for reporting about the effect size 
effsz = function (eff, dic){ 
  esmall=.1
  emed =.3
  elarge=.5
  if (dic == 1) {esmall=.2}
  if (dic == 1) {emed =.5}
  if (dic == 1) {elarge =.8}
  btto =" less than small"
  if (abs(eff)>esmall) btto=" and a small effect size"
  if (abs(eff)>emed) btto=" and a medium effect size"
  if (abs(eff)>elarge) btto=" and a large effect size"
  return(btto) 
}

effszm = function (dic){  
  ess=' (r = '
  if (dic == 1)ess=' (d = '
  return(ess)
}

# Helper function: rounding values 
rround  = function (num,dig){
  #DAK added two more options; not needed here but used for other programs.
  if (dig==5) return(sprintf("%.3f", round(num,5)))
  if (dig==4) return(sprintf("%.3f", round(num,4)))
  
  if (dig==2) return(sprintf("%.2f", round(num,2)))
  if (dig==3) return(sprintf("%.3f", round(num,3)))
  if (dig==1) return(sprintf("%.1f", round(num,1)))
}
rrround  = function (num,dig){
  if (substring(as.character(rround(num,dig)),1,2)=="0.") return(substring(as.character(rround(num,dig)),2,5))
  if (substring(as.character(rround(num,dig)),1,2)=="-0") return(gsub("-0","-", rround(num,dig)))
  if (substring(as.character(rround(num,dig)),1,2)=="1.") return(substring(as.character(rround(num,dig)),1,5))
  if (substring(as.character(rround(num,dig)),1,2)=="-1") return(substring(as.character(rround(num,dig)),1,6))
}

# Helper function: p-values 
pval = function(num) { 
  if (num>=.001) xxx= paste0("p = ",rrround(num,3))
  if (num<.001) xxx= "p < .001"
  return(xxx)}
pval2 = function(num) { 
  if (num>=.001) xxx= rrround(num,3)
  if (num<.001) xxx= "<.001"
  return(xxx)}
pval3 = function(pvv,alph) { 
  xxx=""
  if (pvv < alph) xxx="*"
  return(xxx)}

# report the p-values 
knt = function(pvv,alph) { 
  xxx=""
  if (pvv >= alph) xxx=" not"
  return(xxx)}
pvalchi = function(num,dff,pvv) {
  if (pvv>=.001) xxx= paste0(", p = ",rrround(pvv,3),")")
  if (pvv<.001) xxx= ", p < .001)"
  xxx=paste0("(chi-square(",dff,") = ",rround(num,2),xxx) 
  return(xxx)}
onetailp = function(p) {  
  if (p <  .5) { 
    px=p/2
    return (px) }
  px = (1 + p)/2 
  return (px)} 

vlisty = function(nn,vname) 
{  
  vvlist=""
  for (ii in 1:nn)
  { 
    if (nn > 2 & ii < nn - 1) vvlist = paste0(vvlist, vname[ii],", ") 
    if (nn > 2  & ii  == nn -1) vvlist = paste0(vvlist,vname[ii],", and ") 
    if (ii == 1 & nn ==2) vvlist = paste0(vvlist, vname[ii]," and ") 
    if (ii == nn) vvlist=paste0(vvlist,vname[ii])
  }
  return(vvlist)
}

# Report about k. 
kcon = function (loww,upp) {
  if (loww <  -1.0 & upp > 1.0) 
  { txtx="the confidence interval for k is very wide and it cannot be determined what model is the most likely."
    return(txtx) } 
  if (loww <  0.0 & upp > 1.0) 
  { txtx="the contrast model (k = -1) is implausible and that the couple (k = 1) and the actor-only models (k = 0) are plausible."
    return(txtx) } 
  if (loww <  -1.0 & upp > 0.0) 
  { txtx="the couple model (k = 1) is implausible and that the contrast (k = -1) and the actor-only models (k = 0) are plausible."
    return(txtx)} 
  if (loww <  1 & upp > 1) 
  { txtx="the contrast (k = -1) and the actor-only (k = 0) models are implausible and that the couple model (k = 1) is plausible."
    return(txtx)} 
  if (loww <  0 & upp > 0) 
  { txtx="the contrast (k = -1) and the couple (k = 1) models are implausible and that the actor-only model (k = 0) is plausible."
    return(txtx)} 
  if (loww <  -1 & upp > -1) 
  { txtx="the couple (k = 1) and the actor-only (k = 0) models are implausible and that the contrast model (k = -1) is plausible."
    return(txtx)} 
  if (loww > 0 & upp < 1) 
  { txtx="the model is in between the actor-only (k = 0) and the couple (k = 1) models."
    return(txtx)} 
  if (loww > -1 & upp < 0) 
  { txtx="the model is in between the contrast (k = -1) and the actor-only (k = 0) models."
    return(txtx)} 
  if (upp < -1) 
  { txtx="the model is more extreme than the contrast model (k = -1)."
    return(txtx)} 
  if (loww > 1) 
  { txtx="the model is more extreme than the couple model (k = 1)."
    return(txtx)} 
}

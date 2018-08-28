#============================================================================
# Name        : Statistics
# Author      : Apoorv Saxena
# Version     :
# Copyright   : 
# Description : performance analysis
#============================================================================
source('Systemdefinition.R');
source('Solverfunctions.R');
source('matrixcoefficients.R');
source('performance measures.R');
#source('C:/Users/sapoorv/Downloads/latex/Data Backups/Code/Optimize.R');

# DEFINE THE SYTEM VARIABLES HERE
type = "mod"; #low, mod,high
if(type == "high")
{
  load = 0.99;
}else if( type =="mod")
{
  load = 0.6;
}else
{
  load = 0.2;
}
c=10;
th = 30;
alpha = {}; #these are the probabilites of start based on the system content.
for(i in 1:(th-1))
{
  alpha[i] = i/(3*c); #equally likely to start at every pointprint(lambda);
  #alpha[i]=0.5+i*i/(2*th*th);
}
len=20;
vacation_len = seq(16, 45, length.out=len);

# Computation of the system (queue+server) content based on the load of the system.
ESystCon={};
EServerCon={};
EQCon={};
EServerConinBatc={};
ProbServerBusi={};
EQConVacatio={};
rpoUBv={};
NumVacationFPv={};
NewService = {};
sdrpol = {};
sdrpou = {};
sdrpoapprox = {};
sdrrpo_new = {}
lambda = 0.51#load*c/ServiceGFDerivative(1,c) ; #0.35
#ArrivalGFDerivative(1,lambda)*ServiceGFDerivative(1,c)/c;
for(j in 1:len)
{
  vacation_probabilities=function()
  {
    vacat_len=vacation_len[j];
    #y = floor(vacat_len/4);
    temp = rep_len(0,vacat_len+3);
    # for(j in (vacat_len+2-y):(vacat_len+1))
    # {
    #   temp[j] = 1/(y);
    # }
    temp[vacat_len-1] = 0.2
    temp[vacat_len]   = 0.2
    temp[vacat_len+1] = 0.2
    temp[vacat_len+2] = 0.2
    temp[vacat_len+3] = 0.2
    vacation_probabilities=temp;
  }  
  
  Consts =              Constants(c,lambda,alpha,th);
  ESystCon[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));  
  EServerCon[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
  EQCon[j]=             ESystCon[j]-EServerCon[j];
  EServerConinBatc[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
  ProbServerBusi[j]=    ProbServerBusy(th,alpha,Consts);
  EQConVacatio[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
  rpoUBv[j]=             rpoUB(th,lambda,alpha);
  NumVacationFPv[j]=    NumVacationFP(th,lambda,alpha);
  temp=                 Varrpob(th,lambda,alpha);
  NewService[j]=        StartNewService(th,lambda,alpha,Consts);
  sdrpou[j]=           sqrt(temp[2]);
  sdrpol[j]=           sqrt(temp[1]);
  sdrpoapprox[j]=        temp[3];
  sdrrpo_new[j]=        sdv_age_nu(th,lambda,alpha);
}  

for(i in 1:(th-1))
{
  alpha[i]=0
}

ESystCon2={};
EServerCon2={};
EQCon2={};
EServerConinBatc2={};
ProbServerBusi2={};
EQConVacatio2={};
rpoUBv2={};
NumVacationFPv2={};
NewService2 ={};
evacat={};
sdrpol2 = {};
sdrpou2 = {};
sdrpoapprox2={};
sdrrpo_new2 = {};
for(j in 1:(len))
{
  vacation_probabilities=function()
  {
    vacat_len=vacation_len[j];
    #y = floor(vacat_len/4);
    temp = rep_len(0,vacat_len+3);
    # for(j in (vacat_len+2-y):(vacat_len+1))
    # {
    #   temp[j] = 1/(y);
    # }
    temp[vacat_len-1] = 0.2
    temp[vacat_len]   = 0.2
    temp[vacat_len+1] = 0.2
    temp[vacat_len+2] = 0.2
    temp[vacat_len+3] = 0.2
    vacation_probabilities=temp;
  } 
  evacat[j]=             VacationGFDerivative(1);
  Consts =               Constants(c,lambda,alpha,th);
  ESystCon2[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
  EServerCon2[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
  EQCon2[j]=             ESystCon2[j]-EServerCon2[j];
  EServerConinBatc2[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
  ProbServerBusi2[j]=    ProbServerBusy(th,alpha,Consts);
  EQConVacatio2[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
  rpoUBv2[j]=            rpoUB(th,lambda,alpha);
  NumVacationFPv2[j]=    NumVacationFP(th,lambda,alpha);
  temp=                 Varrpob(th,lambda,alpha);
  NewService2[j]=       StartNewService(th,lambda,alpha,Consts);
  sdrpou2[j]=           sqrt(temp[2]);
  sdrpol2[j]=           sqrt(temp[1]);
  sdrpoapprox2[j]=        temp[3];
  sdrrpo_new2[j]=        sdv_age_nu(th,lambda,alpha);
}
for(i in 1:(th-1))
{
  alpha[i]= 1;
}
ESystCon3={};
EServerCon3={};
EQCon3={};
EServerConinBatc3={};
ProbServerBusi3={};
EQConVacatio3={};
rpoUBv3={};
NumVacationFPv3={};
NewService3={};
sdrpol3 = {};
sdrpou3 = {};
sdrpoapprox3={};
sdrrpo_new3 = {}

for(j in 1:(len))
{
  vacation_probabilities=function()
  {
    vacat_len=vacation_len[j];
    #y = floor(vacat_len/4);
    temp = rep_len(0,vacat_len+3);
    # for(j in (vacat_len+2-y):(vacat_len+1))
    # {
    #   temp[j] = 1/(y);
    # }
    temp[vacat_len-1] = 0.2
    temp[vacat_len]   = 0.2
    temp[vacat_len+1] = 0.2
    temp[vacat_len+2] = 0.2
    temp[vacat_len+3] = 0.2
    vacation_probabilities=temp;
  } 
  Consts =               Constants(c,lambda,alpha,th);
  ESystCon3[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
  EServerCon3[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
  EQCon3[j]=             ESystCon3[j]-EServerCon3[j];
  EServerConinBatc3[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
  ProbServerBusi3[j]=    ProbServerBusy(th,alpha,Consts);
  EQConVacatio3[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
  rpoUBv3[j]=            rpoUB(th,lambda,alpha);
  NumVacationFPv3[j]=    NumVacationFP(th,lambda,alpha);
  temp=                 Varrpob(th,lambda,alpha);
  NewService3[j]=       StartNewService(th,lambda,alpha,Consts);
  sdrpou3[j]=           sqrt(temp[2]);
  sdrpol3[j]=           sqrt(temp[1]);
  sdrpoapprox3[j]=        temp[3];
  sdrrpo_new3[j]=        sdv_age_nu(th,lambda,alpha);
}

system('taskkill /f /im AcroRd32.exe');
#####
#remove inf and nan from the vectors
evacat = evacat[is.finite(evacat)];
ESystCon= ESystCon[is.finite(ESystCon)];
ESystCon2 = ESystCon2[is.finite(ESystCon2)];
ESystCon3= ESystCon3[is.finite(ESystCon3)];
EServerCon= EServerCon[is.finite(EServerCon)];
EServerCon2= EServerCon2[is.finite(EServerCon2)];
EServerCon3= EServerCon3[is.finite(EServerCon3)];
EQCon= EQCon[is.finite(EQCon)];
EQCon2= EQCon2[is.finite(EQCon2)];
EQCon3= EQCon3[is.finite(EQCon3)];
EServerConinBatc= EServerConinBatc[is.finite(EServerConinBatc)];
EServerConinBatc2= EServerConinBatc2[is.finite(EServerConinBatc2)];
EServerConinBatc3= EServerConinBatc3[is.finite(EServerConinBatc3)];
ProbServerBusi= ProbServerBusi[is.finite(ProbServerBusi)];
ProbServerBusi2= ProbServerBusi2[is.finite(ProbServerBusi2)];
ProbServerBusi3= ProbServerBusi3[is.finite(ProbServerBusi3)];
EQConVacatio= EQConVacatio[is.finite(EQConVacatio)];
EQConVacatio2= EQConVacatio2[is.finite(EQConVacatio2)];
EQConVacatio3= EQConVacatio3[is.finite(EQConVacatio3)];
rpoUBv= rpoUBv[is.finite(rpoUBv)];
rpoUBv2= rpoUBv2[is.finite(rpoUBv2)];
rpoUBv3= rpoUBv3[is.finite(rpoUBv3)];
NumVacationFPv= NumVacationFPv[is.finite(NumVacationFPv)];
NumVacationFPv2= NumVacationFPv2[is.finite(NumVacationFPv2)];
NumVacationFPv3= NumVacationFPv3[is.finite(NumVacationFPv3)];
NewService = NewService[is.finite(NewService)];
NewService2 = NewService2[is.finite(NewService2)];
NewService3 = NewService3[is.finite(NewService3)];
sdrpol = sdrpol[is.finite(sdrpol)];
sdrpol2 = sdrpol2[is.finite(sdrpol2)];
sdrpol3 = sdrpol3[is.finite(sdrpol3)];
sdrpou = sdrpou[is.finite(sdrpou)];
sdrpou2 = sdrpou2[is.finite(sdrpou2)];
sdrpou3 = sdrpou3[is.finite(sdrpou3)];
sdrpoapprox = sdrpoapprox[is.finite(sdrpoapprox)];
sdrpoapprox2 = sdrpoapprox2[is.finite(sdrpoapprox2)];
sdrpoapprox3 = sdrpoapprox3[is.finite(sdrpoapprox3)];
saveGraph=function(
  variable_name,
  type,
  y_values1,
  y_values2,
  y_values3,
  x_values,
  y_axis_name,
  x_axis_name,
  min_offset,
  max_offset,
  legend_position
)
{
  filename = q=paste("C:/Users/sapoorv/Downloads/latex/Data Backups/figs/plot_", variable_name,"_vacation_",type,".pdf", sep = "");
  pdf(filename);
  
  plot_colors = c("blue", "green", "red");
  names = c(expression(alpha[i] == paste(i/3*c)) , expression(alpha[i] == paste(0)), expression(alpha[i] == paste(1)));
  
  ylim2 = max(y_values1,y_values2,y_values3 );
  ylim2 = ylim2 + max_offset;
  
  ylim1 = min(y_values1,y_values2,y_values3);
  ylim1 = max(0, ylim1 - min_offset)
  ylim = c(ylim1,ylim2);
  if( max_offset < 0)
  {
    ylim = c(ylim1,0);
  }
  if( min_offset < 0)
  {
    ylim = c(0,ylim2);
  }
  par(mar=c(5,5,5,2)); 
  plot(x_values,y_values1, type = "l", xlab = x_axis_name, ylab = y_axis_name,col="blue", cex.lab=1.5, lwd=3.5,ylim=ylim );
  lines(x_values,y_values2,col="green",type="l", lty=2, lwd=2);
  lines(x_values,y_values3,type="l", lty=3, lwd=3.5, col="red");
  legend(legend_position, names, cex=1.4, col=plot_colors, lty=1:3, lwd=3.5, bty="n");
  dev.off();
}
####
#System content at random slot
saveGraph("syst_cont", type, ESystCon,ESystCon2,ESystCon3,evacat, 'E[System Content]', 'E[Vacation Length]', 25, 25, "topleft")

#Queue content at random slot
saveGraph("queue_cont", type, EQCon,EQCon2,EQCon3,evacat, 'E[Backlog Size]', 'E[Vacation Length]', 0, 0, "topleft")

# # Server content in a random batch
saveGraph("serv_cont_rand_batch", type, EServerConinBatc, EServerConinBatc2, EServerConinBatc3,evacat, 'E[Server Content in Random Batch]', 'E[Vacation Length]', 0, 1, "topleft")

# # Queue Content when server is in vacation
saveGraph("queue_cont_vac", type, EQConVacatio, EQConVacatio2, EQConVacatio3,evacat, 'E[Queue Content in Vacation]', 'E[Vacation Length]', -1, 0, "topleft");

# # Probability that server is busy
saveGraph("prob_server_busy", type, ProbServerBusi, ProbServerBusi2, ProbServerBusi3,evacat, 'Probability Server is busy', 'E[Vacation Length]', 0.1, 0.1, "topleft");

## Num vacations
saveGraph("num_vacat_fp", type, NumVacationFPv, NumVacationFPv2, NumVacationFPv3,evacat, 'E[Number of vacations for first packet]', 'E[Vacation Length]', -1, 0, "topright");

## expected age nu
saveGraph("rpo", type, rpoUBv,rpoUBv2,rpoUBv3,evacat, expression(E(Age[nu])), 'E[Vacation Length]', 0, 0.5, "bottomright");

## Probability of start of new service
saveGraph("prob_new_service", type, NewService,NewService2,NewService3,evacat, 'Pr( New Service at Service Initiation Opportunity )', 'E[Vacation Length]', 0, 0.01, "bottomright");

## standart deviation  of rpo
if( prob_poisson()==1)
{
  saveGraph("sd_rpo", type,sdrrpo_new,sdrrpo_new2,sdrrpo_new3,evacat, 'SD[Age_nu]', 'E[Vacation Length]', 0, 0.5, "topleft");
}
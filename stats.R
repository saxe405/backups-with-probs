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

type = "mod"; #low,mod,high
# vacation_probabilities=function()
# {
#   if(type == "mod")
#   {
#     vacat_len=20;
#     y = floor(vacat_len/4);
#     temp = rep_len(0,vacat_len+1);
#     for(j in (vacat_len+2-y):(vacat_len+1))
#     {
#       temp[j] = 1/(y);
#     }
#     vacation_probabilities=temp;
#    # vacation_probabilities = c(0,1)
#   }
#   else if(type =="low")
#   {
#     vacation_probabilities=c(0,0,0,0,0,0.1,0.7,0.1,0.1); #low
#   }  
#   else
#   {  
#     #high
#     vacat_len = 25;
#     temp =rep_len(0, vacat_len+1);
#     temp[vacat_len-2]=0.1;
#     temp[vacat_len-1]=0.1;
#     temp[vacat_len]=0.7;
#     temp[vacat_len+1]=0.1;
#     vacation_probabilities= temp;
#   }
# }  
# DEFINE THE SYTEM VARIABLES HER
c=10;
th =30; #threshold denoted by l in notebook.
offset = 0;
alpha = {}; #these are the probabilites of start based on the system content.
for(i in 1:(th-1))
{
  #alpha[i] = 0.5+i*i/(2*th*th); #equally likely to start at every pointprint(lambda);
  alpha[i]= i/(3*c);
}
len=20;
#load=seq(0.2, 0.85, length.out=len);

# Computation of the system (queue+server) content based on the load of the system.
ESystCon={};
EServerCon={};
EQCon={};
EServerConinBatc={};
ProbServerBusi={};
EQConVacatio={};
rpoUBv = {};
NumVacationFPv = {};
EAgenu2 = {}
NewService = {};
sdrpol = {};
sdrpou = {};
sdrpoapprox = {};
sdrrpo_new = {}
poissonarrivalrates = seq(0.2, 0.70, length.out = len)
load = {}
for(j in 1:(len-offset))
{
  lambda =              poissonarrivalrates[j];
  gamma =               gammaparam();
  load[j]=              ArrivalGFDerivative(1,lambda)*ServiceGFDerivative(1,c)/c;
  Consts =              Constants(c,lambda,alpha,th);
  ESystCon[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
  EServerCon[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
  EQCon[j]=             ESystCon[j]-EServerCon[j];
  EServerConinBatc[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
  ProbServerBusi[j]=    ProbServerBusy(th,alpha,Consts);
  EQConVacatio[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
  rpoUBv[j]=            rpoUB(th,lambda,alpha);
  NumVacationFPv[j]=    NumVacationFP(th,lambda,alpha);
  EAgenu2[j] = expected_Agenu2(th,lambda,alpha);
  temp=                 Varrpob(th,lambda,alpha);
  NewService[j]= StartNewService(th,lambda,alpha,Consts);
  sdrpou[j]=           sqrt(temp[2]);
  sdrpol[j]=           sqrt(temp[1]);
  sdrpoapprox[j]=      temp[3];
  sdrrpo_new[j] = sdv_age_nu(th,lambda,alpha);
}
#print(sdrrpo_new)
# for(j in 1:(len-offset))
# {
#   lambda =              load[j]*c/EServicet(c);
#   Consts =              Constants(c,lambda,alpha,th);
#   
#   ESystCon[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));  
#   EServerCon[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
#   EQCon[j]=             ESystCon[j]-EServerCon[j];
#   EServerConinBatc[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
#   ProbServerBusi[j]=    ProbServerBusy(th,alpha,Consts);
#   EQConVacatio[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
#   rpoUBv[j]=            rpoUB(th,lambda,alpha);
#   NumVacationFPv[j]=    NumVacationFP(th,lambda,alpha);
#   EAgenu2[j] = expected_Agenu2(th,lambda,alpha);
#   temp=                 Varrpob(th,lambda,alpha);
#   NewService[j]= StartNewService(th,lambda,alpha,Consts);
#   sdrpou[j]=           sqrt(temp[2]);
#   sdrpol[j]=           sqrt(temp[1]);
#   sdrpoapprox[j]=      temp[3];
#   sdrrpo_new[j] = sdv_age_nu(th,lambda,alpha);
# }  

for(i in 1:(th-1))
{
  alpha[i]=0;
}
#th=c;
ESystCon2={};
EServerCon2={};
EQCon2={};
EServerConinBatc2={};
ProbServerBusi2={};
EQConVacatio2={};
rpoUBv2 ={};
NumVacationFPv2={};
EAgenu22 = {}
NewService2={};
sdrpol2={};
sdrpou2={};
sdrpoapprox2 = {};
sdrrpo_new2 = {}


for(j in 1:(len-offset))
{
  lambda =              poissonarrivalrates[j];
  Consts =              Constants(c,lambda,alpha,th);
  ESystCon2[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
  EServerCon2[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
  EQCon2[j]=             ESystCon2[j]-EServerCon2[j];
  EServerConinBatc2[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
  ProbServerBusi2[j]=    ProbServerBusy(th,alpha,Consts);
  EQConVacatio2[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
  rpoUBv2[j]=            rpoUB(th,lambda,alpha);
  NumVacationFPv2[j]=    NumVacationFP(th,lambda,alpha);
  EAgenu22[j] = expected_Agenu2(th,lambda,alpha);
  temp=                  Varrpob(th,lambda,alpha);
  NewService2[j]=       StartNewService(th,lambda,alpha,Consts);
  sdrpou2[j]=           sqrt(temp[2]);
  sdrpol2[j]=           sqrt(temp[1]);
  sdrpoapprox2[j]=      temp[3];
  sdrrpo_new2[j] = sdv_age_nu(th,lambda,alpha)
}

# for(j in 1:(len-offset))
# {
#   lambda =              load[j]*c/EServicet(c);
#   Consts =              Constants(c,lambda,alpha,th);
#   ESystCon2[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
#   EServerCon2[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
#   EQCon2[j]=             ESystCon2[j]-EServerCon2[j];
#   EServerConinBatc2[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
#   ProbServerBusi2[j]=    ProbServerBusy(th,alpha,Consts);
#   EQConVacatio2[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
#   rpoUBv2[j]=            rpoUB(th,lambda,alpha);
#   NumVacationFPv2[j]=    NumVacationFP(th,lambda,alpha);
#   EAgenu22[j] = expected_Agenu2(th,lambda,alpha);
#   temp=                  Varrpob(th,lambda,alpha);
#   NewService2[j]=       StartNewService(th,lambda,alpha,Consts);
#   sdrpou2[j]=           sqrt(temp[2]);
#   sdrpol2[j]=           sqrt(temp[1]);
#   sdrpoapprox2[j]=      temp[3];
#   sdrrpo_new2[j] = sdv_age_nu(th,lambda,alpha)
# }

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
EAgenu23 = {}
NewService3 = {};
sdrpou3 ={};
sdrpol3 ={};
sdrpoapprox3={};
sdrrpo_new3 = {}

for(j in 1:(len-offset))
{
  lambda =              poissonarrivalrates[j];
  Consts =              Constants(c,lambda,alpha,th);
  ESystCon3[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
  EServerCon3[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
  EQCon3[j]=             ESystCon3[j]-EServerCon3[j];
  EServerConinBatc3[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
  ProbServerBusi3[j]=    ProbServerBusy(th,alpha,Consts);
  EQConVacatio3[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
  rpoUBv3[j]=            rpoUB(th,lambda,alpha);
  NumVacationFPv3[j]=    NumVacationFP(th,lambda,alpha);
  EAgenu23[j] = expected_Agenu2(th,lambda,alpha);
  temp=                  Varrpob(th,lambda,alpha);
  NewService3[j]=       StartNewService(th,lambda,alpha,Consts);
  sdrpou3[j]=           sqrt(temp[2]);
  sdrpol3[j]=           sqrt(temp[1]);
  sdrpoapprox3[j]=      temp[3];
  sdrrpo_new3[j] = sdv_age_nu(th,lambda,alpha)
}

# for(j in 1:(len-offset))
# {
#   lambda =              load[j]*c/EServicet(c);
#   Consts =              Constants(c,lambda,alpha,th);
#   ESystCon3[j]  =        Re(ESystCont(c,lambda,alpha,Consts,th));
#   EServerCon3[j]=        Re(EServerCont(c,lambda,alpha,Consts,th));
#   EQCon3[j]=             ESystCon3[j]-EServerCon3[j];
#   EServerConinBatc3[j]=  Re(EServerContinBatch(c,lambda,alpha,Consts,th));
#   ProbServerBusi3[j]=    ProbServerBusy(th,alpha,Consts);
#   EQConVacatio3[j]=      Re(EQConVacation(th,lambda,alpha,Consts));
#   rpoUBv3[j]=            rpoUB(th,lambda,alpha);
#   NumVacationFPv3[j]=    NumVacationFP(th,lambda,alpha);
#   EAgenu23[j] = expected_Agenu2(th,lambda,alpha);
#   temp=                  Varrpob(th,lambda,alpha);
#   NewService3[j]=       StartNewService(th,lambda,alpha,Consts);
#   sdrpou3[j]=           sqrt(temp[2]);
#   sdrpol3[j]=           sqrt(temp[1]);
#   sdrpoapprox3[j]=      temp[3];
#   sdrrpo_new3[j] = sdv_age_nu(th,lambda,alpha)
# }
system('taskkill /f /im AcroRd32.exe');



#System content at random slot
saveGraph("syst_cont", type, ESystCon,ESystCon2,ESystCon3,load[1:(len-offset)], 'E[System Content]', 'load', 0, 0, "topleft")

#Queue content at random slot
saveGraph("queue_cont", type, EQCon,EQCon2,EQCon3,load[1:(len-offset)], 'E[Backlog Size]', 'load', -1, 0, "topleft")

# # Server content in a random batch
saveGraph("serv_cont_rand_batch", type, EServerConinBatc, EServerConinBatc2, EServerConinBatc3,load[1:(len-offset)], 'E[Server Content in Random Batch]', 'load', 0, 1, "topleft")

# # Queue Content when server is in vacation
saveGraph("queue_cont_vac", type, EQConVacatio, EQConVacatio2, EQConVacatio3,load[1:(len-offset)], 'E[Queue Content in Vacation]', 'load', 0, 1, "topleft");

# # Probability that server is busy
saveGraph("prob_server_busy", type, ProbServerBusi, ProbServerBusi2, ProbServerBusi3,load[1:(len-offset)], 'Probability Server is busy', 'load', 0.1, 0.1, "topleft");

## Num vacations
saveGraph("num_vacat_fp", type, NumVacationFPv, NumVacationFPv2, NumVacationFPv3,load[1:(len-offset)], 'E[Number of vacations for first packet]', 'load', 0, 0.5, "topright");

## expected age nu
saveGraph("rpo", type, rpoUBv,rpoUBv2,rpoUBv3,load[1:(len-offset)], expression(E(Age[nu])), 'load', 0, 0.5, "topright");

## Probability of start of new service
saveGraph("prob_new_service", type, NewService,NewService2,NewService3,load[1:(len-offset)], 'Pr( New Service at Service Initiation Opportunity )', 'load', 0, 0.01, "bottomright");

## standart deviation  of rpo
if( prob_poisson()==1)
{
  saveGraph("sd_rpo", type,sdrrpo_new,sdrrpo_new2,sdrrpo_new3,poissonarrivalrates[1:(len-offset)], 'SD[Age_nu]', 'load', 0, 0.5, "topright");
}
require(Deriv)
require(BB)
####################################################################
# Expected system connent using the generating functions defined
####################################################################
ESystCont=function(c,lambda,alpha,constants,th)
{  
  ESystCont=0;
  
  d1A=ArrivalGFDerivative(1,lambda);
  d1Gc = ServiceGFDerivative(1,c);
  d2Gc = ServiceGF2d(1,c);
  d1T= VacationGFDerivative(1);
  d2T= vacationGF2d(1);
  d2A = ArrivalGF2d(1,lambda);
  #d1Gm = ServiceGFDerivative(1,m);
  #d2Gm = ServiceGF2d(1,m);
  
  for(j in 1:(c-1))
  {
    m=j;
    d1Gm = ServiceGFDerivative(1,m);
    d2Gm = ServiceGF2d(1,m);
    temp=0;
    temp = (d2T-d2Gm)*d1A/2;
    temp = temp + m*(d1T-d1Gm);
    temp = temp + c*d1Gc*(m+d1T*d1A-d1Gm*d1A)/(c-d1Gc*d1A);
    temp = temp + d2Gc*d1A*(m+d1T*d1A-d1Gm*d1A)/(2*(c-d1Gc*d1A));
    temp = temp + d1Gc*(2*m*d1T*d1A+d2T*(d1A^2)+d1T*d2A + m*(m-1)-d2Gm*(d1A^2)-d1Gm*d2A)/(2*(c-d1Gc*d1A));
    temp = temp - d1Gc*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)*(m+d1T*d1A-d1Gm*d1A)/(2*(c-d1Gc*d1A)^2);
    # temp = temp + ((m+d1T*d1A-d1Gm*d1A)*(2*c*d1Gc+d2Gc*d1A)+d1Gc*(2*m*d1T*d1A + d2T*(d1A^2)+d1T*d2A + m*(m-1)-d2Gm*(d1A^2)-d1Gm*d2A))/(2*(c-d1Gc*d1A));
    # temp = temp - d1Gc*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)*(m+d1T*d1A-d1Gm*d1A)/(2*(c-d1Gc*d1A)^2);
    temp = temp*constants[j]*(1-alpha[j]);
    
    ESystCont=ESystCont + temp;
    #print(temp)
  }
  if(th>c)
  {
    for(j in c:(th-1))
    {
      m=j;
      d1Gm=ServiceGFDerivative(1,m);
      d2Gm=ServiceGF2d(1,m);
      temp=0;
      temp = c*d1T*(m+d1Gc*d1A)/(c-d1Gc*d1A);
      temp = temp + (d1T*c*(c-1)+c*d2T*d1A)/(2*(c-d1Gc*d1A));
      #temp1 = temp1 + (c*d1T*d2A/(d1A*(c-d1Gc*d1A)));
      temp = temp - c*d1T*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)/(2*((c-d1Gc*d1A)^2)); 
      
      #temp = temp + d2T*d1A/2;
      #temp = temp + (d1T*d2Gc*(d1A^2)+d1Gc*(2*c*d1T*d1A+d2T*(d1A^2)+d1T*d2A))/(2*(c-d1Gc*d1A));
      #temp = temp - (d1T*d1A*d1Gc*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A))/(2*(c-d1Gc*d1A)^2);
      
      temp=temp*constants[j]*(1-alpha[j]);
      #print(temp);
      ESystCont = ESystCont + temp;
      #print(temp)
    }
  }
  #print("after d_0(m)")
  #print(ESystCont)
  #for d(0)
  temp1=0;
  temp1 = (2*c*d1T*d1Gc*d1A+d1T*c*(c-1)+c*d2T*d1A)/(2*(c-d1Gc*d1A));
  #temp1 = temp1 + (c*d1T*d2A/(d1A*(c-d1Gc*d1A)));
  temp1 = temp1 - c*d1T*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)/(2*((c-d1Gc*d1A)^2)); 
  temp1 = temp1*constants[th];
  ESystCont = ESystCont + temp1;
  #print("after d(0)")
  #print(ESystCont)
  for(j in (th+1):(c+th-1))
  {
    m = j-th;
    d1Gm = ServiceGFDerivative(1,m);
    d2Gm = ServiceGF2d(1,m);
    temp2=0;
    
    temp2 = m*(d1Gm - d1Gc);
    temp2 = temp2 + d1A*d2Gm/2;
    temp2 = temp2 + d2Gc*d1A*(d1Gm*d1A-m)/(2*(c-d1Gc*d1A));
    temp2 = temp2 + d1Gc*(2*c*d1Gm*d1A+d2Gm*d1A^2+d1Gm*d2A-m*(m-1)-2*m*d1Gc*d1A)/(2*(c-d1Gc*d1A));
    temp2 = temp2 - d1Gc*(d1Gm*d1A-m)*(c*(c-1)-d2Gc*d1A^2 - d1Gc*d2A)/(2*(c-d1Gc*d1A)^2);
    
    # temp2 = temp2 + (d2Gc*d1A*(d1Gm*d1A-m)+d1Gc*(2*c*d1Gm*d1A+d2Gm*(d1A^2)+d1Gm*d2A-m*(m-1)-2*m*d1Gc*d1A))/(2*(c-d1Gc*d1A));
    # temp2 = temp2 + (d1Gc*(d1Gm*d1A-m)*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A))/(2*(c-d1Gc*d1A)^2);
    # temp2 = temp2 + (c-m+d1Gm*d1A-d1Gc*d1A)*d2Gc*d1A/(2*(c-d1Gc*d1A));
    
    # temp2 = temp2 + (d2Gm-d2Gc)*d1A/2;
    # temp2 = temp2 + d2A*(d1Gm-d1Gc)/(2*(c-d1Gc*d1A)*d1A);
    # temp2 = temp2 + (d1Gm-d1Gc)*(c*(c-1)-d2Gm*(d1A^2)-d1Gm*d2A)/(2*(c-d1A*d1Gc));
    # temp2 = temp2 + d1Gc*(c*(c-1)+2*c*d1Gm*d1A + d2Gm*(d1A^2)+d1Gm*d2A-m*(m-1)-2*m*d1Gc*d1A - d2Gc*(d1A^2)-d1Gc*d2A)/(2*(c-d1Gc*d1A));
    # temp2 = temp2 + d1Gc*(c+d1Gm*d1A-m-d1Gc*d1A)*(d2Gc*(d1A^2)+d1Gc*d2A-c*(c-1))/(2*(c-d1Gc*d1A)^2);
    
    temp2 = temp2*constants[j];
    ESystCont = ESystCont + temp2;
    #print(temp2)
  }
  #print('final')
  #print(ESystCont)
  ESystCont=ESystCont;
}
## This is a temporary function for analysis. most of the code is just copy of Esystemcont
## returns the multipliers to be used for constants in esystcont
TempMultiplierconsts = function(c,lambda,alpha,constants,th)
{  
  TempMultiplierconsts = {}
  
  d1A=ArrivalGFDerivative(1,lambda);
  d1Gc = ServiceGFDerivative(1,c);
  d2Gc = ServiceGF2d(1,c);
  d1T= VacationGFDerivative(1);
  d2T= vacationGF2d(1);
  d2A = ArrivalGF2d(1,lambda);
  
  for(j in 1:(c-1))
  {
    m=j;
    d1Gm = ServiceGFDerivative(1,m);
    d2Gm = ServiceGF2d(1,m);
    temp=0;
    temp = (d2T-d2Gm)*d1A/2;
    temp = temp + m*(d1T-d1Gm);
    temp = temp + c*d1Gc*(m+d1T*d1A-d1Gm*d1A)/(c-d1Gc*d1A);
    temp = temp + d2Gc*d1A*(m+d1T*d1A-d1Gm*d1A)/(2*(c-d1Gc*d1A));
    temp = temp + d1Gc*(2*m*d1T*d1A+d2T*(d1A^2)+d1T*d2A + m*(m-1)-d2Gm*(d1A^2)-d1Gm*d2A)/(2*(c-d1Gc*d1A));
    temp = temp - (c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)*(m+d1T*d1A-d1Gm*d1A)/(2*(c-d1Gc*d1A)^2);
    TempMultiplierconsts = c(TempMultiplierconsts,temp)
    temp = temp*constants[j]*(1-alpha[j]);
    
  }
  if(th>c)
  {
    for(j in c:(th-1))
    {
      m=j;
      d1Gm=ServiceGFDerivative(1,m);
      d2Gm=ServiceGF2d(1,m);
      temp=0;
      temp = c*d1T*(m+d1Gc*d1A)/(c-d1Gc*d1A);
      temp = temp + (d1T*c*(c-1)+c*d2T*d1A)/(2*(c-d1Gc*d1A));
      temp = temp - c*d1T*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)/(2*((c-d1Gc*d1A)^2)); 
      TempMultiplierconsts = c(TempMultiplierconsts,temp)
      temp=temp*constants[j]*(1-alpha[j]);
    }
  }
  #for d(0)
  temp1=0;
  temp1 = (2*c*d1T*d1Gc*d1A+d1T*c*(c-1)+c*d2T*d1A)/(2*(c-d1Gc*d1A));
  #temp1 = temp1 + (c*d1T*d2A/(d1A*(c-d1Gc*d1A)));
  temp1 = temp1 - c*d1T*(c*(c-1)-d2Gc*(d1A^2)-d1Gc*d2A)/(2*((c-d1Gc*d1A)^2)); 
  TempMultiplierconsts = c(TempMultiplierconsts,temp1)
  temp1 = temp1*constants[th];
  for(j in (th+1):(c+th-1))
  {
    m = j-th;
    d1Gm = ServiceGFDerivative(1,m);
    d2Gm = ServiceGF2d(1,m);
    temp2=0;
    
    temp2 = m*(d1Gm - d1Gc);
    temp2 = temp2 + d1A*d2Gm/2;
    temp2 = temp2 + d2Gc*d1A*(d1Gm*d1A-m)/(2*(c-d1Gc*d1A));
    temp2 = temp2 + d1Gc*(2*c*d1Gm*d1A+d2Gm*d1A^2+d1Gm*d2A-m*(m-1)-2*m*d1Gc*d1A)/(2*(c-d1Gc*d1A));
    temp2 = temp2 - d1Gc*(d1Gm*d1A-m)*(c*(c-1)-d2Gc*d1A^2 - d1Gc*d2A)/(2*(c-d1Gc*d1A)^2);
    TempMultiplierconsts = c(TempMultiplierconsts,temp2)
    temp2 = temp2*constants[j];
  }
  TempMultiplierconsts
}
####################################################################
# For generating function V(z,y,x) it returns the value of V(1,1,0)
####################################################################
v110=function(c,lambda,alpha,constants,th)
{
  v110=1;
  d1A=ArrivalGFDerivative(1,lambda);
  d1Gc = ServiceGFDerivative(1,c);
  d2Gc = ServiceGF2d(1,c);
  d1T= VacationGFDerivative(1);
  d2T= vacationGF2d(1);
  d2A = ArrivalGF2d(1,lambda);
  # for( m in 1:(c-1))
  # {
  #   d1Gm = ServiceGFDerivative(1,m);
  #   #for d_0(m)
  #   v110 = v110 + (m+d1T*d1A-d1Gm*d1A)*(1-alpha[m])*constants[m]/(c-d1Gc*d1A); 
  # }
  # for( m in c:(th-1))
  # {
  #   v110 = v110 + (c+d1T*d1A-d1Gc*d1A)*(1-alpha[m])*constants[m]/(c-d1Gc*d1A);
  # }
  # #for d(0)
  # v110 = v110 + (c+d1T*d1A-d1Gc*d1A)*constants[th]/(c-d1Gc*d1A);
  # for( j in (th+1):(c+th-1))
  # {
  #   #for d(m), m>0
  #   m=j-th;
  #   d1Gm = ServiceGFDerivative(1,m);
  #   v110 = v110 + (c-m+d1Gm*d1A-d1Gc*d1A)*constants[j]/(c-d1Gc*d1A);
  # }
  for( m in 1:(c-1))
  {
    d1Gm = ServiceGFDerivative(1,m);
    v110 = v110 - (d1T-d1Gm)*(1-alpha[m])*constants[m];
  }
  if(th>c)
  {
    for( m in c:(th-1))
    {
      v110 = v110 - (d1T - d1Gc)*(1-alpha[m])*constants[m];
    }
  }
  v110 = v110 - (d1T-d1Gc)*constants[th];
  for( j in (th+1):(c+th-1))
  {
    m = j-th;
    d1Gm = ServiceGFDerivative(1,m);
    v110 = v110 - (d1Gm-d1Gc)*constants[j];
  }
  v110=v110/d1Gc;
}

###################################################################
# Expected server content using the generating function
###################################################################
EServerCont=function(c,lambda,alpha,constants,th)
{
  # v110=v110(c,lambda,alpha,constants,th);
  # print("new");
  # print(v110);
  # d1Gc= ServiceGFDerivative(1,c);
  # EServerCont = c*d1Gc*(v110-constants[th]);
  # print(EServerCont);
  # for(j in 1:(c-1))
  # {
  #   m=j;
  #   d1Gm = ServiceGFDerivative(1,m);
  #   EServerCont=EServerCont-m*d1Gm*(1-alpha[j])*constants[j];
  # }
  # print(EServerCont);
  # for(j in c:(th-1))
  # {
  #   EServerCont = EServerCont -c*d1Gc*(1-alpha[j])*constants[j];
  # }
  # print(EServerCont);
  # for(j in (th+1):(c+th-1))
  # {
  #   m=j-th;
  #   d1Gm= ServiceGFDerivative(1,m);
  #   EServerCont=EServerCont+ (m*d1Gm-c*d1Gc)*constants[j];
  # }
  # print(EServerCont);
  d1T = VacationGFDerivative(1);
  d1Gc = ServiceGFDerivative(1,c);
  EServerCont = c;
  
  for( j in 1:(c-1))
  {
    m=j;
    d1Gm = ServiceGFDerivative(1,m);
    EServerCont = EServerCont - (c*d1T-c*d1Gm+m*d1Gm )*constants[j]*(1-alpha[j]);
  }
  if(th>c)
  {  
    for( j in c:(th-1))
    {
      EServerCont = EServerCont -c*d1T*(1-alpha[j])*constants[j];
    }
  }
  EServerCont = EServerCont - c*d1T*constants[th];
  for( j in (th+1):(c+th-1))
  {
    m=j-th;
    d1Gm = ServiceGFDerivative(1,m);
    EServerCont = EServerCont +(m-c)*d1Gm*constants[j];
  }
  EServerCont;
}

####################################################################
# for generating function V(z,y,x) returns the value V(1,0,1)
####################################################################
V101 = function(th,alpha,constants)
{
  V101=constants[th];
  for(j in 1:(th-1))
  {
    V101=V101+(1-alpha[j])*constants[j];
  }
  V101=V101*(VacationGFDerivative(1));
}

####################################################################
# To find the number of customers served in a random batch
####################################################################

EServerContinBatch=function(c,lambda,alpha,constants,th)
{
  EServerContinBatch=EServerCont(c,lambda,alpha,constants,th);
  V101=V101(th,alpha,constants);
  EServerContinBatch=EServerContinBatch/(1-V101);
}
####################################################################
# Returns Expected Queue Content when server is in vacation
###################################################################
EQConVacation = function(th,lambda,alpha,constants)
{
  d1A=ArrivalGFDerivative(1,lambda);
  d1T=VacationGFDerivative(1);
  d2T=vacationGF2d(1);
  EQConVacation=d2T*d1A*constants[th]/2;
  for(j in 1:(th-1))
  {
    EQConVacation= EQConVacation + (j*d1T+d2T*d1A/2)*(1-alpha[j])*constants[j];
  }
  EQConVacation=EQConVacation/V101(th,alpha,constants);
}

####################################################################
# Returns probability of the server being busy
###################################################################

ProbServerBusy=function(th,alpha,constants)
{
  ProbServerBusy=1-Re(V101(th,alpha,constants));
}

###################################################################
# V100
###################################################################
v100=function(th,alpha,constants)
{
  v100=constants[th];
  for(j in 1:(th-1))
  {
    v100 = v100+ (1-alpha[j])*constants[j];
  }
  v100;
}

####################################################################
# New service initiated
# gives us P(SK+1=1, Sk=0)
####################################################################

StartNewService=function(th,lambda,alpha,constants)
{
  StartNewService=1;
  v100=v100(th,alpha,constants);
  for(j in 1:(th-1))  
  {
    StartNewService = StartNewService+ (alpha[j]-1)*constants[j]/v100;
  }
  StartNewService = StartNewService - VacationGF(ArrivalGF(0,lambda))*constants[th]/v100;
  StartNewService = Re(StartNewService);
}

####################################################################
# Recovery time objective (RT0)
# this gives the upper bound on the rpo computed in the paper using
# recursion of the expected values.
####################################################################

rpoUB=function(th,lambda,alpha)
{
  d1T = VacationGFDerivative(1);
  A0 = ArrivalGF(0,lambda);
  rpoUB = NumVacationFP(th,lambda,alpha)*d1T;
  rpoUB = rpoUB + d1T/(1- VacationGF(A0))-(1/(1-A0));
  rpoUB;
}

######################################################################
# NumVacationFP
# Number of vacations seen by the first packets which arrives in an empty system
######################################################################

NumVacationFP = function(th,lambda,alpha)
{
  M = matrix(0,nrow=th,ncol=th);
  tAvec = {};
  for(i in 1:th)
  {
    tAvec[i]=Tai(i-1,lambda); #tAvec[i] =ta(i-1) because we dont have 0 subscript in R
  }
  for(j in 2:th)
  {
    M[1,j] = (1-alpha[j-1])*tAvec[j]/(1-tAvec[1]);
  }
  for( i in 2:th)
  {
    for(j in i:th)
    {
      M[i,j]=(1-alpha[j-1])*tAvec[j-i+1]; #corresponds to tA(j-i);
    }
  }
  beta = matrix(0,nrow=1,ncol=th);
  beta[1,1]=1;
  M0 = 1-rowSums(M);
  invM = solve(diag(th)-M);
  NumVacationFP = beta%*%invM;
  NumVacationFP = NumVacationFP%*%invM;
  NumVacationFP = NumVacationFP%*%M0;
  NumVacationFP = NumVacationFP-1;
}

#####################################
## second method of computing expected Age nu
#######################################

expected_Agenu2 = function( th, lambda, alpha)
{
  rivec = {};
  
  d1T = VacationGFDerivative(1);
  A0 = ArrivalGF(0,lambda);
  for(i in 1:th)
  {
    rivec[i]=Tai(i-1,lambda)/(1- VacationGF(A0)); #tAvec[i] =ta(i-1) because we dont have 0 subscript in R
  }
  
  expected_Agenu2 = d1T/(1- VacationGF(A0))- 1/(1-A0);
  gip = gip(th,lambda,alpha)
  for( i in 1:(th-1))
  {
    expected_Agenu2 = expected_Agenu2 + (1-alpha[i])*rivec[i+1]*gip[i];
  }
  # exp_1 = rpoUB(th,lambda,alpha)
  # if( abs( exp_1 - expected_Agenu2) > 0.0001 )
  # {
  #   print("Problem")
  #   print(exp_1[1,1])
  #   print(expected_Agenu2)
  # }
  expected_Agenu2
}
#######################################################################
# return a vector of first moment of gip 
#######################################################################
gip = function(th, lambda, alpha)
{
  gip = rep_len(0,th-1);
  tivec = {}
  
  A0 = ArrivalGF(0,lambda);
  for(i in 1:th)
  {
    tivec[i] = Tai(i-1,lambda);#because we dont have 0 subscript in R
  }
  
  for( i in (th-1):1)
  {
    gip[i] = VacationGFDerivative(1);
    if(i< (th-1) )
    {  
      for( j in (i+1):(th-1))
      {  
        gip[i] = gip[i] + (1-alpha[j])*gip[j]*tivec[j-i+1];
      }
    }
    gip[i] = gip[i]/(1-(1-alpha[i])*tivec[1])
  }
  gip
}
#######################################################################
# returns second moments of gips 
#######################################################################
gip2 = function(th, lambda, alpha)
{
  
  tivec = {};
  tipvec = {};
  
  
  d1T = VacationGFDerivative(1);
  d2T = vacationGF2d(1)
  A0 = ArrivalGF(0,lambda);
  for(i in 1:th)
  {
    #ta(i-1) because we dont have 0 subscript in R
    tivec[i] = Tai(i-1,lambda);
    tipvec[i]= Aip(i-1,lambda);
  }
  
  gip = gip(th,lambda,alpha);
  gip2 = rep_len(0,th-1);
  for( i in (th-1):1)
  {
    gip2[i] = d2T + 2*gip[i]*(1-alpha[i])*tipvec[1];
    if(i<(th-1))
    {
      for( j in (i+1):(th-1))
      {
        gip2[i] = gip2[i] + (1-alpha[j])*(gip2[j]*tivec[j-i+1]+2*gip[j]*tipvec[j-i+1]);
      }
    }
    gip2[i] = gip2[i]/(1-(1-alpha[i])*tivec[1]);
  }
  gip2
}
########################################################################
# sdv agenu
# exact computation of sdv of age nu
########################################################################

sdv_age_nu = function(th, lambda, alpha)
{
  ### NOTE THAT THIS IS VALID ONLY FOR POISSON ARRIVAL
  p = prob_poisson();
  if( p < 1)
  {
    stop("can not do this calculation for non-Poisson arrival distribution")
  }
  tivec = {};
  tip = {};
  rivec = {};
  rip = {};
  
  d1T = VacationGFDerivative(1);
  A0 = ArrivalGF(0,lambda);
  for(i in 1:th)
  {
    tivec[i] = Tai(i-1,lambda); #tAvec[i] =ta(i-1) because we dont have 0 subscript in R
    tip[i]   = Tip(i-1,lambda);
    rivec[i] = tivec[i]/(1- VacationGF(A0));
    rip[i]   = Rip(i-1,lambda);
  }
  rip[1]=0;
  print(tivec)
  print(tip)
  print(rivec)
  print(rip)
  
  gip = gip(th,lambda,alpha);
  gip2 = gip2(th,lambda,alpha)  

  sdv_age_nu = (vacationGF2d(1)/(1- VacationGF(A0))) - (2*VacationGFDerivative(1)/((1- VacationGF(A0))*(1-A0))) + (2/((1-A0)^2));
  
  for( i in 1:(th-1))
  {
    sdv_age_nu = sdv_age_nu + (1-alpha[i])*(2*rip[i+1]*gip[i] + rivec[i+1]*gip2[i]);
  }
  eagenu = expected_Agenu2(th,lambda,alpha);
  print(eagenu)
  sdv_age_nu = sdv_age_nu + eagenu - eagenu^2; ##this is variance
  #print(sdv_age_nu)
  sdv_age_nu = sqrt(sdv_age_nu);
  sdv_age_nu;
}
########################################################################
# VarN
# Variance of number of vacations seen by first packet in an empty system
########################################################################
VarN = function(th,lambda,alpha)
{
  M = matrix(0,nrow=th,ncol=th);
  tAvec = {};
  for(i in 1:th)
  {
    tAvec[i]=Tai(i-1,lambda); #tAvec[i] =ta(i-1) because we dont have 0 subscript in R
  }
  for(j in 2:th)
  {
    M[1,j] = (1-alpha[j-1])*tAvec[j]/(1-tAvec[1]);
  }
  for( i in 2:th)
  {
    for(j in i:th)
    {
      M[i,j]=(1-alpha[j-1])*tAvec[j-i+1]; #corresponds to tA(j-i);
    }
  }
  beta = matrix(0,nrow=1,ncol=th);
  beta[1,1]=1;
  M0 = 1-rowSums(M);
  invM = solve(diag(th)-M);
  NumVacation = beta%*%invM;
  NumVacation = NumVacation%*%invM;
  NumVacation = NumVacation%*%M0;
  NumVacation = NumVacation - 1;
  N21 = beta%*%invM;
  N21 = N21%*%invM;
  N21 = N21%*%M;
  N21 = N21%*%invM;
  N21 = N21%*%M0;
  N21 = 2*N21;
  VarN = N21 - NumVacation - NumVacation^2;
}
############################################################################
# save the graph with the parameters and filename!
############################################################################
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
  filename = q=paste("C:/Users/sapoorv/Downloads/latex/Data Backups/figs/plot_", variable_name,"_load_",type,".pdf", sep = "");
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

########################################################################
# VarR
# Variance of the remaining vacation time seen by the first packet entering in empty system
########################################################################
VarR = function(th,lambda,alpha)
{
  A0 = ArrivalGF(0,lambda);
  d1T = VacationGFDerivative(1);
  d2T = vacationGF2d(1);
  ET2 = d2T+d1T;
  ER2 = ET2/(1- VacationGF(A0));
  ER2 = ER2 - (2*d1T+1)/((1- VacationGF(A0))*(1-A0));
  ER2 = ER2 + (2-(A0+1)*VacationGF(A0))/((1- VacationGF(A0))*(1-A0)^2);
  ER1 = d1T/(1- VacationGF(A0))-(1/(1-A0));
  VarR = ER2-ER1^2;
}

########################################################################
# VarSp
# Variance of S' as defined in the paper
#######################################################################
VarSp = function(th,lambda, alpha)
{
  EN = NumVacationFP(th,lambda,alpha);
  d1T = VacationGFDerivative(1);
  d2T = vacationGF2d(1);
  VarT = d2T + d1T - d1T^2;
  VarSp = EN*VarT + VarR(th,lambda,alpha);
}

#######################################################################
# Varrpob
# this computes the upper and lower bounds on the variance of rpo
#######################################################################
Varrpob = function(th,lambda,alpha)
{
  d1T = VacationGFDerivative(1);
  Varrpob = {};
  VarSp = VarSp(th,lambda,alpha);
  VarN = VarN(th,lambda,alpha);
  Varrpob[1] = (sqrt(VarSp)-d1T*sqrt(VarN))^2;
  Varrpob[2] = (sqrt(VarSp)+d1T*sqrt(VarN))^2;
  Varrpob[3] = max(sqrt(VarSp), d1T*sqrt(VarN));
  Varrpob;
}

########################################################################
# approximation of derivative of T(zA(y))
########################################################################
deriv_approximation<- function(m,lambda,z_root)
{
  deriv_approximation = 0;
  sign = 1;
  h = 10^-5;
  for( r in 0:m)
  {
    deriv_approximation = deriv_approximation + sign*choose(m,r)*VacationGF(z_root+(m-r)*h)
    sign = sign*(-1);
  }
  deriv_approximation = deriv_approximation/(h^m)
}
#######################################################################
# age_tail
# this compute the log of tail of Pr(age > 30 \times E(T)) i.e. log10(pr(age)>30E(T))
#######################################################################
Age_tail <- function(th,lambda,alpha)
{
  min_abs = 100;
  
  A0 = ArrivalGF(0,lambda);
  alpha_smallest = min(alpha)
  z_root = -1;
  while(z_root < 0)
  {  
    result1 = BBsolve(par=1.9*abs(z_root), fn=myFunction_sing,control = list(maxit = 5000, tol = 10^-15), quiet = TRUE, params=c(alpha_smallest,lambda))
    z_root = result1$par; 
  }  
  # z_root is the dominant singularity!
  t_kzd = {};
  d1A0 = ArrivalGFDerivative(0,lambda)
  
  r_1zd = d1A0*(VacationGF(z_root*A0)-VacationGF(A0))/((z_root-1)*A0*(1- VacationGF(A0)));
  
  for(m in 1:(th-1))
  {
    new_val <- deriv_approximation(m,lambda,z_root);
    t_kzd <- c(t_kzd, new_val);
  }
  Tzd = VacationGF(z_root);
  Gk = rep_len(0,length.out = th-1);
  Gk[th-1] = (Tzd-(1-alpha[th-1])*VacationGF(z_root*A0))/(1-(1-alpha[th-1])*VacationGF(z_root*A0));
  for( q in 2:(th-2))
  {
    q2 = th-q; 
    val = (Tzd-(1-alpha[q2])*VacationGF(z_root*A0)); 
    for( r in (q2+1):(th-1))
    {
      val = val + (1-alpha[r])*(Gk[r]-1)*t_kzd[r-q2];
    }
    val = val /(1-(1-alpha[q2])*VacationGF(z_root*A0));
    Gk[q2] = val;
  }
  ord = sum(new_alpha == alpha_smallest);
  n = 30*VacationGFDerivative(1);
  b_k = (Tzd-(1-alpha[1])*VacationGF(z_root*A0))
  for( q in (1+ord):(th-1))
  {
    b_k = b_k + (1-alpha[q])*(Gk[q]-1)*t_kzd[q-ord];
  }
  b_k = b_k*r_1zd*t_kzd[1]^(ord-1)/(A0*VacationGFDerivative(z_root*A0))^ord; 
  
  y_n = b_k/z_root^n;
  Age_tail = log10(y_n)
}
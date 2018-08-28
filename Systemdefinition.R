gammaparam <- function()
{
  gammaparam = 3.5 # has to more than 3 for finite mean and variance
}

starting_powerlaw <- function()
{
  starting_powerlaw = 25
}
prob_poisson <- function()
{
  prob_poisson = 1 #0.995
}
#############################################
# Return A(z) for a give value of z
#############################################
ArrivalGF <- function(z,lambda)
{
  
  #Poisson Arrival
  p = prob_poisson();
  ArrivalGF = p*exp(lambda*(z-1)); 
  
  # poisson + power law distribution  
  
  #lambda is the probability of non zero arrivals
  gamma = gammaparam()
  ArrivalGF = ArrivalGF + (1-p)*Liz(z,gamma)/Liz(1,gamma)
  
}

############################################
# Returns Li_gamma (z) the zeta function for power law distribution
############################################
Liz <- function(z,gamma)
{
  Liz = 0
  a = starting_powerlaw()
  for( i in (a:200)) #should go to infinity rather than 100 but 100 should be enought as abs(z) < 1
  {
    Liz = Liz + z^i/(i^gamma)
  }
  Liz
}
#############################################
# Lizprime derivative of Liz
#############################################
Lizprime <- function(z,gamma)
{
  Lizprime = 0
  a = starting_powerlaw()
  for( i in (a:200)) #should go to infinity rather than 100 but 100 should be enought as abs(z) < 1
  {
    Lizprime = Lizprime + i*z^(i-1)/(i^gamma)
  }
  Lizprime
}
#############################################
# Liz2prime second derivative of Liz
#############################################
Liz2prime <- function(z,gamma)
{
  Liz2prime = 0
  a = starting_powerlaw();
  for( i in (a:200)) #should go to infinity rather than 100 but 100 should be enought as abs(z) < 1
  {
    Liz2prime = Liz2prime + i*(i-1)*z^(i-2)/(i^gamma)
  }
  Liz2prime
}
#############################################
# Return first derivative of arrival gf
#############################################
ArrivalGFDerivative <- function(z, lambda)
{
  #Poisson Arrival
  p = prob_poisson()
  gamma = gammaparam();
  ArrivalGFDerivative = p*lambda*exp(lambda*(z-1))+(1-p)*Lizprime(z,gamma)/Liz(1,gamma)
}
#############################################
# Return second derivative of arrival gf
#############################################

ArrivalGF2d <- function(z,lambda)
{
  #Poisson
  gamma = gammaparam();
  p = prob_poisson();
  ArrivalGF2d=p*lambda^2*exp(lambda*(z-1))+(1-p)*Liz2prime(z,gamma)/Liz(1,gamma)
}  

#############################################
# Return factorial of argument n
#############################################
factorial = function(n)
{
  if(n==0||n==1)
    factorial= 1
  else
    factorial = (n*factorial(n-1));
}

#########################################################
#Returns probability of n arrivals in t slots for arrival rate lambda
#########################################################
ArrivalDist = function(n,t,lambda)
{
  #Poisson
  #ArrivalDist = (exp(-lambda*t)*((lambda*t)^n))/(factorial(n)); #poisson
  a = starting_powerlaw();
  p = prob_poisson();
  ArrivalDist = (p^t)*(exp(-lambda*t)*((lambda*t)^n))/(factorial(n));
  gamma = gammaparam();
  if(n >=a && p < 1)
  {
    for( i in a:n)
    {
      ArrivalDist = ArrivalDist + t*(p^(t-1))*(1-p)*(exp(-lambda*(t-1)))*(lambda*(t-1))^(n-i)/(factorial(n-i)*Liz(1,gamma)*i^gamma)
    }
  }
  ArrivalDist
}

############################################################
# Returns the probability of n arrivals in a vacation period
############################################################
Tai=function(n,lambda)
{
  
  vacation_probs = vacation_probabilities();
  Tai=0;
  for(j in 1:length(vacation_probs))
  {
    if(vacation_probs[j]>0)
      Tai = Tai + ArrivalDist(n,j-1,lambda)*vacation_probs[j];
  }
  Tai;
}
###############################################################
# expected length of the vacation which sees n arrivals
##############################################################
Aip = function(n,lambda)
{
  # THIS IS VALID ONLY FOR POISSION ARRIVALS, change the ArrivalGF appropriately before calling this
  p = prob_poisson();
  if( p < 1)
  {
    stop("can not do this calculation for non-Poisson arrival distribution")
  }
  vacation_probs = vacation_probabilities();
  Aip = 0;
  for( j in 1:length(vacation_probs))
  {
    if(vacation_probs[j]>0)
      Aip = Aip + ArrivalDist(n,j-1,lambda)*(j-1)*vacation_probs[j];
  }
  Aip
}
  
###############################################################
# expected length of the vacation which sees n arrivals
##############################################################
Tip = function(n,lambda)
{
  # THIS IS VALID ONLY FOR POISSION ARRIVALS, change the ArrivalGF appropriately before calling this
  p = prob_poisson();
  if( p < 1)
  {
    stop("can not do this calculation for non-Poisson arrival distribution")
  }
  Tip = Aip(n,lambda)
}
###############################################################
# expected length of the vacation in which nu arrives and sees n arrivals
##############################################################
Rip = function(n,lambda)
{
  p = prob_poisson();
  if( p < 1)
  {
    stop("can not do this calculation for non-Poisson arrival distribution")
  }
  # THIS IS VALID ONLY FOR POISSION ARRIVALS, change the ArrivalGF appropriately before calling this
  vacation_probs = vacation_probabilities();
  Rip = Tip(n,lambda);
  A0 = ArrivalGF(0,lambda);
  powersums = rep_len(0,length(vacation_probs));
  powersums[1] = 1;
  for( j in 2:length(powersums))
  {
    powersums[j] = powersums[j-1] + j^n;
  }
  
  for( j in 1:length(vacation_probs))
  {
    if( vacation_probs[j]>0)
    {
      t =j-1
      Rip = Rip - exp(-lambda*t)*(lambda^n)*vacation_probs[j]*powersums[t]/factorial(n);  
    }
  }
  Rip = Rip/(1 - VacationGF(A0));
 
  Rip
}
#################################################################
## Rvacat(z)
## the generating function of the remaining vacation time in which nu arrives
################################################################
Rvacat = function(z,lambda)
{
  A0 = ArrivalGF(0,lambda);
  Rvacat = (1-A0)*(VacationGF(A0)-VacationGF(z))/((1- VacationGF(A0))*(A0-z));
}
#############################################################
# degraded : seems to have not been used anywhere in the code.
##############################################################
Lai=function(n,lambda)
{
  stop('do not execute this function');
  vacation_probs = vacation_probabilities();
  Lai = 0;
  for(j in 1:length(vacation_probs))
  {
    if(vacation_probs[j]>0)
      Lai = Lai + ArrivalDist(n,j,lambda)*exp(lambda)*(j-1)*vacation_probs[j]; #length of vacation = j-1
  }
  Lai;
}

###########################################################
# prob of success of Geo(p)
##############################################################
psuccess=function(m)
{
  constant = 1;
  psuccess = 1/(m+constant);
}  
############################################################
#Returns the service distribution of m packets
###########################################################
EServicet = function(m)
{
  EServicet=ServiceGFDerivative(1,m);
}

############################################################
#Generating function of the service time of m packets
############################################################
ServiceGF=function(z,m)
{
  # IF YOU CHANGE THIS MAKE SURE YOU CHANGE THE SOLVER FUNCTION TO GET ROOTS OF z^c - G_c(A(z))
  am = m  #bm = m+2 so that avg = m+1
  
  #ServiceGF = z^(m) 
  ServiceGF = (z^am + z^(am+1) + z^(am+2))/3
  #U(a,b)
  # p=psuccess(m);
  # ServiceGF = (p*z/(1-((1-p)*z)));
}

#############################################################
#Derivative of function above
#############################################################
ServiceGFDerivative=function(z,m)
{
  # am = m
  # bm = m+2 #so that avg is m+1
   ServiceGFDerivative = ((m)*z^(m-1) + (m+1)*z^m + (m+2)*z^(m+1))/3
  # p=psuccess(m);
  # ServiceGFDerivative = p/((1-(1-p)*z)^2);
}

ServiceGF2d=function(z,m)
{
  am = m
  ServiceGF2d = 0
  #if(m>1)
    ServiceGF2d = (m*(m+1)*z^(m-1) + (m+2)*(m+1)*z^m )/3
  # ServiceGF2d
  #  
  if( am >=2) {
    ServiceGF2d = ServiceGF2d + (am*(am-1)*z^(am-2))/3
  }
  ServiceGF2d
  # p=psuccess(m);
  # ServiceGF2d=2*p*((1-p))/((1-(1-p)*z)^3);
}
###########################################################################
# Vacation[i] --> probability that the vacation is of length i
#                 ith element quantifies the vacation prob is of length i-1
###########################################################################
vacation_probabilities=function()
{
  #vacation_probabilities=c(0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2);  #U[8,12]
  len = 31;
  vacation_probabilities = rep_len(0,len)
  vacation_probabilities[len-2]= 0.2
  vacation_probabilities[len-1]= 0.2
  vacation_probabilities[len]= 0.2
  vacation_probabilities[len+1]= 0.2
  vacation_probabilities[len+2]= 0.2
  vacation_probabilities
  #vacation_probabilities=integer(10);
  #vacation_probabilities[10]=1;
  #vacation_probabilities=c(0,1);
}  

###############################################################
# Generating funciton of the length of the vacation
###############################################################
VacationGF=function(z)
{
  vacation_prob = vacation_probabilities();
  VacationGF= 0;
  for(j in 1:length(vacation_prob))
  {
    VacationGF = VacationGF + vacation_prob[j]*z^(j-1);
  }  
  VacationGF;
}  

###############################################################
# Derivative of the function above
###############################################################
VacationGFDerivative=function(z)
{
  vacation_prob = vacation_probabilities();
  VacationGFDerivative= 0;
  if( length(vacation_prob) > 1)
  {
    for(j in 2:length(vacation_prob))
    {
      VacationGFDerivative = VacationGFDerivative + vacation_prob[j]*(j-1)*z^(j-2);
    }  
  }
  VacationGFDerivative;
}

vacationGF2d = function(z)
{
  vacation_prob = vacation_probabilities();
  VacationGF2d=0;
  if( length(vacation_prob) > 2)
  {
    for(j in 3:length(vacation_prob))
    {
      VacationGF2d = VacationGF2d + vacation_prob[j]*(j-1)*(j-2)*z^(j-3);
    }   
  }
  VacationGF2d;
}
B= function(lambda,service_parameters,n,z)
{
  y = ArrivalGF(z,lambda);
  B = ServiceGF(y,n);
}
dB = function(lambda,service_parameters,n,z)
{
  y = ArrivalGF(z,lambda);
  dB = ServiceGFDerivative(y,n)*ArrivalGFDerivative(z,lambda);
}
source("Systemdefinition.R");
source("performance measures.R");
c=10;
th =30; 

len=20;
poissonarrivalrates = seq(0.2, 0.70, length.out = len);
load = {};
sdrpo1 = {};
sdrpo2 = {};
sdrpo3 = {};
alpha1 = {}; 
alpha2 = {}; 
alpha3 = {}; 
for(i in 1:(th-1))
{
  alpha1[i]= i/(3*c);
  alpha2[i]= 0;
  alpha3[i]= 1;
}

p = prob_poisson();
if( p < 1)
{
  stop("can not do this calculation for non-Poisson arrival distribution")
}

for(lambda in poissonarrivalrates)
{
  load =  c(load,ArrivalGFDerivative(1,lambda)*ServiceGFDerivative(1,c)/c);
  sdrpo1 = c(sdrpo1,sdv_age_nu(th,lambda,alpha1));
  sdrpo2 = c(sdrpo2,sdv_age_nu(th,lambda,alpha2));
  sdrpo3 = c(sdrpo3,sdv_age_nu(th,lambda,alpha3));
}
#alpha[i] = 0.5+i*i/(2*th*th);

system('taskkill /f /im AcroRd32.exe');
saveGraph("sd_rpo", "",sdrpo1, sdrpo2,sdrpo3,load, 'SD[Age_nu]', 'load', 0, 0.5, "topright")
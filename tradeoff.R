#install.packages(c("alabama","nloptr"), lib = getwd(),repos = "http://cran.us.r-project.org")
require(alabama)
require(nloptr)
#library(alabama)
#library(nloptr)
rm(list = ls(all.names = TRUE))
source("Systemdefinition.R");
source('Solverfunctions.R');
source('matrixcoefficients.R');
source('performance measures.R');
system('taskkill /f /im AcroRd32.exe');

c=10;

max_th = 60;

lambda = 0.6
alpha = {}
for( i in 1:max_th)
{
  alpha<- c(alpha,min(1,i/100))
}

d1A  = ArrivalGFDerivative(1,lambda);
d1Gc = ServiceGFDerivative(1,c);
d2Gc = ServiceGF2d(1,c);
d1T  = VacationGFDerivative(1);
d2T  = vacationGF2d(1);
d2A  = ArrivalGF2d(1,lambda);

coeff1 = 0.5;
coeff2 = 20;
coeff3 = 0.25;
coeff4 = 30;

Cost=function(x) #  : age_nu, probability of new connection, backlog, probability server busy
{
  age=x[1];
  p_new_con=x[2];
  bl=x[3];
  psb=x[4];
  #print(x)
  val1 = coeff1*age;
  val2 = coeff2*exp(p_new_con);
  val3 = coeff3*bl;
  val4 = coeff4*exp(psb);
  Cost = val1 + val2 + val3 + val4;
  #cat(sprintf("%f %f %f %f\n", val1,val2, val3,val4));
  #print(Cost)
  Cost
}
fn <- function(th)
{
  Consts    = Constants(c,lambda,alpha,th);
  psb       = ProbServerBusy(th,alpha,Consts);
  p_new_con = StartNewService(th,lambda,alpha,Consts);
  ESystCon  = Re(ESystCont(c,lambda,alpha,Consts,th));
  EServerCon= Re(EServerCont(c,lambda,alpha,Consts,th));
  bl        = ESystCon-EServerCon;
  age       = rpoUB(th,lambda,alpha);
  fn        = Cost(c(age,p_new_con,bl,psb))
}

result = {}
th_vals = 11:max_th#c(11,20,30,40,50,max_th)
for( th in th_vals)
{
  result <- c(result, fn(th))
}

filename = q=paste("plot_tradeoff_load","",".pdf", sep = "");
pdf(filename);
ylim2 = max(result)+2;
ylim1 = min(result)-5;
par(mar=c(5,5,5,2));
plot(th_vals,result, type = "l", xlab = 'restarting threshold (l)', ylab = 'Cost',col="blue", lwd=3.5,cex.lab=1.5,ylim=c(ylim1,ylim2) );
dev.off()




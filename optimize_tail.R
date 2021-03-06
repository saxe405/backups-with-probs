#install.packages(c("alabama","nloptr"), lib = getwd(),repos = "http://cran.us.r-project.org")
require(alabama)
require(nloptr)
rm(list = ls(all.names = TRUE))
source("Systemdefinition.R");
source('Solverfunctions.R');
source('matrixcoefficients.R');
source('performance measures.R');

c=10;

th = 15;

arrivalrates = seq(0.3, 0.7, length.out = 10);
result_alphas ={};
Utz = {};
Utz2 = {};
rhos = {};
count = 1;
for(lambda in arrivalrates)
{
  
  load = ArrivalGFDerivative(1,lambda)*ServiceGFDerivative(1,c)/c;
  rhos = c(rhos, load)
  
  d1A=ArrivalGFDerivative(1,lambda);
  d1Gc = ServiceGFDerivative(1,c);
  d2Gc = ServiceGF2d(1,c);
  d1T= VacationGFDerivative(1);
  d2T= vacationGF2d(1);
  d2A = ArrivalGF2d(1,lambda);
  
  coeff1 = 20;
  coeff2 = 15;
  
  alpha_0 = rep_len(0,th-1);
  alpha_1 = rep_len(1,th-1);
  
  Utility=function(x) #  : age_nu, probability of new connection, backlog, probability server busy
  {
    
    p_new_con=x[1];
    psb=x[2];
    utility = coeff1*exp(p_new_con) +coeff2*exp(psb);
  }
  
  fn <- function(x) 
  {
    alpha = x;
    if(is.null(alpha))
    {
      debug_val = 1;
    }
    Consts =    Constants(c,lambda,alpha,th);
    psb    =    ProbServerBusy(th,alpha,Consts);
    p_new_con = StartNewService(th,lambda,alpha,Consts);
    fn     =    Utility(c(p_new_con,psb));
  }
  hin <- function(x)
  {
    alpha = x # hin >= 0
    hin = {}
    Consts = Constants(c,lambda,alpha,th);
    ESystCon   =   Re(ESystCont(c,lambda,alpha,Consts,th));  
    EServerCon =   Re(EServerCont(c,lambda,alpha,Consts,th));
    bl         =   ESystCon-EServerCon;
    age_tail   =   Age_tail(th,lambda,alpha);
    age        =   rpoUB(th,lambda,alpha);
    hin        <-  c(hin, 20 - bl)  #backlog <= 200
    hin        <-  c(hin, 30 - age)
    hin        <-  c(hin, -15 - age_tail) # log10(pr(age)>30E(T)) <= -10 
    hin
  }
  max_val = Inf;
  max_digit=5;
  for(digit in 2:(th))
  {
    new_alpha=c(rep_len(0,digit-1),rep_len(1,th-digit));
    new_val=fn(new_alpha);
    
    if(new_val<max_val && sum(hin(new_alpha) >=0)==length(hin(new_alpha)))
    {
      max_val=new_val;
      max_digit = digit;
    }
  }
  best_threshold = c(rep_len(0,max_digit-1),rep_len(1,th-max_digit));
  alpha_start = c(rep_len(0.1,max_digit-1),rep_len(1,th-max_digit));
  print("running auglag")
  abcd = auglag(alpha_start, fn, gr = NULL, lower = rep_len(0,th-1), upper = rep_len(1,th-1), hin = hin,localtol = 1e-10, control = list(maxeval = 1000))
  cache = abcd$par
  Utz[count] = fn(cache);
  Utz2[count] = fn(best_threshold);
  print("inequalities : optimal")
  print(cache)
  print(hin(cache))
  print("inequalities : threshold")
  print(best_threshold)
  print(hin(best_threshold))
  count = count + 1;
}
filename = q=paste("plot_Utility_tail_load","",".pdf", sep = "");
pdf(filename);
ylim2 = max(Utz,Utz2 );
ylim2 = ylim2 + 0.5;

ylim1 = min(Utz,Utz2);
par(mar=c(5,5,5,2));

plot(rhos,Utz, type = "l", xlab = 'load', ylab = 'Cost',col="blue", lwd=3.5,cex.lab=1.5,ylim=c(ylim1,ylim2) );
lines(rhos,Utz2,col="green",type="l", lty=1, lwd=2);
plot_colors = c("blue", "green");

names2 = c("optimal probabilistic policy" , "optimal pure threshold policy");
legend("topleft", names2, cex=1.4, col=plot_colors,lty=1, lwd=3.5, bty="n");
dev.off()




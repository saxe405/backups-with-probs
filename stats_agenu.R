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
type = "mod"; #low,mod,high
vacation_probabilities=function()
{
  if(type == "mod")
  {
    #vacation_probabilities=c(0,0,0,0,0,0,0,0,0,0,0,0,0.1,0,0.7,0.1,0,0.1);# mod
    vacat_len=20;
    y = floor(vacat_len/4);
    temp = rep_len(0,vacat_len+1);
    for(j in (vacat_len+2-y):(vacat_len+1))
    {
      temp[j] = 1/(y);
    }
    vacation_probabilities=temp;
  }
  else if(type =="low")
  {
    vacation_probabilities=c(0,0,0,0,0,0.1,0.7,0.1,0.1); #low
  }  
  else
  {  
    #high
    vacat_len = 25;
    temp =rep_len(0, vacat_len+1);
    temp[vacat_len-2]=0.1;
    temp[vacat_len-1]=0.1;
    temp[vacat_len]=0.7;
    temp[vacat_len+1]=0.1;
    vacation_probabilities= temp;
  }
}  
# DEFINE THE SYTEM VARIABLES HERE
c=10;
th = c+5; #threshold denoted by l in notebook.
offset = 0;
alpha = {}; #these are the probabilites of start based on the system content.
for(i in 1:(th-1))
{
  alpha[i] = 0.5+i*i/(2*th*th); #equally likely to start at every pointprint(lambda);
  #alpha[i]=1;
}
len=16;
load=seq(0.2, 0.95, length.out=len);

# Computation of the system (queue+server) content based on the load of the system.

sdrpol = {};
sdrpou = {};
sdrpoapprox = {};
utilityn = {};
for(j in 1:(len-offset))
{
  lambda =              load[j]*c/EServicet(c);
  temp=                 Varrpob(th,lambda,alpha);
  sdrpou[j]=           sqrt(temp[2]);
  sdrpol[j]=           sqrt(temp[1]);
  sdrpoapprox[j]=      temp[3];
  #utilityn[j] = Utility(c(rpoUBv[j],NewService[j],EQCon[j],ProbServerBusi[j]));
}  
#plot(load[1:(len-offset)],NumVacationFPv, type = "l", xlab = 'load', ylab = 'NumVacat',col="blue");

for(i in 1:(th-1))
{
  alpha[i]=0;
}
#th=c;
sdrpol2={};
sdrpou2={};
sdrpoapprox2 = {};
for(j in 1:(len-offset))
{
  lambda =              load[j]*c/EServicet(c);
  temp=                  Varrpob(th,lambda,alpha);
  sdrpou2[j]=           sqrt(temp[2]);
  sdrpol2[j]=           sqrt(temp[1]);
  sdrpoapprox2[j]=      temp[3];
}

for(i in 1:(th-1))
{
  alpha[i]= 1;
}

sdrpou3 ={};
sdrpol3 ={};
sdrpoapprox3={};
for(j in 1:(len-offset))
{
  lambda =              load[j]*c/EServicet(c);
  temp=                  Varrpob(th,lambda,alpha);
  sdrpou3[j]=           sqrt(temp[2]);
  sdrpol3[j]=           sqrt(temp[1]);
  sdrpoapprox3[j]=      temp[3];
}
plot_colors = c("blue", "green", "red");
names = c(expression(alpha[i] == paste(i/2*c)) , expression(alpha[i] == paste(0)), expression(alpha[i] == paste(1)));
system('taskkill /f /im AcroRd32.exe');

## standart deviation  of rpo
filename = q=paste("C:/Users/sapoorv/Downloads/latex/Data Backups/figs/plot_sd_rpo_load_",type,".pdf", sep = "");
pdf(filename);

ylim2 = max(sdrpou,sdrpou2,sdrpou3 );
ylim2 = ylim2 + 0.5;

ylim1 = min(sdrpol,sdrpol2,sdrpol3);
par(mar=c(5,5,5,2));

plot(load[1:(len-offset)],sdrpol, type = "l", xlab = 'load', ylab = 'SD[max-Age]',col="blue", cex.lab=1.5, lwd=3.5,ylim=c(ylim1,ylim2) );
lines(load[1:(len-offset)],sdrpou,col="blue");
lines(load[1:(len-offset)],sdrpol2,col="green",type="l", lty=2, lwd=2);
lines(load[1:(len-offset)],sdrpou2,col="green",type="l", lty=2, lwd=2);
lines(load[1:(len-offset)],sdrpol3,type="l", lty=3, lwd=3.5, col="red");
lines(load[1:(len-offset)],sdrpou3,type="l", lty=3, lwd=4, col="red");
legend("topright", names, cex=1.4, col=plot_colors, lty=1:3, lwd=3.5, bty="n");
dev.off()
 print(sdrpol)
 print(sdrpou)
 print(sdrpol2)
 print(sdrpou2)
 print(sdrpol3)
print(sdrpou3)
#print(paste(sdrpoapprox, collapse = ','))
#print(paste(sdrpoapprox2,collapse = ','))
#print(paste(sdrpoapprox3,collapse = ','))

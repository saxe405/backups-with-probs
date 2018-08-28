ServerCont={};
for(j in 1:10)
{
  lambda=load[j]*c/EServicet(c);
  Consts = Constants(c,lambda,alpha);
  ServerCont[j]=Re(EServerCont(c,lambda,alpha,Consts));
}
#plot(load,ServerCont)

ServerContinBatch={};
for(j in 1:10)
{
  lambda=load[j]*c/EServicet(c);
  Consts = Constants(c,lambda,alpha);
  ServerContinBatch[j]=Re(EServerContinBatch(c,lambda,alpha,Consts));
}
#plot(load,ServerContinBatch)

ServerBusyProb={};
for(j in 1:10)
{
  lambda=load[j]*c/EServicet(c);
  Consts = Constants(c,lambda,alpha);
  ServerBusyProb[j]=Re(ProbServerBusy(c,alpha,Consts));
}
plot(load,ServerBusyProb)
lambda=0.1*c/EServicet(c);
Consts = Constants(c,lambda,alpha);
ServerBusyP=Re(ProbServerBusy(c,alpha,Consts))
val=ServiceGFDerivative(1,c)*v110(c,lambda,alpha,Consts)-ServiceGFDerivative(1,c)*sum(Consts[c:(2*c-1)])+VacationGFDerivative(1)*Consts[c];
for(j in 1:(c-1))
{
  val=val+ServiceGFDerivative(1,j)*Consts[j+c]+VacationGFDerivative(1)*(1-alpha[j])*Consts[j]-ServiceGFDerivative(1,j)*(1-alpha[j])*Consts[j];  
}
val
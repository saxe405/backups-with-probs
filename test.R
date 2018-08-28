#check if the coefficients of the normality condition are +ve
source('Systemdefinition.R');

poissonarrivalrates = seq(0.2, 0.80, length.out = len)
c=10;
th = 30; #threshold denoted by l in notebook.
len=20;

d1Gc = ServiceGFDerivative(1,c);
d1T = VacationGFDerivative(1);
for(j in 1:(c-1))
{
  m = j;
  d1Gm = ServiceGFDerivative(1,m);
  print(m*d1Gc+c*d1T-c*d1Gm);
}  
for(j in c:(th-1))
{
  print( c*d1T);
}

for(j in (th+1):(c+th-1))
{
  m=j-th;
  d1Gm = ServiceGFDerivative(1,m);
  print(c*d1Gm-m*d1Gc);
}


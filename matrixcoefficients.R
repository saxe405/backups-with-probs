#install.packages("BB", lib = getwd(),repos = "http://cran.us.r-project.org")
require("BB")
#######################################################
# Coefficients of d(m)'s
#######################################################
Bmz=function(z,m,c,lambda)
{
  y = ArrivalGF(z,lambda); #A(z)
  if(m==0)
  {
    Bmz=ServiceGF(y,c)*(VacationGF(y)-1)*(1-z^c);
  }
  else
  {
    gmaz = ServiceGF(y,m);
    gcaz = ServiceGF(y,c);
    # Bmz = gmaz*gcaz*(z^m-z^c);
    # Bmz = Bmz + gcaz*(z^m)*(z^c-1);
    # Bmz = Bmz + gmaz*(z^c)*(1-z^m);
    Bmz = (z^c*gmaz-z^m*gcaz)*(1-gcaz) - (z^c-gcaz)*(gmaz-gcaz)*z^m;
  }
  Bmz;
}

#########################################################
# Derivative of the coefficients of d(m)'s
#########################################################
Bmzprime=function(z,m,c,lambda)
{
  y=ArrivalGF(z,lambda);
  if(m==0)
  {
    Bmzprime=ServiceGFDerivative(y,c)*ArrivalGFDerivative(z,lambda)*(1-z^c)*(VacationGF(y)-1);
    Bmzprime=Bmzprime-ServiceGF(y,c)*(c*z^(c-1))*(VacationGF(y)-1);
    Bmzprime=Bmzprime+ServiceGF(y,c)*(1-z^c)*VacationGFDerivative(y)*ArrivalGFDerivative(z,lambda);
  }
  else
  {
    Bmzprime=(c*z^(c-1)-(m+c)*z^(m+c-1))*ServiceGF(y,m)+z^c(1-z^m)*ServiceGFDerivative(y,m)*ArrivalGFDerivative(z,lambda);
    Bmzprime=Bmzprime+((m+c)*z^(m+c-1)-m*z^(m-1))*ServiceGF(y,c)+z^m*(z^c-1)*ServiceGFDerivative(y,c)*ArrivalGFDerivative(z,c);
    Bmzprime=Bmzprime+(m*z^(m-1)-c*z^(c-1))*ServiceGF(y,m)*ServiceGF(y,c)+(z^m-z^c)*(ServiceGFDerivative(y,c)*ServiceGF(y,m)+ServiceGFDerivative(y,m)*ServiceGF(y,c))*ArrivalGFDerivative(z,lambda);
  }
  Bmzprime;
}

###############################################################
# Coefficients of d_0(m)'s
#################################################################
Wmz=function(z,m,c,lambda,alpham) #alpha here is the alpha corresponding to m
{
  y=ArrivalGF(z,lambda);
  gmaz = ServiceGF(y,m);
  gcaz = ServiceGF(y,c);
  taz = VacationGF(y);
  # Wmz = gmaz*gcaz*(z^c-z^m);
  # Wmz = Wmz + gmaz*(z^c)*(z^m-1);
  # Wmz = Wmz + gcaz*taz*(z^m)*(1-z^(c));
  # Wmz = Wmz*(1-alpham);
  Wmz = z^c*(z^m*taz - gmaz)*(1-gcaz)-(z^c-gcaz)*(taz-gmaz)*z^m;
  Wmz = Wmz*(1-alpham);
  Wmz;
}

#######################################################################
# Coefficients of d_0(m) for m>=c and <l
#######################################################################
Rmz=function(z,m,c,lambda,alpham)
{
  y = ArrivalGF(z,lambda);
  gcaz = ServiceGF(y,c);
  taz = VacationGF(y);
  Rmz = (1-gcaz)*((z^c)*taz-gcaz)-((z^c)-gcaz)*(taz-gcaz);
  Rmz = Rmz*(1-alpham)*z^m;
  Rmz;
}

###############################################################
# Derivative of the Coefficients of d_0(m)'s
#################################################################
Wmzprime=function(z,m,c,lambda,alpham)
{
  y=ArrivalGF(z,lambda);
  Wmzprime=(-c*z^(c-1)+(m+c)*z^(m+c-1))*ServiceGF(y,m)+z^c(-1+z^m)*ServiceGFDerivative(y,m)*ArrivalGFDerivative(z,lambda);
  Wmzprime=Wmzprime+(-m*z^(m-1)+c*z^(c-1))*ServiceGF(y,m)*ServiceGF(y,c)+(-z^m+z^c)*(ServiceGFDerivative(y,c)*ServiceGF(y,m)+ServiceGFDerivative(y,m)*ServiceGF(y,c))*ArrivalGFDerivative(z,lambda);
  Wmzprime=Wmzprime-(m+c)*z^(m+c-1)*ServiceGF(y,c)*VacationGF(y)+(1-z^(m+c))*VacationGF(y)*ServiceGFDerivative(y,c)*ArrivalGFDerivative(z,lambda);
  Wmzprime=Wmzprime+ServiceGF(y,c)*VacationGFDerivative(y)*ArrivalGFDerivative(z,lambda)*(1-z^(m+c));
  Wmzprime=Wmzprime+(z^(m+c)-z^c)*VacationGFDerivative(y)*ArrivalGFDerivative(z,lambda)+VacationGF(y)*((m+c)*z^(m+c-1)-c*z^(c-1));
  Wmzprime = Wmzprime*(1-alpham);
  Wmzprime;
}

###################################################################
# This returns the coefficients of the generating function equation
###################################################################
Constants=function(c,lambda,alpha,th)
{  
  #ro=0.75; defined as lambda*E(G_c)/c
  results ={};
  results2={};
  
  #x0=2.0 + 1.0i;
  #x1=3.0 + 2.0i;
  #x2=1.0 + 3.0i;
  #x3=0.0 + 0.0i;
  #guesses<-c(x0,x1,x2,x3);
  guesses = complex(real=runif(4, 1, 3),imaginary = runif(4, 3, 6));
  
  k=1;
  #print("end")
  while(k<c)
  {
    guess1 =complex(modulus=1.0,argument = 2*pi*k/c);
    #result1=newt(guess1,myFunction,myFunctionsDerivative,(10^-10),lambda,c,k, 200);
    #result1=mullerfn(myFunction, guess1,(10^-100),(10^5),c,lambda,k);
    result1 = BBsolve(par=c(Re(guess1),Im(guess1)), fn=myFunction2,control = list(maxit = 5000, tol = 10^-15), quiet = TRUE, params=c(c,lambda,k) )
    #print(result1)
    #results[k]=result1;
    #print("next")
    results[k]=complex( real = result1$par[1], imaginary = result1$par[2]);
    #print(myFunction(c,result1,lambda,k));
    #result1=mullerfn(myFunction, guesses,(10^-100),(10^5),c,lambda,k);
    #results2[k]=result1;
    k=k+1;
  }
  #print("end2")
  #plot(results);
  for(j in 1:(c-1))
  {
    if( abs(results[j]^c - ServiceGF(ArrivalGF(results[j],lambda),c)) > 10^-3)
    {
      print('we got a problem!')
      print(results[j])
    }
  }
  A = matrix( 0, nrow=c+th-1,ncol=c+th-1);
  b = matrix( 0, nrow=c+th-1,ncol=1);
  b[c+th-1]=(c- ServiceGFDerivative(1,c)*ArrivalGFDerivative(1,lambda) ); #corresponding to normalization condition
  
  tAvec = {};
  
  for(i in 1:th)
  {
    tAvec[i]=Tai(i-1,lambda); #tAvec[i] =ta(i-1) because we dont have 0 subscript in R
  }  
  #tAvec[c]=Tai(c-1,lambda);
  
  for (i in 1:(th-1)) #skip the normalization row and d(m)s
  {
    temp=matrix( 0, nrow=c+th-1,ncol=1);
    for(j in 1:i)
    {
      if( is.na(alpha[j]) || is.na(tAvec[i-j+1]))
      {
        debug_val = 1;
      }
      temp[j]=(-1*(1-alpha[j])*tAvec[i-j+1]);
      if(i==j)
        temp[j]= temp[j]+1;
    }
    temp[th]=-1*tAvec[i+1] #corresponding to d(0) t_A(i) is tAvec(i+1)
    A[i,]=temp;
  }
  for(i in th:(c+th-2))
  {
    temp=matrix( 0, nrow=c+th-1,ncol=1);
    for( j in 1:(c-1))
    {
      temp[j]=Wmz(results[i-th+1],j,c,lambda,alpha[j]);  
    }  
    for( j in c:(th-1))
    {
      temp[j]=Rmz(results[i-th+1],j,c,lambda,alpha[j]);
    }
    for(j in th:(c+th-1))
    {
      temp[j]=Bmz(results[i-th+1],j-th,c,lambda);#(z,m,c,lambda)
    }
    A[i,]=temp;
  }  
  #this gives us the coefficients of the unknowns.
  #A[2*c-1][j]
  
  d1Gc = ServiceGFDerivative(1,c);
  d1T = VacationGFDerivative(1);
  for(j in 1:(c-1))
  {
    m = j;
    d1Gm = ServiceGFDerivative(1,m);
    temp[j] = m*d1Gc+c*d1T-c*d1Gm;
    temp[j] = temp[j]*(1-alpha[j]);
  }  
  for(j in c:(th-1))
  {
    temp[j] = c*d1T*(1-alpha[j]);
  }
  temp[th]=VacationGFDerivative(1)*c;
  for(j in (th+1):(c+th-1))
  {
    m=j-th;
    d1Gm = ServiceGFDerivative(1,m);
    temp[j] = c*d1Gm-m*d1Gc;
  }
  A[c+th-1,]=temp;
  #print(tAvec)
  constants = solve(A,b);
  
  #constants = Re(constants);
  #print(A%*%constants - b)
}
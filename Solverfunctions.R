
myFunction2 = function(x, params){
  c = params[1]
  lambda = params[2]
  k = params[3]
  z <- complex( real=x[1], imaginary=x[2])
  y = ArrivalGF(z,lambda)
  epsilon = complex(modulus=1.0,argument = 2*pi*k/c);
  f <- z - 1.0*y*((1/3+y/3+y^2/3)^(1/c))*epsilon;
  myFunction2 = c(Re(f), Im(f))
}
######################################################################
# The function for which we want to compute zeros of
######################################################################
myFunction=function(c,xn,lambda,k) #k defines the root of unity to use
{
  yn = ArrivalGF(xn,lambda);
  epsilon = complex(modulus=1.0,argument = 2*pi*k/c);
  myFunction = xn-1.0*(ServiceGF(yn,c)^(1.0/c))*epsilon; 
}
######################################################################
# myfunction_sing : the function (1-alpha)t0(z) = 1 to find the dominant singularity
######################################################################
myFunction_sing = function(x, params){
  alpha = params[1]
  A0 = ArrivalGF(0,params[2])
  myfunction_sing <- VacationGF(x*A0)-1/(1-alpha);
}
######################################################################
# Derivative of the function we want to optimize
######################################################################
myFunctionsDerivative=function(c,xn,lambda,k)
{
  yn = ArrivalGF(xn,lambda);
  epsilon = complex(modulus=1.0,argument = 2*pi*k/c);
  myFunctionsDerivative = (1.0 - (1.0/c)*(ServiceGF(yn,c)^(-1.0+(1.0/c)))*ServiceGFDerivative(yn,c)*ArrivalGFDerivative(xn,lambda)*epsilon);
}

######################################################################
# Newton solver
######################################################################
newt=function(xn,f,fDash,accu,lambda,c,k, max_iterations)
{
  xnplus1 = xn - ( f(c,xn,lambda,k) / fDash(c,xn,lambda,k) );
  if (Mod(xn - xnplus1) > accu && (Mod(f(c,xnplus1,lambda,k))>0) && (max_iterations > 0)) {
    newt= newt(xnplus1, f, fDash, accu,lambda,c,k,max_iterations-1)
  } else {
    newt=xnplus1;
  }
}

######################################################################
# Muller solver
######################################################################
mullerfn=function(f,x,allerr,maxmitr,c,lambda,k)
{
  itr=1;
  #fxk, fxk1, fxk2, a, b;
  while(itr<=maxmitr)
  {
    fxk  = f(c,x[3],lambda,k);
    fxk1 = f(c,x[2],lambda,k);
    fxk2 = f(c,x[1],lambda,k);
    
    b= (fxk-fxk1)/(x[3]-x[2]);
    a= (fxk1-fxk2)/(x[2]-x[1]);
    a= a-((b-a)/(x[3]-x[1]));
    b= b+a*(x[3]-x[2]);
    
    val=sqrt(b*b-4.0*a*fxk);
    x[4] = x[3]-2.0*fxk/(b+val); 
    if(Mod(f(c,x[4],lambda,k))<allerr)
    {
      mullerfn = x[4];
      break;
    }  
    x[1]=x[2];
    x[2]=x[3];
    x[3]=x[4];
    if(Mod(x[2]-x[3])==0 || Mod(x[1]-x[2])==0)
      x[2]=complex(real=runif(1,1,4),imaginary = runif(1,3,5));
    if(Mod(x[3]-x[1])==0)
      x[1]=complex(real=runif(1,1,2),imaginary = runif(1,4,6));
    itr = itr+1;
  }
  mullerfn=x[4];
}
calculate_zeroes=function(lambda,service_parameters,n)
{
  tol = 10^(-9); # tolerance
  k = 0; #ANDERS!k = 1; number of zeros sought ?
  number_of_zeroes_found = 0;  #ANDERS!!!!number_of_zeroes_found = 1; holds some zeros found at
  #zeroes(1) = 1+0*i; #The first zero is not 1.
  zeroes = integer(n);
  #MORE SO NOW !!!!!!!!!!!!!
    
  stepmax = 0.1; #Max step size; apparently has its uses in part 2
  rand_max = 32767;
  
  startvalue = 0+0*i;
  
  #-------------------------------------------------------------------------
  # PART 1: the n zeros search by vglen n = B z (z) ^ (1 / n) exp (i * 2 * pi * k / n)
  # 0 <= k <n to solve (for k = 0, the solution must be equal to 1)
  # Any such equation has exactly one solution and is such a solution
  # Corresponds to exactly one of the n solutions of z ^ n-B (z) = 0; 1-1, there is a
  # Image between the n solutions of the equations and the n
  # Of solutions z ^ n-B (z)
  #-------------------------------------------------------------------------
    
  while(k < n)
  {
    rotor = exp(1i*2*pi*k/n);
    iterations = 1; 
    if( number_of_zeroes_found == 0)
    {
      z_curr = 1.0+0.0*i;
    }
    else
    {
      z_curr = zeroes[number_of_zeroes_found]; #START value for the algorithm
    }
  
    z_last = 2.0+2.0*i; 
    while (abs(z_curr-z_last) > tol & iterations < 500)
    {  
      z_last = z_curr;
      Bz = B(lambda,service_parameters,n,z_curr); #B(z_curr)
      dBz = dB(lambda,service_parameters,n,z_curr); #B'(z_curr)
      delta1 = z_curr-Bz^(1/n)*rotor; #stemt match f (z_curr) from formula NR
      delta2 = 1.0 - Bz^(1/n)*rotor*dBz/(Bz*n); #stemt corresponds to f '(z_curr) from formula NR
      #testen or the derivative is or is not equal to 0
      if (abs(delta2) < tol) #Select a random direction in order to get away from
      {  
        mod = 0.001; #will but very small
        arg = 2*pi*random('unid',rand_max,1,1)/rand_max;
        z_curr = z_curr + mod*exp(1i*arg);
      }
      else
      {
        z_curr = z_curr -(delta1/delta2);
      }
  
      iterations = iterations+1;
      print('pass1');
      print(delta1);
      print(delta2);
    }
    print('pass11');
    if (abs((z_curr^n)-B(lambda,service_parameters,n,z_curr)) <= tol) #In this case z_curr zero
    {
      zeroes[number_of_zeroes_found+1] = z_curr;
      number_of_zeroes_found = number_of_zeroes_found+1;
    }
    else #proberen with acceleration Aitken; this can sometimes provide relief
    {
      #vergeet also not to work here with the rotor, otherwise
      #riskeer find your roots that have been found !!!
      iterations = 1;
    }  
    fz_0 = abs((z_curr^n)-B(lambda,service_parameters,n,z_curr));
    while( (abs(fz_0) > tol) & (iterations <= 1000) )
    {
      fz_0 = z_curr - B(lambda,service_parameters,n,z_curr)^(1/n)*rotor;
      z_1 = z_curr + fz_0;
      fz_1 = z_1 - (B(lambda,service_parameters,n,z_1))^(1/n)*rotor;
      z_2 = z_1 + fz_1;
      z_curr = z_2 - (z_2-z_1)^2/(z_curr-2*z_1+z_2);
      fz_0 = z_curr - (B(lambda,service_parameters,n,z_curr))^(1/n)*rotor;
      iterations = iterations + 1;
      print('pass2');
    }
  
    if ((abs((z_curr^n)-B(lambda,service_parameters,n,z_curr)) <= tol && abs(z_curr)<=1)) #in dit geval is z_curr het nulpunt
    {
      zeroes[1+number_of_zeroes_found] = z_curr;
      number_of_zeroes_found = number_of_zeroes_found+1;
    }
    else #de startwaarde van z_curr opslaan en opnieuw gebruiken in het
    {
      #tweede deel
      startvalue = rotor; 
    }
  }
  k = k+1;
  
 
  
# ------------------------------------------------- --------------
#Second PART: THE OTHER zeros SEARCH: ALSO THROUGH Newton-Raphson
#MAAR PREVIOUS zeros NOW BE LEFT TO PAY
# ------------------------------------------------- --------------
 
    startvalue = 0;
    #Opgelet: Here must not be worked with a rotor
    #want last piece: n vglen each with one solution that you are looking for
    #hier 1 cf. n solutions of which you have already j !!!!!
    #indien you would continue working with the j vglen that have not been resolved,you do not 
    #kun other zeros road sections because no zeros
    #zijn that separate cf.!
    
    while( number_of_zeroes_found < n) #zolang found not all zeros
    {
      iterations = 0;
      z_curr = startvalue;#0+0*i;
      z_last = 2.0+2.0*i;
  
  
      while (abs(z_curr-z_last)>tol & iterations < 5000)
      {
        z_last = z_curr;
  
        delta1 = z_curr^n - B(lambda,service_parameters,n,z_curr);
        delta2 = n*z_curr^(n-1)-dB(lambda,service_parameters,n,z_curr);
  
        #displays in these five lines that the previous zeros are left handed
        l = 0;
        while (l < number_of_zeroes_found )
          delta2 = delta2 - delta1/(z_curr-zeroes[l+1]);
          l = l+1;
        end
  
          #testen or the derivative is or is not equal to 0
        if( abs(delta2)<tol) #Select a random direction in order to get away from
        {
          mod = stepmax*0.01;
          arg = 2*pi*random('unid',rand_max,1,1)/rand_max;
          z_curr = z_curr + mod*exp(1i*arg);
        }
        else
        {
          p = -delta1/delta2;
        }
        if(abs(p)>stepmax)
        {
          p = p*stepmax/abs(p); 
        }
        z_curr = z_curr+p;
      }
      iterations = iterations + 1;
      print('pass3');
    }
  
    if( abs(z_curr^n-B(lambda,service_parameters,n,z_curr))<=tol&&abs(z_curr)<=1) #maw de opl voldoet
    {
      zeroes[number_of_zeroes_found+1] = z_curr;
      number_of_zeroes_found = number_of_zeroes_found+1;
    }
  
  #testen found or the zero points lie inside the unit circle
  index = 1;
  while (index <= n)
  {
    if (abs(zeroes[index]) > 1 && abs(abs(zeroes[index])-1)> tol)
    {
      zeroes[index]
      abs(zeroes[index])
      warning('not all zeros lie on the unit circle')
    }
    index = index+1;
  }
  
  #testen or found zero or differ
  #this slow the program considerably
  index1 = 1;
  while( index1 <= n-1)
  {
    index2 = index1 + 1;
    while(index2 <= n)
    {
        if(abs(zeroes[index1]-zeroes[index2]) < 10^(-6))            
        warning('some zeros are frequently found')
    }
    index2 = index2 + 1;
  }
  index1 = index1 + 1;

  
  
  #de volgende lijn code schrijft de nulpunten weg naar een file
  #dlmwrite('zeroes.data',[zeroes']);
                           
calculate_zeroes = zeroes; 
}                         
                           
                           
                           
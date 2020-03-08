

function [x,y] = nonhomopp(intens,pos,T)
  % example of generating a 
  % nonhomogeneousl poisson process on [0,T] with intensity function intens

  if length(T)==1
    x = 0:.1:T;
  else
    x = T;
  end
  
  m = intens(x);              %% process-intensity at timepoints x
  m2 = max(m);                %% generate homogeneouos poisson process (with rate = maximum rate)
  
  u = rand(1,ceil(1.5*T(end)*m2)); %% generate a sufficient amount of random variables to cover the whole time
  
  y = cumsum(-(1/m2)*log(u)); %% generate points of homogeneous pp
  y = y(y<T(end));                 %% select those points less than T
  n=length(y);
  
  m = intens(y);              %% evaluates intensity function at homogeneous pp points
  
  y = y(rand(1,n)<m/m2);      %% filter out some points with prob according to ratio of nonhomo/homo pp
  
  figure
  hold on
  plot(x,intens(x),'k')
  y
  intens
  size(y)
  size(intens(y))
  scatter(y,intens(y),'ro','filled')
%    hist(y,10)

end
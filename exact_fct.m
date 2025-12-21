function [val,dval,valp,dvalp]=exact_fct(xx) 
%val-->u, ddval-->Delta u=p
x=xx(1);
y=xx(2); 
R=(x^2+y^2)^(1/2);
val=R^3; 
dval=3*[x,y]*R; 
valp=9*R;
dvalp=9*[x,y]/R;

 
  
 



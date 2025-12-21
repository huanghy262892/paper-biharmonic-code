function val=source_fct(xx)  
%val-->u, ddval-->Delta u=p
x=xx(1);
y=xx(2);
R=(x^2+y^2)^(1/2); 
val=9/R;

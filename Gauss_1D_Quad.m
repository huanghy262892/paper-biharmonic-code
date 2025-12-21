function [wg,xg] = Gauss_1D_Quad(n)
% GL_RULE 1D Gauss-Legendre quadrature rule.
%
%   QUADRULE = GAULEG(a,b,n) computes the N-point Gauss-Legendre
%   quadrature rule on the interval [a,b] up to mashine precision. 
%
%   Note that all quadrature rules obtained from GAULEG are of order 2*N-1.
%
%   The struct QUADRULE contains the following fields:
%    W N-by-1 matrix specifying the weights of the quadrature rule.
%    X N-by-1 matrix specifying the abscissae of the quadrature rule.
%   
%   Example:
%
%   QuadRule = gauleg(0,1,10);

if (n==1)
    xg = 0; 
    wg = 2;
else
    c = zeros(n-1,1);
    for i=1:(n-1)
        c(i)=i/sqrt(4*i*i-1); 
    end
    J=diag(c,-1)+diag(c,1); [ev,ew]=eig(J);
    xg=diag(ew); 
    wg=(2*(ev(1,:).*ev(1,:)))';
end
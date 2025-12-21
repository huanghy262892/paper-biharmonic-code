function [wg,xg] = Duffy(QuadRule)
% DUFFY Duffy transformation.
%
%   QUADRULE = DUFFY(QUADRULE) computes a Gauss quadrature on the reference
%   triangle using the duffy trick and the 1D Gauss quadrature QUADRULE.
%
%   The struct QUADRULE should at least contain the following fields:
%    W N-by-1 matrix specifying the weights of a 2D quadrature on the unit
%      square.
%    X N-by-2 matrix specifying the abscissae of a 2D quadrature rule on
%      the unit square.
%
%   Example:
%
%   QuadRule = Duffy(TProd(gauleg(0,1,4)));

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nGauss = size(QuadRule.w,1);
 
  % Preallocate memory

  wg = zeros(nGauss,1);
  xg = zeros(nGauss,2);
  
  % Compute Gauss quadrature on reference element
  
  k = 0;
  for i = 1:nGauss
    wg(i) = QuadRule.w(i)*(1-QuadRule.x(i,1));
    xg(i,1) = QuadRule.x(i,1);
    xg(i,2) = QuadRule.x(i,2)*(1-QuadRule.x(i,1));
  end
  
  % Assign output arguments
  
%   QuadRule.w = wg;
%   QuadRule.x = xg;
  
return
% paper title"Symmetric Direct Discontinuous Galerkin Method for the Biharmonic Equation with
%             Non-homogeneous Boundary Conditions"
% Authors: Hongying Huang, Huanhuan Wang, Lin Zhang
% Journal name: Computers and Mathematics with Applications
%%%%  Example 5.3  %%%%%%
% Direct discontinuous Galerkin method for Biharmonic equation: Delta^2 u=f
% exact solution: u=(x^2+y^2)^(3/2)
%domain: [AX,BX]*[CY,DY]=[0,1]*[0,1]
% val_penal---> [b_1, b_2] in the paper
% Beta2---> [b_3, b_4] in the paper
% K---> degree of piecewise polynomial
% errL2---> error of approximate solution in L^2-norm
% errenergy---> error of approximate solution in energy norm


%function [errL2,errH1,errenp,errenergy]=DDG_Biharmonic(NX,K,val_penal1,beta)
NX=20;
K=3;
val_penal1=[0.9,0.9]; 
beta=[0.1,0.1]*0;
NY=NX;
AX=0;
BX=1.0;
CY=0;
DY=1;
EN=2*NX*NY; 
Nloc=(K+2)*(K+1)/2;  
NDOFS=EN*Nloc; 

val_penal=[K^6,K^2].*val_penal1;
Beta2=beta.*[K^2,1/K^2];
xy=coordinate(AX,BX,CY,DY,NX,NY);
elem= element(NX,NY,EN);  
ete=element_element(EN,NDOFS,elem);  
[G_ran,G_col,G_matr,rvec]=globalmat_Biharmonic(EN,K,val_penal,Beta2,ete,xy,elem);
matr_sparse=sparse(G_ran,G_col,G_matr,NDOFS,NDOFS); 
uapprox=matr_sparse\rvec';  
[errmax,errL2,errH1,errenp,errenergy]=compute_approxerror_norm_Biharmonic(K,elem,xy,uapprox,val_penal,ete);
NX
errL2
errenergy
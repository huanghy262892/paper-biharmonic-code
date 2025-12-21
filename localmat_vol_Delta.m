 function [Aloc,Floc]=localmat_vol_Delta(K,xy)
 % int_E Delta u Delta v dxdy, \int_E fvdxdy
 	% local variables
Nloc=(K+2)*(K+1)/2; 
[BE,inv_BE,det1]=compute_triangle_jacobimatrix(xy); 
Aloc=zeros(1,Nloc*Nloc); 
Floc=zeros(1,Nloc); % loop over the quadrature points
%[wg,xg]=Gauss_triangle_nodes_weights_ref(13);
 [wg,xg]=Duffy(TProd(gauleg(0,1,2*K))); % gauleg:[0,1]=[a,b]
 ng = length(wg);
for ig=1:ng		% compute values and derivatives of basis functions and determinant	
    [val_basis,~,Lbasis,~]=elem_basisvalue2(K,xg(ig), xg(ig + ng),inv_BE);
	% compute global coordinates of quadrature point
    xx=xy(1,1:2)'+BE*[xg(ig);xg(ig+ng)];	 	 
	% get source function
	source=source_fct(xx);
	% compute the entries of local matrix
    for idofs=1:Nloc
        for jdofs=1:Nloc 
            rr=Lbasis(jdofs)*Lbasis(idofs)*det1*wg(ig);
			Aloc((idofs-1)*Nloc + jdofs)=Aloc((idofs-1)*Nloc + jdofs)+ rr;
        end
			Floc(idofs )=Floc(idofs)+ source*val_basis(idofs) * det1*wg(ig);
    end
end
	 

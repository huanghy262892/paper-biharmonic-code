function [Bloc11,Bloc12,Floc]=localmat_face_Biharmonic(K,val_penal,Beta2,E1,index_E1,xy1,normal_vec,E2,xy2)
% val_penal=[beta_0,beta_1]; Beta2=[beta_2,beta_3]
Nloc=(K+2)*(K+1)/2;  
normal_vec=normal_vec';
he=norm(normal_vec);
val_penal2=val_penal/(he^2);
[~,inv_BE2,~]=compute_triangle_jacobimatrix(xy2);  
[BE1,inv_BE1,det1]=compute_triangle_jacobimatrix(xy1); 	 
Bloc11=zeros(1,Nloc*Nloc);
Bloc12=zeros(1,Nloc*Nloc);
Floc=zeros(1,Nloc); 
%[xg,wg]= Gauss_nodes_weights_1D(9);
%[xg,wg]= Gauss_1D_nodes_weights_ref(6);
[wg,xg] = Gauss_1D_Quad(2*K); 
ng = length(wg);
for ig=1:ng 
	wg(ig) = wg(ig) / 2.0;
	xg(ig) = 0.5 + 0.5*xg(ig);
end
if (E1 == E2) 
	bn_val = 1.0;  
	bn_par = 0.0;  
	bn_par1 = 0.0;  
else 
	bn_val = 1.0/2;
	bn_par = 1.0/2;
	bn_par1 = 1.0;
end	     
% loop over the quadrature points
for ig=1:ng 
	% compute local coordinates of quadrature point		 
	[ss1,ss2]=loc_coor_quad(index_E1, xy1,BE1,xy2,inv_BE2, xg(ig));
	% compute normal vector to E1 and length of edge		
	% compute values and derivatives of basis functions and determinant
    xx=xy1(1,1:2)'+BE1*ss1(1:2)';  
    [val_basis1,der_basis1,Lbasis1,dLbasis1]=elem_basisvalue2(K,ss1(1), ss1(2), inv_BE1);
    [val_basis2,der_basis2,Lbasis2,dLbasis2]=elem_basisvalue2(K,ss2(1), ss2(2), inv_BE2); 
	% compute the entries of local matrix Bloc11
    for idofs=1:Nloc %for v
        for jdofs=1:Nloc % for u
            du=[der_basis1(jdofs),der_basis1(Nloc + jdofs)]*normal_vec; 
            dv=[der_basis1(idofs), der_basis1(Nloc + idofs)]* normal_vec;
            dLu=[dLbasis1(jdofs),dLbasis1(Nloc + jdofs)]*normal_vec; 
            dLv=[dLbasis1(idofs),dLbasis1(Nloc + idofs)]* normal_vec;
            r0=val_penal2(1)*val_basis1(idofs)* val_basis1(jdofs) * wg(ig);
            r1=val_penal2(2)*du*dv *wg(ig);            
			r2=(dLv*val_basis1(jdofs)+dLu*val_basis1(idofs))*wg(ig) *bn_val;
            r3=-(Lbasis1(jdofs)*dv+Lbasis1(idofs)*du)*wg(ig)*bn_val;
            r4=Beta2(1)*(Lbasis1(idofs)*val_basis1(jdofs)+Lbasis1(jdofs)*val_basis1(idofs))*wg(ig)*bn_par1;
			r5=Beta2(2)*(dLu*dv+dLv*du)*wg(ig)*bn_par1; 
            nn=(idofs-1)*Nloc + jdofs;
            Bloc11(nn) = Bloc11(nn)+r0+r1+r2+r3+r4+r5; 
        end
    end 	 
    % compute the entries of local matrix Bloc12
    for idofs=1:Nloc
        for jdofs=1:Nloc 
            du2=[der_basis2(jdofs),der_basis2(Nloc + jdofs)]*normal_vec; 
            dv=[der_basis1(idofs), der_basis1(Nloc + idofs)]* normal_vec;
            dLu2=[dLbasis2(jdofs),dLbasis2(Nloc + jdofs)]*normal_vec; 
            dLv=[dLbasis1(idofs),dLbasis1(Nloc + idofs)]* normal_vec;
            r0=-val_penal2(1)*val_basis1(idofs)* val_basis2(jdofs) * wg(ig)*bn_par1;
            r1=-val_penal2(2)*du2*dv *wg(ig)*bn_par1;            
			r2=(-dLv*val_basis2(jdofs)+dLu2*val_basis1(idofs))*wg(ig) *bn_par;
            r3=(-Lbasis2(jdofs)*dv+Lbasis1(idofs)*du2)*wg(ig)*bn_par;
            r4=-Beta2(1)*(Lbasis1(idofs)*val_basis2(jdofs)+Lbasis2(jdofs)*val_basis1(idofs))*wg(ig)*bn_par1;
			r5=-Beta2(2)*(dLu2*dv+dLv*du2)*wg(ig)*bn_par1; 
            nn=(idofs-1)*Nloc + jdofs;
            Bloc12(nn) = Bloc12(nn)+r0+r1+r2+r3+r4+r5;
        end
    end
    if (E1 == E2) 
		xx=xy1(1,1:2)'+BE1*ss1(1:2)';  
        [val,dval,~,~]=exact_fct(xx);
        for idofs=1:Nloc
            dv=[der_basis1(idofs),der_basis1(Nloc + idofs)]* normal_vec;
            dLv=[dLbasis1(idofs),dLbasis1(Nloc + idofs)]* normal_vec;
            g1=dval*normal_vec;
            r0=val_penal2(1)*val_basis1(idofs)* val * wg(ig) ; 
            r1=val_penal2(2)*g1*dv*wg(ig);
            r2=val*dLv* wg(ig)-g1*Lbasis1(idofs)*wg(ig) ;            
            Floc(idofs) =Floc(idofs)+r0+r1+r2;
        end
    end
end



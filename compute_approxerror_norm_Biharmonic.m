function [errmax,errL2,errH1,errenp,errenergy]=compute_approxerror_norm_Biharmonic(K,elem,xy,uapprox,val_penal,ete)
EN=length(elem);
Nloc=(K+1)*(K+2)/2;
xy_loc = [0,0;1,0;0,1];  
%[wg,xg]=Gauss_triangle_nodes_weights_ref(13);
%[wg,xg]=P73O19();
[wg,xg]=Duffy(TProd(gauleg(0,1,2*K))); % gauleg:[0,1]=[a,b] 
%[xg1,wg1]= Gauss_nodes_weights_1D(9);
%[xg1,wg1]= Gauss_1D_nodes_weights_ref(6);
[wg1,xg1] = Gauss_1D_Quad(2*K);  
ng = length(wg);
ng1=length(wg1); 
INV_BE1=zeros(2,2); 
errmaxnode = 0.0;
facen=[1,2;2,3;3,1];
for ii=1:EN
    tt=elem(ii,1:3);
    xy_vex=xy(tt,1:2); 
    for jj=1:3
		[val,~]=exact_fct(xy_vex(jj,:)); 
        [vbasis,dbasis]=elem_basisvalue(K,xy_loc(jj,1), xy_loc(jj,2), INV_BE1);
		for  mm=1:Nloc   
            val = val - uapprox((ii-1)*Nloc + mm) *vbasis(mm); 
        end
		if  abs(val)>errmaxnode
            errmaxnode = abs(val);
        end
    end
end
errL2 = 0.0;	
errH1 = 0.0;
errmax = 0.0;
errenp=0;
err_jump_u=0; 
err_jump_un=0; 
for ii=1:EN
    tt=elem(ii,1:3);
    xy_vex=xy(tt,1:2); 
    [BE1,INV_BE1,det1]=compute_triangle_jacobimatrix(xy_vex);
	for jj=1:ng  
        xx=xy_vex(1,1:2)'+BE1*[xg(jj);xg(jj+ng)];          
        [val,dval,valp,dvalp]=exact_fct(xx); 		 
        [vbasis,dbasis,Lbasis,dLbasis]=elem_basisvalue2(K,xg(jj),xg(jj+ng), INV_BE1);
	 	for mm=1:Nloc
			val = val - uapprox((ii-1)*Nloc + mm) * vbasis(mm);
			dval(1) = dval(1) - uapprox((ii-1)*Nloc + mm) * dbasis(mm);
			dval(2) = dval(2) - uapprox((ii-1)*Nloc + mm) * dbasis(mm + Nloc);
            valp=valp- uapprox((ii-1)*Nloc + mm) * Lbasis(mm);
        end
		if (abs(val)>errmax) 
             errmax = abs(val);
        end
		det1 = abs(det1);
		errL2 = errL2 + wg(jj) * det1*val*val;
		errH1 = errH1 + wg(jj) * det1*(dval(1)*dval(1)+dval(2)*dval(2));
        errenp=errenp+ wg(jj) * det1*valp*valp;
    end
     %  
    for jj=1:3
        tt = elem(ete(ii,jj),1:3);
        xy_vex1=xy(tt,1:2); 
        [~,INV_BE2,~]=compute_triangle_jacobimatrix(xy_vex1);         
        for ig=1:ng1
            lam=0.5*(1+xg1(ig));
            [ss1,ss2]=loc_coor_quad(jj, xy_vex,BE1,xy_vex1,INV_BE2, lam);             
            xx=xy_vex(1,1:2)'+BE1*ss1'; 
            [vbasis1,dbasis1]=elem_basisvalue(K,ss1(1),ss1(2), INV_BE1);
            [vbasis2,dbasis2]=elem_basisvalue(K,ss2(1),ss2(2), INV_BE2);
            val1=0;
            val2=0;
            adval1=[0,0];
            adval2=[0,0];
            for mm=1:Nloc
                t1=(ii-1)*Nloc;
                t2=(ete(ii,jj)-1)*Nloc;
			    val2 =val2+uapprox(t2 + mm) * vbasis2(mm);
                val1=val1+uapprox(t1 + mm) * vbasis1(mm); 
                adval1(1) = adval1(1) - uapprox(t1 + mm) * dbasis1(mm);
			    adval1(2) = adval1(2) - uapprox(t1 + mm) * dbasis1(mm + Nloc);
                adval2(1) = adval2(1) - uapprox(t2 + mm) * dbasis2(mm);
			    adval2(2) = adval2(2) - uapprox(t2 + mm) * dbasis2(mm + Nloc);
            end  
            line_xy=xy_vex1(facen(jj,2),1:2)-xy_vex1(facen(jj,1),1:2);
            normal_vec=[line_xy(2);-line_xy(1)];
            he=norm(normal_vec);
            val_penal2=val_penal/he/he;
            dvaln=(adval2-adval1)*normal_vec;
            if ii==ete(ii,jj)
                [val2,~,~,~]=exact_fct(xx);                 
            else
                err_jump_u=err_jump_u+(val1-val2)^2*val_penal2(1)*wg1(ig)*0.5;
                err_jump_un=err_jump_un+(dvaln)^2*val_penal2(2)*wg1(ig)*0.5;
            end
        end        
    end   
    
end 
errH1=sqrt(errH1); %seminorm
errL2=sqrt(errL2);
errenergy=sqrt(errenp+err_jump_u+err_jump_un); 
errenp=sqrt(errenp);

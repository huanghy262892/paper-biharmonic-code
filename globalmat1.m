function [G_ran,G_col,G_matr,rvec]=globalmat1(EN,K,val_penal,Beta2,ete,xy,elem)
%compute the matrix about \int_K\nable u\cdot\nable v\dd x and the right vector \int_K fv dx

%facen= [0, 1; 1, 2 ; 2, 0]; 
tt=zeros(1,3);
xy_vex=zeros(3,2);
loc_xy2=xy_vex;
Nloc=(K+2)*(K+1)/2;  
Nloc2=Nloc*Nloc;
matr=zeros(1,EN*Nloc2);
G_matr=zeros(1,4*EN*Nloc2);
G_ran=G_matr;
G_col=G_ran;
rvec=zeros(1,EN*Nloc);
[ran,col]=generate_submatrix_index(Nloc);
nBloc=0;
for ii=1:EN
    tt=elem(ii,1:3);
    xy_vex=xy(tt,1:2); 
	[ Aloc, Floc]=localmat_vol(K,xy_vex);    
	mm = (ii-1)*Nloc2+(1:Nloc2);
	ss = (ii-1)*Nloc+(1:Nloc);
    rvec(ss)=rvec(ss)+Floc;
    matr(mm)=matr(mm)+Aloc; 
    for jj=1:3
		if (jj==3)
            s0 = 1;
        else
            s0=jj+1;
        end
		normal_vec(1) = xy_vex(s0,2) - xy_vex(jj,2);
		normal_vec(2) = xy_vex(jj,1) - xy_vex(s0,1);        
        tt=elem(ete(ii,jj),1:3);
        loc_xy2=xy(tt,1:2);   
		[Bloc11, Bloc12,Floc]=localmat_face(K,val_penal,Beta2,ii,jj,xy_vex, normal_vec,ete(ii,jj), loc_xy2);
		mm = (ii-1)*Nloc2+(1:Nloc2); 
        ndofs=(ii-1)*Nloc+(1:Nloc);
        rvec(ndofs) = rvec(ndofs) - Floc;
        matr(mm) = matr(mm) + Bloc11;
        nBloc=nBloc+1;
        nBlocs=(nBloc-1)*Nloc2+(1:Nloc2);
        G_matr(nBlocs)=Bloc12;  
        G_ran(nBlocs)=ran+(ii-1)*Nloc;
        G_col(nBlocs)=col+(ete(ii,jj)-1)*Nloc;
    end
end
% 把matr组装到全局刚度矩阵中
for ii=1:EN
    nBloc=nBloc+1;
    nBlocs=(nBloc-1)*Nloc2+(1:Nloc2);
    G_matr(nBlocs)=matr((ii-1)*Nloc2+(1:Nloc2));  
    G_ran(nBlocs)=ran+(ii-1)*Nloc;
    G_col(nBlocs)=col+(ii-1)*Nloc;     
end




	 
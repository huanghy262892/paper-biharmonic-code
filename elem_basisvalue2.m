function [vbasis,dbasis,Lbasis,dLbasis]=elem_basisvalue2(K,xhat, yhat,inv_BE)
% vbasis-->varphi(lambda_2,\lambda_3)
%dbasis-->nabla(varphi) gradient
%Lbasis-->Delta(varphi) Laplace operation
%dLbasis-->gradient(Delta(varphi))
Nloc=(K+2)*(K+1)/2; 
B=inv_BE';
[ssn,ttn,dssn,dttn,ddssn,ddttn,d3ssn,d3ttn]=elem_basis2(K,xhat,yhat);
ii = 1;
vbasis=zeros(1,Nloc);
dbasis=zeros(1,2*Nloc);
Lbasis=vbasis;
dLbasis=dbasis; 
tt=zeros(2,1);
for i=1:(K+1)
    for j=1:i
		vbasis(ii) = ssn(i-j+1) * ttn(j);
		ii = ii + 1;
    end
end
ii = 0;
for i=1:(K+1)
    for j=1:i 
		tt(1) = dssn(i - j+1) * ttn(j);
		tt(2) = ssn(i - j+1) * dttn(j);	
        tt=B*tt;
        ii=ii+1;
        dbasis(ii)=tt(1);
        dbasis(ii+Nloc)=tt(2);
    end
end
%Laplace basis functions---Delta basis function
H=zeros(2,2);

ii=0;
for i=1:(K+1) %多项式次数i-1
    for j=1:i
        H(1,1)=ddssn(i-j+1) * ttn(j);
        H(1,2)= dssn(i-j+1) * dttn(j);
        H(2,1)=H(1,2);
        H(2,2)=ssn(i-j+1) * ddttn(j); 
        ii=ii+1;
        Lbasis(ii)=B(1,:)*H*B(1,:)'+B(2,:)*H*B(2,:)';
    end
end 
ii=0;
for i=1:(K+1) %多项式次数i-1
    for j=1:i
        H(1,1)=d3ssn(i-j+1) * ttn(j);
        H(1,2)= ddssn(i-j+1) * dttn(j);
        H(2,1)=H(1,2);
        H(2,2)=dssn(i-j+1) * ddttn(j); 
        ii=ii+1;
        tt(1)=B(1,:)*H*B(1,:)'+B(2,:)*H*B(2,:)';
        H(1,1)=ddssn(i-j+1) * dttn(j);
        H(1,2)= dssn(i-j+1) * ddttn(j);
        H(2,1)=H(1,2);
        H(2,2)=ssn(i-j+1) * d3ttn(j);  
        tt(2)=B(1,:)*H*B(1,:)'+B(2,:)*H*B(2,:)';
        tt=B*tt;
        dLbasis(ii)=tt(1);
        dLbasis(ii+Nloc)=tt(2);
    end
end
        
               

 
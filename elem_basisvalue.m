function [vbasis,dbasis]=elem_basisvalue(K,xhat, yhat,inv_BE)
Nloc=(K+2)*(K+1)/2; 
 [ssn,ttn,dssn,dttn]=elem_basis(K,xhat,yhat);
ii = 1;
vbasis=zeros(1,Nloc);
dbasis=zeros(1,2*Nloc);
for i=1:(K+1)
    for j=1:i
		vbasis(ii) = ssn(i - j+1) * ttn(j);
		ii = ii + 1;
    end
end
ii = 1;
for i=1:(K+1)
    for j=1:i 
		dbasis(ii) = dssn(i - j+1) * ttn(j);
		dbasis(Nloc+ ii) = ssn(i - j+1) * dttn(j);		
		temp = dbasis(ii)*inv_BE(1,1) + dbasis(Nloc+ii)*inv_BE(2,1);
		dbasis(Nloc+ii) = dbasis(ii)*inv_BE(1,2) + dbasis(Nloc+ii)*inv_BE(2,2);
		dbasis(ii)= temp;
		ii = ii + 1;
    end
end
	 

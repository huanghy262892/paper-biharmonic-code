function [ssn,ttn,dssn,dttn]=elem_basis(K,xhat,yhat)
ssn=zeros(1,K+1);
ttn=ssn;
dssn=ssn;
dttn=ssn;
ssn(1) = 1.0;
ttn(1) = 1.0;
dssn(1) = 0.0;
dttn(1) = 0.0;	 
ssn(2) = xhat;
ttn(2) = yhat;
dssn(2) = 1.0;
dttn(2) = 1.0;	 
for i=3:(K+1) 
	ssn(i) = xhat*ssn(i - 1);
	ttn(i) = yhat*ttn(i - 1);
	dssn(i) = xhat*dssn(i - 1);
	dttn(i) = yhat*dttn(i - 1);		 
end
for i=3:(K+1) 
	dssn(i) = dssn(i) * (i-1);
	dttn(i) = dttn(i) * (i-1);		 
end
 
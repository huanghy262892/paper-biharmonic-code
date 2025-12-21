function [ssn,ttn,dssn,dttn,ddssn,ddttn,d3ssn,d3ttn]=elem_basis2(K,xhat,yhat)
%K>=1
ssn=zeros(1,K+1);
ttn=ssn;
dssn=ssn;
dttn=ssn;
ddssn=ssn;
ddttn=ssn;
d3ttn=zeros(1,K+1);
d3ssn=zeros(1,K+1);
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
if K>=2
    ddssn(3)=1.0;
    ddttn(3)=1.0;
    for i=4:(K+1)
        ddssn(i) = xhat*ddssn(i - 1);
	    ddttn(i) = yhat*ddttn(i - 1);
    end
end
if K>=3
    d3ssn(4)=1;
    d3ttn(4)=1;
    for i=5:(K+1)
        d3ssn(i)=d3ssn(i-1)*xhat;
        d3ttn(i)=d3ttn(i-1)*yhat;
    end    
end
for i=3:(K+1) 
	dssn(i) = dssn(i) * (i-1);
	dttn(i) = dttn(i) * (i-1);
    ddssn(i)=ddssn(i)*(i-1)*(i-2);
    ddttn(i)=ddttn(i)*(i-1)*(i-2);
    d3ssn(i)=d3ssn(i)*(i-1)*(i-2)*(i-3);
    d3ttn(i)=d3ttn(i)*(i-1)*(i-2)*(i-3);
end

 
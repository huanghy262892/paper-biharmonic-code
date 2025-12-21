function xy=coordinate(AX,BX,CY,DY,NX,NY)  
hx = (BX - AX) / NX;
hy = (DY - CY) / NY;
xy=zeros((NX+1)*(NY+1),2);
for ii=0:NY
    for jj=0:NX
		ss = ii*(NX+1) + jj+1;
		xy(ss,1) = AX + hx*jj;
		xy(ss,2) = CY + hy*ii;
    end
end

	 
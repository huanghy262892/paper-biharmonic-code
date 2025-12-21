function elem= element(NX,NY,EN) 
elem=zeros(EN,3);
for ii=1:NY
	for jj=1:NX 
		ss = (ii-1) * 2 * NX+(jj-1)*2+1;
		elem(ss,1) = (ii-1)*(NX + 1) + jj;
		elem(ss,2) = (ii-1)*(NX + 1) + jj + 1;
		elem(ss,3) = ii*(NX + 1) + jj;
		elem(ss + 1,1) = (ii-1)*(NX + 1) + jj + 1;
		elem(ss + 1,2) = ii*(NX + 1) + jj + 1;
		elem(ss + 1,3) = ii*(NX + 1) + jj;
    end
end
		 

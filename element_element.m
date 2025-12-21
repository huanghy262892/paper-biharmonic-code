function ete=element_element(EN,NDOFS,elem)
sk=1;
a1=zeros(1,4);
efn=zeros(NDOFS,2);
ete=zeros(EN,3);
etf=zeros(EN,3);
for ii=1:EN
    for jj=1:3
		a1(jj) = elem(ii,jj);
    end
	a1(4) = elem(ii,1);
	for jj=1:3
		if (a1(jj)>a1(jj + 1)) 
				efn(sk,1) = a1(jj+1);
				efn(sk,2) = a1(jj);
        else  
				efn(sk,1) = a1(jj);
				efn(sk,2) = a1(jj + 1);
        end
		sk = sk + 1;
    end
end
for ii=1:EN
    for jj=1:3 
		m1 = (ii-1) * 3 + jj;
		kk = 1;
		ete(ii,jj) = ii;
		etf(ii,jj) = jj;
		for ss=1:EN
            if (kk ~= 2)  
				if (ss ~= ii)	 
					for tt = 1:3  
						m2 = (ss-1) * 3 + tt;
						if ((efn(m2,1) == efn(m1,1)) && (efn(m2,2) == efn(m1,2))) 
							ete(ii,jj) = ss;
							etf(ii,jj) = tt;
							kk = kk + 1;
                        end
                    end
                end
            end
        end
    end 
end


function [ss1,ss2]=loc_coor_quad(index_E1,xy1,BE1, xy2,inv_BE2, xg)
% ²Î¿¼×ø±ê(lambda_2,lambda_3)
xy=zeros(1,2);
ss1=zeros(1,2);
ss2=zeros(1,2);
if (index_E1 == 1) 
	ss1(1) = xg;
	ss1(2) = 0.0;
elseif (index_E1 == 2) 
	ss1(1) =xg  ;
	ss1(2) =1.0-xg;
else 
	ss1(1) = 0.0;
	ss1(2) = xg;
end
for ii = 1:2
	xy(ii) =xy1(1,ii)-xy2(1,ii);
    for jj=1:2
        xy(ii) = xy(ii) + BE1(ii,jj) * ss1(jj);
    end
end 
for ii=1:2
    ss2(ii) = 0.0;
	for jj=1:2
        ss2(ii) = ss2(ii) + inv_BE2(ii,jj) * xy(jj);
    end
end
 
function [tk,Ak]=Gauss_nodes_weights_1D(n)

syms x
p=sym2poly(diff((x^2-1)^(n+1),n+1))/(2^n*factorial(n));
tk=roots(p);
 
% 计算求积系数
tol=1e-16;
Ak=zeros(n+1,1);
for i=1:n+1
    xkt=tk;
    xkt(i)=[];
    pn=poly(xkt);
    fp=@(x)polyval(pn,x)/polyval(pn,tk(i));
    Ak(i)=integral(fp,-1,1,'AbsTol',1e-20,'RelTol',1e-20); % 求积系数
end

% for i=1:(n+1)
%     fprintf('%.20e\n',tk(i));
% end
 

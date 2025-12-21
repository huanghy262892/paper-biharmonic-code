function [ran,col]=generate_submatrix_index(Nloc)
Nloc2=Nloc*Nloc;
ran=zeros(1,Nloc2);
col=ran;
for ii=1:Nloc
    ss=(ii-1)*Nloc+(1:Nloc);
    ran(ss)=ii*ones(1,Nloc);
    col(ss)=1:Nloc;
end
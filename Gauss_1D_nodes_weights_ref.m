function [g,gweights]= Gauss_1D_nodes_weights_ref(n)

g=zeros(1,n);

if n==1

    g(1)=0;

elseif n==2

    g(1)=-0.57735026918962576451;
    g(2)=0.57735026918962576451;

elseif n==3

    g(1)=-0.77459666924148337704;
    g(2)=0;
    g(3)=0.77459666924148337704;

elseif n==4

    g(1)=-sqrt(3+2*sqrt(6/5))/sqrt(7);
    g(2)=-sqrt(3-2*sqrt(6/5))/sqrt(7);
    g(3)=-g(2);
    g(4)=-g(1);

elseif n==5

    g(1)=-0.9061798459386639928;
    g(2)=-0.53846931010568309104;
    g(3)=0;
    g(4)=-g(2);
    g(5)=-g(1);

elseif n == 6
    g(1) = -0.93246951420315202781;
    g(2) = -0.66120938646626451366;
    g(3) = -0.23861918608319690863;
    g(4) = 0.23861918608319690863;
    g(5) = -g(2);
    g(6) = -g(1);
end

gweights=zeros(1, n);

if n == 1

    gweights = 2;

elseif n == 2

    gweights=[1,1];

elseif n == 3

    gweights=[5/9 8/9 5/9];

elseif n == 4

    gweights=[0.347854845137454,0.652145154862546,0.652145154862546,0.347854845137454];

elseif n == 5

    gweights(1)=0.2369268850561890875;
    gweights(2)=0.478628670499366468;
    gweights(3)=0.56888888888888888889;
    gweights(4)=gweights(2);
    gweights(5)=gweights(1);

elseif n == 6
    gweights(1) = 0.17132449237917034504;
    gweights(2) = 0.36076157304813860757;
    gweights(3) = 0.46791393457269104739;
    gweights(4) = gweights(3);
    gweights(5) = gweights(2);
    gweights(6) = gweights(1);
else
    disp('this program works on for n_gnodes <= 6')
    disp('n_gnodes is too large')
    disp('terminate the program by ctr+c')
end
 


return;



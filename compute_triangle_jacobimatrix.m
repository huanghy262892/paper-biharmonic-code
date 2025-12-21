function [BE1,inv_BE1,det1]=compute_triangle_jacobimatrix(xy)
%计算一般三角形的直角坐标(x,y)到标准三角形的坐标(\lambda_2,\lambda_3)
%BE1为标准坐标变换到直角坐标的Jacobi矩阵
%inv_BE1为直角坐标变换到标准坐标的Jacobi矩阵
BE1(1,1) = xy(2,1) - xy(1,1); 
BE1(1,2) = xy(3,1) - xy(1,1); 
BE1(2,1) = xy(2,2) - xy(1,2); 
BE1(2,2) = xy(3,2) - xy(1,2); 
det1 = BE1(1,1) * BE1(2,2) - BE1(1,2) * BE1(2,1);     
inv_BE1(1,1) = BE1(2,2) / det1;
inv_BE1(1,2) = -BE1(1,2) / det1;
inv_BE1(2,1) = -BE1(2,1) / det1;
inv_BE1(2,2) = BE1(1,1) / det1;
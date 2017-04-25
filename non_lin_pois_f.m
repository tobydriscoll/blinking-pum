function [Res,Jac] = non_lin_pois_f(V,M,DX,DY,DXX,DYY,XPf,nu,grid_sq_ind_in,grid_sq_ind_b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[m,~] = size(M);
Res = zeros(m,1);
Jac = zeros(m,length(V));

lap = M(grid_sq_ind_in,:)*(DXX+DYY)*V;

U_in = M(grid_sq_ind_in,:)*V;

Res(grid_sq_ind_in) = -M(grid_sq_ind_in,:)*(DXX+DYY)*V+cosh(M(grid_sq_ind_in,:)*V);

Res(grid_sq_ind_b) = M(grid_sq_ind_b,:)*V;

Jac(grid_sq_ind_in,:) = -M(grid_sq_ind_in,:)*(DXX+DYY)+diag(sinh(M(grid_sq_ind_in,:)*V))*M(grid_sq_ind_in,:);

Jac(grid_sq_ind_b,:) = M(grid_sq_ind_b,:);
end


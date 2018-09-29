domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [41 41];
cheb_struct.cdegs = [11 11];
cheb_struct.split_flag = [true true];

tol_n = [1e-5,1e-4];
parms = [1000,-1,.5,0];

% %Test with 4 patches
Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);
Tree.split(1);
Tree.split(2);

F = PUchebfun(Tree);
F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

nu = 1/200;
f = @(u,leaf)Burgers(u,leaf,nu);
Jac = @(u,leaf) BurgersJacobian(u,leaf,nu);

init = zeros(length(F),1);

% Regular Newton Schwarz Krylov
%
tic;
[sol2,it_hist2, ierr2] = nsoldAS_NK(init,F,f,Jac,[1e-10 1e-9],parms);
toc

% New Schwarz Newton Krylov
%
% tic;
% [sol, it_hist, ierr, x_hist] = nsoldPAR_AS(init,f,Jac,F,[1e-10 1e-9],parms,[1e-7 1e-6]);
% toc

% New two level Schwarz Newton Krylov
%
%
% tic;
% [sol, it_hist3, ierr3, x_hist3] = nsoldPAR_AS_2_level(init,f,Jac,F,[1e-10 1e-10],parms,[1e-7 1e-6],1e-7,2);
% toc


F.sample(sol);
plot(F);
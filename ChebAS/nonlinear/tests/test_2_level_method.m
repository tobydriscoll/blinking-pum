domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [10 10];
cheb_struct.cdegs = [5 5];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);

bound = @(x,y) atan(x.^2+y.^2);

F = PUchebfun(Tree);

tol_n = [1e-10,1e-10];
tol_c = 1e-10;

F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

% f = @(u,leaf) u.^2;
% Jac = @(u,leaf) diag(2*u);
% F.sample(bound);
% init = F.Getvalues();


nu = 20;
bound_f = @(x,y) atan((cos(pi*3/16)*x+sin(pi*3/16)*y)*1);
f = @(u,leaf) Burgers(u,leaf,nu,bound_f);
Jac = @(u,leaf) BurgersJacobian(u,leaf,nu);
F.sample(bound_f);
init = F.Getvalues();


[f0,L,U,p] = SNK_forward_eval(init,F,f,Jac,tol_n);

E = eye(length(F));
JG = zeros(length(F));

for i=1:length(F)
    JG(:,i) = JacobianFowardLU(F,L,U,p,E(:,i));
end

RES = @(sol)SNK_forward_eval(sol,F,f,Jac,tol_n);

AGJ = jacobi(RES,init);

%linear residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [FJv,FJv_hat] = ComputeJacs(F,Jac,init);
% 
% E = eye(length(F));
% JC = zeros(length(F));
% 
% for i=1:length(F)
%     JC(:,i) = ParLinearResidual(E(:,i),F,FJv);
% end
% 
% RES = @(sol) ParResidual(sol,F,f);
% 
% AJC = jacobi(RES,init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init_hat = F.Fine2Coarse(init);
% 
% F.Coarsen();
% 
% E = eye(length(F));
% JC_hat = zeros(length(F));
% 
% for i=1:length(F)
%     JC_hat(:,i) = ParLinearResidual(E(:,i),F,FJv_hat);
% end
% 
% RES = @(sol) ParResidual(sol,F,f);
% 
% AJC_hat = jacobi(RES,init_hat);

F.Refine();
%coarse correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [ r_er,J_v_pls_er ] = CoarseCorrect( F, init,f,Jac,2,tol_c);
% 
% [FJv,FJv_hat] = ComputeJac(F,Jac,init);
% [Lc,Uc,Pc,Qc] = lu(J_v_pls_er);
% 
% E = eye(length(F));
% JC = zeros(length(F));
% 
% for i=1:length(F)
%     JC(:,i) = LinearCoarseCorrect( F, E(:,i),Lc,Uc,Pc,Qc,FJv,FJv_hat,2);
% end
% 
% RES = @(sol) CoarseCorrect(F,sol,f,Jac,2,tol_c);
% 
% AJC = jacobi(RES,init);


%Two level method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [z,L,U,p,J_v_pls_er] = SNK2level_forward_eval(init,F,f,Jac,2,tol_n,tol_c);
% 
% [FJv,FJv_hat] = ComputeJac(F,Jac,init);
% [Lc,Uc,Pc,Qc] = lu(J_v_pls_er);
% 
% E = eye(length(F));
% JC = zeros(length(F));
% 
% for i=1:length(F)
%     JC(:,i) = JacobianFoward2Level(F,L,U,p,Lc,Uc,Pc,Qc,FJv,FJv_hat,E(:,i),2);
% end
% 
% RES = @(sol) SNK2level_forward_eval(sol,F,f,Jac,2,tol_n,tol_c);
% 
% AJC = jacobi(RES,init);






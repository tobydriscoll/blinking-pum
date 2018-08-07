domain = [0 1;0 1];
deg_in = [5 5];
cheb_struct.domain = domain;
cheb_struct.degs = [20 20];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;


Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);

bound = @(x,y) atan(x.^2+y.^2);

F = PUchebfun(Tree);


tol_n = [1e-5,1e-4];

parms = [20,-1,.5,0];

F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

%f = @ SimpNonlinear;
%Jac = @ SimpNonlinearJac;

f = @(u,leaf) CavityFlow(1,u,leaf);
Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);

%F.sample(bound);
%init = F.Getvalues();

num_sols = 3;
init = zeros(num_sols*length(F),1);

% [f0,L,U,p] = ParPreconditionedNewtonForward(init,F,f,Jac);
% 
% E = eye(length(F));
% J = zeros(length(F));
% 
% for i=1:length(F)
%     J(:,i) = JacobianFowardLU(F,L,U,p,E(:,i));
% end
% 
% RES = @(sol) ParPreconditionedNewtonForward(sol,F,f,Jac);
% 
% AJ = jacobi(RES,init);


%Course Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [c_sol,J_v_pls_er ] = CoarseCorrect(F,init,f,Jac);
% 
% T_hat  = CoarseInterfaceInterp(F,1);
% 
% [FJv,FJv_hat] = ComputeJac(F,Jac,init);
% 
% FJv_hat = blkdiag(FJv_hat{:});
% 
% E = eye(length(F));
% J = zeros(length(F));
% 
% for i=1:length(F)
%     J(:,i) = LinearCoarseCorrect( F, E(:,i),J_v_pls_er,T_hat,FJv,FJv_hat);
% end
% 
% RES = @(sol) CoarseCorrect(F,sol,f,Jac);
% 
% AJ = jacobi(RES,init);


% Two level
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f0,L,U,p,J_v_pls_er] = ParPreconditionedTwoLevel(init,F,f,Jac);

T_hat  = CoarseInterfaceInterp(F,num_sols);

[FJv,FJv_hat] = ComputeJac(F,Jac,init);

FJv_hat = blkdiag(FJv_hat{:});

E = eye(num_sols*length(F));
J = zeros(num_sols*length(F));

for i=1:length(F)
    J(:,i) = JacobianFoward2Level(F,L,U,p,J_v_pls_er,T_hat,FJv,FJv_hat,E(:,i));
end

RES = @(sol) ParPreconditionedTwoLevel(sol,F,f,Jac);

AJ = jacobi(RES,init);
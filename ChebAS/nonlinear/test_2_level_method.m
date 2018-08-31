domain = [0 1;0 1];
deg_in = [5 5];
cheb_struct.domain = domain;
cheb_struct.degs = [17 17];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;


Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);

bound = @(x,y) atan(x.^2+y.^2);

F = PUchebfun(Tree);


tol_n = [1e-8,1e-7];

parms = [20,-1,.5,0];

F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

f = @ SimpNonlinear;
Jac = @ SimpNonlinearJac;

%f = @(u,leaf) CavityFlow(1,u,leaf);
%Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);

F.sample(@(x,y)sin(x+y));
init = F.Getvalues();

num_sols = 1;
%init = rand(num_sols*length(F),1);
%init = ones(num_sols*length(F),1);

% 
% [f0,L,U,p] = ParPreconditionedNewtonForward(init,F,f,Jac);
% 
% E = eye(length(F));
% JG = zeros(length(F));
% 
% for i=1:length(F)
%     JG(:,i) = JacobianFowardLU(F,L,U,p,E(:,i));
% end
% 
% RES = @(sol) ParPreconditionedNewtonForward(sol,F,f,Jac);
% 
% AGJ = jacobi(RES,init);


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
% JC = zeros(length(F));
% 
% for i=1:length(F)
%     JC(:,i) = LinearCoarseCorrect( F, E(:,i),J_v_pls_er,T_hat,FJv,FJv_hat);
% end
% 
% RES = @(sol) CoarseCorrect(F,sol,f,Jac);
% AJC = jacobi(RES,init);


%Course Correction 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[ r_er, g_mat,g_p_mat,J_er] = CoarseCorrect2(F,init,f,Jac);

% T_hat = CoarseInterfaceInterp(F,num_sols);
% 
% E = eye(length(F));
% JC = zeros(length(F));
% 
% for i=1:length(F)
%     JC(:,i) = LinearCoarseCorrect2(F,E(:,i),g_mat,g_p_mat,J_er,T_hat);
% end
% 
% RES = @(sol) CoarseCorrect2(F,sol,f,Jac);
% AJC = jacobi(RES,init);


% Two level
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [z,L,U,p,g_mat,g_p_mat,J_er] = ParPreconditionedTwoLevelG(init,F,f,Jac);
% % 
% T_hat  = CoarseInterfaceInterp(F,num_sols);
% 
% [FJv,FJv_hat] = ComputeJac(F,Jac,init);
% 
% FJv_hat = blkdiag(FJv_hat{:});
% 
% E = eye(num_sols*length(F));
% J = zeros(num_sols*length(F));
% 
% for i=1:length(F)
%     J(:,i) = JacobianFoward2LevelG(F,L,U,p,g_mat,g_p_mat,J_er,T_hat,E(:,i));
% end
% 
% RES = @(sol) ParPreconditionedTwoLevelG(sol,F,f,Jac);
% 
% AJ = jacobi(RES,init);

%Course Correction 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = ParResidual(init,F,f);
[J,L,U,p] = ComputeJacs(init,F,Jac); 

E = eye(length(F));
JG = zeros(length(F));

for i=1:length(F)
    JG(:,i) = JacobianFoward(F,J,E(:,i));
end

RES = @(sol) ParResidual(sol,F,f);

AGJ = jacobi(RES,init);





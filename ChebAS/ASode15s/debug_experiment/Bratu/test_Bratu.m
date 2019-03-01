domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [10 10];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

R = 1;

%Test with 2 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
% Tree.split(2);
 
U_t = PUchebfun(Tree);
U_t.sample(@(x,y) zeros(size(x)));

V_t = U_t.copy();

setInterpMatrices(U_t,false);
setInterpMatrices(V_t,false);

[BratuOp,M,y0] = SetBratu(U_t,R,false);

rhs = zeros(size(y0));

Burger_step = @(y,leaf)SteadyStateBurgers(y,leaf,R);
Burger_step_jac = @(y,leaf)SteadyStateBurgersJacobian(y,leaf,R);

[z,L,U,p] = SNK_time_deriv_resid(0,y0,rhs,{U_t},BratuOp,1,M);

AJ = jacobi(@(z)SNK_time_deriv_resid(1,z,rhs,{U_t},BratuOp,1,M),y0);

E = eye(length(y0));

AJ3 = E;

for i=1:length(y0)
    AJ3(:,i) = JacobianFowardLUTime({U_t},L,U,p,E(:,i));
end
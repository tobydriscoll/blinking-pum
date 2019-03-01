domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

R = 80;

odetol = 1e-4;
tspan = [0 1];

%Test with 2 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 Tree.split(2);
 
U_t = PUchebfun(Tree);
U_t.sample(@(x,y) zeros(size(x)));

V_t = U_t.copy();

setInterpMatrices(U_t,false);
setInterpMatrices(V_t,false);

[BurgerOp,M,y0] = SetBurgers(U_t,V_t,R,true);

dt = 0.1;
% pred = y0+dt*ParLocalResidual(0,y0,{U_t,V_t},BurgerOp);
% 
% for i=1:10
%     
% [z,L,U,p] = SNK_time_deriv_resid(dt,pred,y0,{U_t,V_t},BurgerOp,dt,M);
% 
% del = JacobianFowardLUTime({U_t,V_t},L,U,p,-z);
% 
% pred = pred + del;
% 
% [norm(z) norm(del)]
% 
% end
% E = eye(length(y0));
% 
% J = E;
% 
% for i=1:length(y0)
% J(:,i) = JacobianFowardLUTime({U_t,V_t},L,U,p,E(:,i));
% end
% 
% AJ = jacobi(@(z)SNK_time_deriv_resid(0.1,z,y0,{U_t,V_t},BurgerOp,0.1,M),yp0);

opt = odeset('mass',M,'reltol',odetol,'abstol',odetol);
[t,U] = ASode15s(BurgerOp,tspan,y0,{U_t,V_t},opt);
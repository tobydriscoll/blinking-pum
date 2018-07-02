domain = [0 1;0 1];
split_flag = [false false];
degs = [10 10];
cdegs = [9 9];
tol = 1e-8;
R = 80;
num_sols = 2;

f = @(x,y,t) 3/4 - 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
g = @(x,y,t) 3/4 + 1./(4*(1+exp((R*(4*y-4*x-t)/32))));

cheb_struct.domain = domain;
cheb_struct.degs = degs;
cheb_struct.split_flag = split_flag;
cheb_struct.cdegs = cdegs;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);

F = PUchebfun(Tree);
%set mass matricies
for i=1:length(F.leafArray)
   M{i} = BurgersMassMatrix(F.leafArray{i});
end


F.Setvalues(@(x,y)f(x,y,0));
u0 = F.Getvalues();
F.Setvalues(@(x,y)g(x,y,0));
v0 = F.Getvalues();

y0 = [u0;v0];

%set interpolation matrices
setInterpMatrices(F);

odetol = 1e-3;

tspan = [0 0.5];

dh = 0.01;

yp0 = ParLocalResidual(0,y0,F,@(Approx,t,y) BurgersEvaluation(Approx,t,y,R),2);
pred = y0+dh*yp0;
ynew = pred;

for j=1:7
    
    dif1 = reshape(pred-ynew,length(F),num_sols);
    
    RHS = [];
    for i=1:length(F.leafArray)
        tmp = dif1((i-1)*prod(degs)+(1:prod(degs)),:);
        RHS = [RHS;reshape(M{i}*tmp(:),prod(degs),2)];
    end
    
    RHS = RHS(:);
    
    [z,J] = ParPreconditionedNewtonForwardTime(0+dh,ynew,RHS,F,@(Approx,t,y) BurgersEvaluation(Approx,t,y,R),2,dh,@(t,y,approx)BurgersJacobian(t,y,approx,R),M);
    
    for i=1:length(F.leafArray)
        [L{i},U{i},p{i}] = lu(J{i},'vector');
    end

    [del,~,~,~,~] = gmres(@(x)JacobianFowardLU(F,L,U,p,x,2),-z,[],[]);
    
    ynew = ynew+del;
end



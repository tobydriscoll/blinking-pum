domain = [0 1;0 1];
split_flag = [false false];
degs = [32 32];
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

dh = 0.05;

yp0 = ParLocalResidual(0,y0,F,@(Approx,t,y) BurgersEvaluation(Approx,t,y,R),2);
pred = y0+dh*yp0;
ynew = pred;

<<<<<<< HEAD
tic;
for j=1:5
=======
del_old = zeros(size(yp0));

for j=1:4
>>>>>>> parent of 46f15b1... updated burger_be_step with single approximation comparison
    
    dif1 = reshape(pred-ynew,length(F),num_sols);
    
    RHS = [];
    for i=1:length(F.leafArray)
        tmp = dif1((i-1)*prod(degs)+(1:prod(degs)),:);
        RHS = [RHS;reshape(M{i}*tmp(:),prod(degs),2)];
    end
    
    RHS = RHS(:);
    
    [z,J] = ParPreconditionedNewtonForwardTime(0+dh,ynew,RHS,F,@(Approx,t,y) BurgersEvaluation(Approx,t,y,R),2,dh,@(t,y,approx)BurgersJacobian(t,y,approx,R),M);
    
    
    [del,~,~,~,~] = gmres(@(x)JacobianFoward(F,J,x,2),-z,[],[]);
    
    ynew = ynew+del;
    
<<<<<<< HEAD
<<<<<<< HEAD
    d1(j) = norm(del);
=======
    d1(j) = norm(del,inf);
>>>>>>> e25f1768842eb92dcb2072b3176545ae3a612657
    
    norm(del,inf)
   
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

<<<<<<< HEAD
domain = [0 1;0 1];
split_flag = [false false];
degs = [64 64];
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


F = PUchebfun(Tree);

%set mass matricies
for i=1:length(F.leafArray)
   M{i} = BurgersMassMatrix(F.leafArray{i});
end


odetol = 1e-3;

tspan = [0 0.5];

dh = 0.05;

yp0  = BurgersEvaluation(Tree,0,y0,R);
pred = y0+dh*yp0;
ynew = pred;

 J = M{1} - dh * BurgersJacobian(dh,pred,Tree,R);
 
tic; 
for j=1:5
    
    dif1 = pred-ynew;
    
    rhs = dh*BurgersEvaluation(Tree,dh,ynew,R) + M{1}*dif1;
    
    J = dh * BurgersJacobian(dh,ynew,Tree,R)-M{1};

             
    
    del = J\(rhs);
    
    ynew = ynew-del;
    
    d2(j) = norm(del);
    
    norm(del)
=======
    d(j) = norm(del);
>>>>>>> parent of 46f15b1... updated burger_be_step with single approximation comparison
   
end
toc
=======
% domain = [0 1;0 1];
% split_flag = [false false];

% degs = [40 40];

% cdegs = [9 9];
% tol = 1e-8;
% R = 80;
% num_sols = 2;
% 
% f = @(x,y,t) 3/4 - 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
% g = @(x,y,t) 3/4 + 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
% 
% cheb_struct.domain = domain;
% cheb_struct.degs = degs;
% cheb_struct.split_flag = split_flag;
% cheb_struct.cdegs = cdegs;
% cheb_struct.tol = tol;
% 
% Tree = ChebPatch(cheb_struct);
% 
% F = PUchebfun(Tree);
% %set mass matricies
% for i=1:length(F.leafArray)
%    M{i} = BurgersMassMatrix(F.leafArray{i});
% end
% 
% 
% F.Setvalues(@(x,y)f(x,y,0));
% u0 = F.Getvalues();
% F.Setvalues(@(x,y)g(x,y,0));
% v0 = F.Getvalues();
% 
% y0 = [u0;v0];
% 
% %set interpolation matrices
% setInterpMatrices(F);
% 
% 
% F = PUchebfun(Tree);
% 
% %set mass matricies
% for i=1:length(F.leafArray)
%    M{i} = BurgersMassMatrix(F.leafArray{i});
% end
% 
% 
% odetol = 1e-3;
% 
% tspan = [0 0.5];
% 
% dh = 0.05;
% 
% yp0  = BurgersEvaluation(Tree,0,y0,R);
% pred = y0+dh*yp0;
% ynew = pred;

% 
%  J = M{1} - dh * BurgersJacobian(dh,pred,Tree,R);

% for j=1:5
%     
%     dif1 = pred-ynew;
%     
%     rhs = dh*BurgersEvaluation(Tree,dh,ynew,R) + M{1}*dif1;
%     
%     J = dh * BurgersJacobian(dh,ynew,Tree,R)-M{1};
% 
%              
%     
%     del = J\(rhs);
%     
%     ynew = ynew-del;
%     
%     d2(j) = norm(del);
%     
%     norm(del)
%    
% end
>>>>>>> 73a000f
>>>>>>> e25f1768842eb92dcb2072b3176545ae3a612657



function [Blinks,M,y0,result,initial_H] = setBlinks(H_tree,P_tree,pctClosed,BoundaryH,pA,pS,he)
%This function sets up the blink objects for each leaf. Here 'blink' is set
%to the NonlinOp property. adsfasdf

degs = [65 65];

result = blink(pctClosed,degs,[-1 1;-1 1],BoundaryH);

result.pA = pA;
result.pS = pS;



result.n = degs;
result.boundaryH = 13;
result.percentClosed = pctClosed;
result.odetol = 1e-4;
result.h_e = he;   
result.initcond = 'laplace';

domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = degs;
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

initial_H = ChebPatch(cheb_struct);
initial_P = ChebPatch(cheb_struct);

[H,P] = result.initial;

initial_H.sample(H(:));
initial_P.sample(P(:));

for i=1:length(H_tree.leafArray)
    
    %Set blink object and blink motion
    Blinks{i} = blink(pctClosed,H_tree.leafArray{i}.degs,H_tree.leafArray{i}.domain,BoundaryH);
    Blinks{i}.percentClosed = pctClosed;
    
    Blinks{i}.pA = pA;
    Blinks{i}.pS = pS;
    Blinks{i}.h_e = he;
    
    G = H_tree.leafArray{i}.leafGrids();
    
    H = initial_H.evalfGrid(G);
    
    H_tree.leafArray{i}.Setvalues(H(:));
    H_tree.leafArray{i}.sample(H(:));
    
    P = initial_P.evalfGrid(G);
    
    P_tree.leafArray{i}.Setvalues(P(:));
    P_tree.leafArray{i}.sample(P(:));
    
    outer_boundary = H_tree.leafArray{i}.outer_boundary;
    
    Blinks{i}.disc.boundary.loc_outer = outer_boundary;  
        
end

    H_tree.pack();
    
    for i=1:length(H_tree.leafArray)      
        
        Blinks{i}.disc.num.h = length(H_tree.leafArray{i});
        
        %Set mass matrix
        M{i} = Blinks{i}.massMatrix();
    end

    y0 = [H_tree.Getvalues();P_tree.Getvalues()];
    
end


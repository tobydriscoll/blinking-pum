function [Blinks,M,y0] = setBlinks(H_tree,P_tree,param,initstate)
%This function sets up the blink objects for each leaf. Here 'blink' is set
%to the NonlinOp property. adsfasdf

initial_H = ChebPatch(param);
initial_P = ChebPatch(param);

if nargin < 4
    bobj = blink(param.percentClosed,param.degs,param.domain,param.BoundaryH,param.fluxvolume);

    bobj.n = param.degs;
    bobj.h_e = param.he;
    bobj.initvolume = param.initvolume;
    bobj.initcond = 'laplace';
    [H,P] = bobj.initial;
    initial_H.sample(H(:));
    initial_P.sample(P(:));
else
    initial_H.sample(initstate.H);
    initial_P.sample(initstate.P);
end

for i=1:length(H_tree.leafArray)
    
    %Set blink object and blink motion
    Blinks{i} = blink(param.percentClosed,H_tree.leafArray{i}.degs,H_tree.leafArray{i}.domain,param.BoundaryH,param.fluxvolume);
    Blinks{i}.percentClosed = param.percentClosed;
    
    Blinks{i}.pA = param.pA;
    Blinks{i}.pS = param.pS;
    Blinks{i}.h_e = param.he;
    Blinks{i}.initvolume = param.initvolume;
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
    M{i} = Blinks{i}.massMatrix();
end

y0 = [H_tree.Getvalues();P_tree.Getvalues()];

end


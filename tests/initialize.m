function [Blinks,M,y0,dy0] = initialize(H,P,param,initstate)
%This function sets up the blink objects for each leaf. Here 'blink' is set
%to the NonlinOp property. adsfasdf

initial_H = ChebPatch(param);
initial_P = ChebPatch(param);

if nargin < 4
    initstate = [];
end

for i=1:length(H.leafArray)    
    %Set blink object and blink motion
    Blinks{i} = blink(param.percentClosed,H.leafArray{i}.degs,H.leafArray{i}.domain,param.BoundaryH,param.fluxvolume);
    Blinks{i}.percentClosed = param.percentClosed;
    Blinks{i}.pA = param.pA;
    Blinks{i}.pS = param.pS;
    Blinks{i}.h_e = param.h_e;
    Blinks{i}.initvolume = param.initvolume;
    
    outer_boundary = H.leafArray{i}.outer_boundary;
    Blinks{i}.disc.boundary.loc_outer = outer_boundary;    
end

if isempty(initstate)
    bobj = blink(param.percentClosed,param.degs,param.domain,param.BoundaryH,param.fluxvolume);
    bobj.n = param.degs;
    bobj.initvolume = param.initvolume;
    bobj.initcond = 'laplace';
    [valH,valP] = bobj.initial;
    initial_H.sample(valH(:));
    initial_P.sample(valP(:));
    dy0 = [];
else   
    dH = copy(H);
    dP = copy(P);
    initial_H.sample(initstate.dH);
    transfer(initial_H,dH)
    initial_P.sample(initstate.dP);
    transfer(initial_P,dP)
    dH.pack();
    dy0 = [ dH.Getvalues(); dP.Getvalues() ];    
    initial_H.sample(initstate.H);
    initial_P.sample(initstate.P);
end

transfer(initial_H,H)
transfer(initial_P,P)
H.pack();
y0 = [H.Getvalues();P.Getvalues()];

for i=1:length(H.leafArray)    
    Blinks{i}.disc.num.h = length(H.leafArray{i});
    M{i} = Blinks{i}.massMatrix();
end


end

function transfer(From,To)
    for i=1:length(To.leafArray)   
        G = To.leafArray{i}.leafGrids();
        val = From.evalfGrid(G);
        To.leafArray{i}.Setvalues(val(:));
        To.leafArray{i}.sample(val(:));
    end
end
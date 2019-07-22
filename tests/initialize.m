function [Blinks,M,y0,dy0] = initialize(H,P,param,initstate)
%This function sets up the blink objects for each leaf. Here 'blink' is set
%to the NonlinOp property. adsfasdf

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
    pi = param;
    pi.degs = [60,60];

    bobj = blink(pi.percentClosed,pi.degs,pi.domain,pi.BoundaryH,pi.fluxvolume);
    bobj.n = pi.degs;
    bobj.initvolume = pi.initvolume;
    bobj.initcond = 'laplace';
    [valH,valP] = bobj.initial;
    H.sample(chebfun2(valH));
    P.sample(chebfun2(valP));
    H.pack();
    y0 = [H.Getvalues();P.Getvalues()];
    dy0 = [];
else   
    dH = copy(H);
    dP = copy(P);
    dH.sample(initstate.dH);
    dP.sample(initstate.dP);
    dH.pack();
    dy0 = [ dH.Getvalues(); dP.Getvalues() ];    
    H.sample(initstate.H);
    P.sample(initstate.P);
    H.pack();
    y0 = [ H.Getvalues(); P.Getvalues() ];    
end


for i=1:length(H.leafArray)    
    Blinks{i}.disc.num.h = length(H.leafArray{i});
    M{i} = Blinks{i}.massMatrix();
end


end

function [h,p,dh,dp] = evaluate(sol,t,H,P)

nh = length(H);
np = length(P);

n = 4*ceil(sqrt(np));
x = chebpts(n);

[ufinal,dufinal] = deval(sol,t);
H.sample(ufinal(1:nh));
P.sample(ufinal(nh+(1:np)));

h = chebfun2(H.evalfGrid({x,x})');
p = chebfun2(P.evalfGrid({x,x})');

if nargout > 2
  dP = copy(P);
  dP.sample(dufinal(nh+(1:np)));
  dp = chebfun2(dP.evalfGrid({x,x})');
  
  dH = copy(H);
  for i=1:length(dH.leafArray)  
      dH.leafArray{i}.SetBoundaryValues(0);
  end

  dH.sample(dufinal(1:nh));
  U = dH.evalfGrid({x,x})';
  %U([1 end],:) = 0;
  %U(:,[1 end]) = 0;
  dh = chebfun2(U);
end
 
end
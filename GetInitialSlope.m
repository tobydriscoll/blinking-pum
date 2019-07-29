function [y,yp] = GetInitialSlope(M,y,yp,t0,PUApproxArray,NonLinOps,reltol,alpha)
%Assume M is diagonal
if ~iscell(PUApproxArray)
	PUApproxArray = {PUApproxArray};
end

num_leaves = length(PUApproxArray{1}.leafArray);
num_sols = length(PUApproxArray);

[ ~,lens_k ] = unpackPUvecs(y,PUApproxArray);

for k=1:num_leaves
	D_k = diag(M{k});
	D_k = mat2cell(D_k,lens_k{k});
	%Account for interfance boundary condition
	for i=1:num_sols
		interface = PUApproxArray{i}.leafArray{k}.inner_boundary;
		D_k{i}(interface) = 0;
	end
	D_k = cell2mat(D_k);
	D{k} = D_k;
	
	locAlg{k} = D_k==0;
end

%get right ordering of indicies
TotalAlg = logical(packPUvecs(locAlg,PUApproxArray));

D = packPUvecs(D,PUApproxArray);

TotalDiff = ~TotalAlg;

f = ParLocalResidual(t0,y,1,PUApproxArray,NonLinOps,alpha);

if norm(f(TotalAlg)) <= 1e-3*reltol*norm(f)
	yp(TotalDiff) = f(TotalDiff) ./ D(TotalDiff);
	return;
end

for iter=1:15
	[J,L,U,p] = ComputeJacsTime(t0,y,PUApproxArray,NonLinOps,1,0,alpha,locAlg);
	
	Md = @(x)ASPreconditionerTime(PUApproxArray,L,U,p,x,TotalAlg,locAlg);
	[delY,~,~,~,gmhist] = gmres(@(x)LinearResidual(PUApproxArray,J,x,alpha,TotalAlg,locAlg),-f(TotalAlg),[],1e-16,200,Md);
	
	res = norm(delY);
	% Weak line search with affine invariant test.
	lambda = 1;
	ynew = y;
	for probe = 1:3
		ynew(TotalAlg) = y(TotalAlg) + lambda*delY;
		
		fnew = ParLocalResidual(t0,ynew,1,PUApproxArray,NonLinOps,alpha);
		
		if norm(fnew(TotalAlg)) <= 1e-3*reltol*norm(fnew)
			y = ynew;
			f = fnew;
			yp(TotalDiff) = fnew(TotalDiff) ./ D(TotalDiff);
			
			return;
		end
		
		resnew = norm(gmres(@(x)LinearResidual(PUApproxArray,J,x,alpha,TotalAlg,locAlg),-fnew(TotalAlg),[],1e-5,200,Md));
		
		if resnew < 0.9*res
			break;
		else
			lambda = 0.5*lambda;
		end
	end
	
	Ynorm = max(norm(y(TotalAlg)),norm(ynew(TotalAlg)));
	
	if Ynorm == 0
		Ynorm = eps;
	end
	
	y = ynew;
	
	f = fnew;
	if resnew <= 1e-3*reltol&Ynorm
		yp(TotalDiff) = f(TotalDiff) ./ D(TotalDiff);
		
		return;
	end
end
end
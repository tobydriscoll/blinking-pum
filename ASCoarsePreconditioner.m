function [ output ] = ASCoarsePreconditioner(PUApprox,domain,rhs)

LEAVES = PUApprox.collectLeaves({});

%initialize solution to zero
PUApprox.sample(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

PUApprox.sample(@(x)zeros(length(x),1));

output = [];



for j=1:4*length(LEAVES)
    
    step_n = 0;
    
    for k=1:length(LEAVES)
        
        dim = LEAVES{k}.degs;
        
        rhs_k = rhsc(step_n+(1:prod(dim)));
        

        
        
       [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
       pointsl = LEAVES{k}.points();
       approx = PUApprox.evalfZone(pointsl(in_border,:));
       rhs_k(in_border) = approx;
        
%         [out_border, in_border] = FindBoundaryGridIndex2D(dim,LEAVES{k}.domain(),domain);
%         grids = LEAVES{k}.leafGrids();
%         for i=1:4
%             if any(in_border{i,1}) && any(in_border{i,2})
%                 approx = PUApprox.evalfZoneGrid({grids{1}(in_border{i,1}) grids{2}(in_border{i,2})});
%                 [X_in,Y_in] = ndgrid(in_border{i,1},in_border{i,2});
%                 IND = X_in & Y_in;
%                 IND = IND(:);
%                 rhs_k(IND) = approx(:);
%             end
%         end
        
        step_n = step_n + prod(dim);
        
        LEAVES{k}.sample((LEAVES{k}.ClinOp\rhs_k));
    end
end

PUApprox.Refine();

output = PUApprox.Getvalues();

end
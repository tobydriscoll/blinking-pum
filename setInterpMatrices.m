% This method sets up the linear operators to be used in the patches for
% the theta method.
%
%   input:
%    Tree: Tree that has the current time step stored
%       L: Operator of PDE u_t = L u
%       B: cell array of boundary condition operator in order NORTH SOUTH EAST WEST
%  border: rhs of boundary condution
%   theta: parameter used in theta method
%    dt,t: time step and current time
%  output:
%     rhs: returns RHS for the theta method.
function setInterpMatrices(PUfun,coarse_too)

if nargin==1
    coarse_too = false;
end

for k=1:length(PUfun.leafArray)
    
    points = PUfun.leafArray{k}.points;
    
    if ~PUfun.ChebRoot.is_leaf
        
                [~,out_border,in_border,~] = FindBoundaryIndex2DSides(PUfun.leafArray{k});
                PUfun.leafArray{k}.Binterp = PUfun.ChebRoot.interpSparseMatrixZone(points(in_border,:));
                
                
    end
    
end

if coarse_too
    
    
    for k=1:length(PUfun.leafArray)
        
        Grid = PUfun.leafArray{k}.leafGrids();
        
        PUfun.leafArray{k}.Coarsen();
        
        PUfun.leafArray{k}.C2FinterpMat = PUfun.leafArray{k}.interpMatrixGrid(Grid);  
        
        PUfun.leafArray{k}.Refine();
    end
    
    PUfun.Coarsen();
   
    for k=1:length(PUfun.leafArray)
    
    points = PUfun.leafArray{k}.points;
    
    if ~PUfun.ChebRoot.is_leaf
        
                %degs = PUfun.leafArray{k}.degs;
                [~,out_border,in_border,~] = FindBoundaryIndex2DSides(PUfun.leafArray{k});
                PUfun.leafArray{k}.CBinterp = PUfun.ChebRoot.interpSparseMatrixZone(points(in_border,:));
                
                %PUfun.leafArray{k}.coarse_outer_boundry = out_border;
                %PUfun.leafArray{k}.coarse_inner_boundry = in_border;
                
    end
    
    
    end

    PUfun.Refine();
    
    
end


end




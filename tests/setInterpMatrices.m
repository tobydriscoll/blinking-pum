function setInterpMatrices(PUfun,coarse_too)

if nargin==1
    coarse_too = false;
end

leaf = PUfun.leafArray;
for k=1:length(leaf)  
    points = leaf{k}.points;    
    if ~PUfun.ChebRoot.is_leaf        
        [~,out_border,in_border,~] = FindBoundaryIndex2DSides(leaf{k});
        leaf{k}.Binterp = PUfun.ChebRoot.interpSparseMatrixZone(points(in_border,:));
    end   
end

if coarse_too
       
    for k=1:length(leaf)
        Grid = leaf{k}.leafGrids();       
        leaf{k}.Coarsen();       
        leaf{k}.C2FinterpMat = leaf{k}.interpMatrixGrid(Grid);     
        leaf{k}.Refine();
    end
    
    PUfun.Coarsen();    
    for k=1:length(leaf)       
        points = leaf{k}.points;
        if ~PUfun.ChebRoot.is_leaf           
            %degs = leaf{k}.degs;
            [~,out_border,in_border,~] = FindBoundaryIndex2DSides(leaf{k});
            leaf{k}.CBinterp = PUfun.ChebRoot.interpSparseMatrixZone(points(in_border,:));           
            %leaf{k}.coarse_outer_boundry = out_border;
            %leaf{k}.coarse_inner_boundry = in_border;
        end        
    end
    
    PUfun.Refine();
end


end




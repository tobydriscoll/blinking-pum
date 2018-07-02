function [out_border,in_border] = FindBoundaryIndex2D(dim,domain,out_domain,points)

out_border = false(dim);
in_border = false(dim);

South = false(dim); South(:,1) = true;

North = false(dim); North(:,end) = true;

East = false(dim); East(1,:) = true;

West  = false(dim); West(end,:) = true;

border = South | North | East | West;

if domain(1,1)==out_domain(1,1)
    out_border = out_border | East;
end

if domain(1,2)==out_domain(1,2)
    out_border = out_border | West;
end

if domain(2,1)==out_domain(2,1)
    out_border = out_border | South;
end

if domain(2,2)==out_domain(2,2)
    out_border = out_border | North;
end

in_border(border) = ~out_border(border);

out_border = out_border(:);

in_border = in_border(:);
end


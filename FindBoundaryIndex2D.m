function [out_border,in_border] = FindBoundaryIndex2D(dim,domain,out_domain,points)

out_border = false(prod(dim),1);
in_border = false(prod(dim),1);

%Here I chop off North and South so I'm not double counting
South = false(prod(dim),1); South(2:dim(1)-1)=true;

North = false(prod(dim),1); North(prod(dim)+2-dim(1):prod(dim)-1)=true;

East  = false(prod(dim),1); East(1:dim(1):prod(dim))=true;

West  = false(prod(dim),1); West(dim(2):dim(1):prod(dim))=true;

if domain(1,1)==out_domain(1,1)
    out_border = out_border | East;
else
    in_border = in_border | East;
end

if domain(1,2)==out_domain(1,2)
    out_border = out_border | West;
else
    in_border = in_border | West;
end

if domain(2,1)==out_domain(2,1)
    out_border = out_border | South;
else
    in_border = in_border | South;
end

if domain(2,2)==out_domain(2,2)
    out_border = out_border | North;
else
    in_border = in_border | North;
end

end


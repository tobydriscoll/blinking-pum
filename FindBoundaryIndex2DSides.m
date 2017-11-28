function [out_border_c,out_border,in_border,in_border_c,in_border_g] = FindBoundaryIndex2DSides(dim,domain,out_domain)

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
out_border_c = cell(4,1);

out_border_c{1} = out_border & East(:);
out_border_c{2} = out_border & West(:) & ~ out_border_c{1};

out_border_c{3} = out_border & South(:) & ~(out_border_c{1} | out_border_c{2});
out_border_c{4} = out_border & North(:) & ~(out_border_c{1} | out_border_c{2} | out_border_c{3});

in_border_c{1} = in_border & East(:);
in_border_c{2} = (in_border & West(:)) & ~ in_border_c{1};

in_border_c{3} = in_border & South(:) & ~(in_border_c{1} | in_border_c{2});
in_border_c{4} = in_border & North(:) & ~(in_border_c{1} | in_border_c{2} | in_border_c{3});


B = reshape(in_border_c{1},dim);
in_border_g{1} = {[true;false(dim(1)-1,1)] B(1,:)};
B = reshape(in_border_c{2},dim);
in_border_g{2} = {[false(dim(1)-1,1);true] B(end,:)};
B = reshape(in_border_c{3},dim);
in_border_g{3} = {B(:,1) [true;false(dim(2)-1,1)]};
B = reshape(in_border_c{4},dim);
in_border_g{4} = {B(:,end) [false(dim(2)-1,1);true]};
end
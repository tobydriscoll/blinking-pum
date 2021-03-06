% This method returns the indicies of the inner boundary and outer
% boundary. Some care is given to not double count indicies.
% 
%  out_border, in_border: logical array for points in the outer boundary
%  and inner boundary.
%
%  out_border_s, in_border_s: cell array for west,east,south,north sides (in
%  that order) for the outer boundary and interface.
function [out_border_s,out_border,in_border,in_border_s,border,border_s] = FindBoundaryIndex2DSides(leaf)

degs = leaf.degs;

domain = leaf.domain;

out_domain = leaf.outerbox;

out_border = false(degs);
in_border = false(degs);

South = false(degs); South(:,1) = true;

North = false(degs); North(:,end) = true;

West = false(degs); West(1,:) = true;

East  = false(degs); East(end,:) = true;

border = South | North | East | West;

border_s{4} = North;
border_s{3} = South & ~ North;
border_s{1} = West & ~(North | South);
border_s{2} = East & ~(North | South + West);



%Might have to make this more robust!
if domain(1,1)==out_domain(1,1)
    out_border = out_border | West;
end

if domain(1,2)==out_domain(1,2)
    out_border = out_border | East;
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
out_border_s = cell(4,1);

out_border_s{1} = out_border & West(:);
out_border_s{2} = out_border & East(:) & ~ out_border_s{1};

out_border_s{3} = out_border & South(:) & ~(out_border_s{1} | out_border_s{2});
out_border_s{4} = out_border & North(:) & ~(out_border_s{1} | out_border_s{2} | out_border_s{3});

in_border_s{1} = in_border & West(:);
in_border_s{2} = (in_border & East(:)) & ~ in_border_s{1};

in_border_s{3} = in_border & South(:) & ~(in_border_s{1} | in_border_s{2});
in_border_s{4} = in_border & North(:) & ~(in_border_s{1} | in_border_s{2} | in_border_s{3});

end
function [out_border,in_border] = FindBoundaryGridIndex2D(dim,domain,out_domain)

no_border = {false(dim(1),1) false(dim(2),1)};

out_border = cell(4,2);

out_border(1,:) = no_border;
out_border(2,:) = no_border;
out_border(3,:) = no_border;
out_border(4,:) = no_border;

in_border = out_border;

east = {[true;false(dim(1)-1,1)] true(dim(2),1)};
west = {[false(dim(1)-1,1);true] true(dim(2),1)};

south = {true(dim(1),1) [true;false(dim(2)-1,1)]};
north = {true(dim(1),1) [false(dim(2)-1,1);true]};

is_out_east = domain(1,1)==out_domain(1,1);
is_out_west = domain(1,2)==out_domain(1,2);

is_out_south = domain(2,1)==out_domain(2,1);
is_out_north = domain(2,2)==out_domain(2,2);


%First go through and figure out the outer boundary.
%The nonsense below is to make sure all outer boundary
%points are in the boundary index. I use logical indexing
%to make my life easier.
if is_out_east
   out_border(1,:) = east;
   
   %Don't double count
   north{1}(1) = false;
   south{1}(1) = false;  
end

if is_out_west
   out_border(2,:) = west;
   
%Don't double count
   north{1}(end) = false;
   south{1}(end) = false;
end

if is_out_south
   out_border(3,:) = south;
   
%Don't double count
   east{2}(1) = false;
   west{2}(1) = false;
end

if is_out_north
   out_border(4,:) = north;
   
   %Don't double count
   east{2}(end) = false;
   west{2}(end) = false;
end

%Now go throught and set the in_border
if ~is_out_east
   in_border(1,:) = east;
   
   %Don't double count
   north{1}(1) = false;
   south{1}(1) = false;
end

if ~is_out_west
   in_border(2,:) = west;
   
%Don't double count
   north{1}(end) = false;
   south{1}(end) = false;
end

if ~is_out_south
   in_border(3,:) = south;
   
end

if ~is_out_north
   in_border(4,:) = north;
end

end


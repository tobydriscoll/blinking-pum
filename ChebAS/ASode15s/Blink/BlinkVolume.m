function [ v ] = BlinkVolume(Blinks,H,t)

num_leaves = length(Blinks);

v = 0;

for i=1:num_leaves
    
    U = reshape(Getvalues(H.leafArray{i},true),H.leafArray{i}.degs);
    
    [xc,wx] = chebpts(H.leafArray{i}.degs(1),H.leafArray{i}.zone(1,:));
    [yc,wy] = chebpts(H.leafArray{i}.degs(2),H.leafArray{i}.zone(2,:));
    
    v = v + volume_fun_subgrid(Blinks{i},t,U,xc,yc,wx,wy);
    
end

end


function [ M ] = BurgersMassMatrix(Approx)
    degs = Approx.degs;
    %Determine outer border
    [~,out_border,in_border,~,border] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    border = out_border;
    M1 = eye(prod(degs));
    M2 = M1;
    M1(border,:) = zeros(sum(border),prod(degs));
    M2(border,:) = zeros(sum(border),prod(degs));
    M = blkdiag(M1,M2);
    %(border,:) = zeros(sum(border),prod(degs));
end
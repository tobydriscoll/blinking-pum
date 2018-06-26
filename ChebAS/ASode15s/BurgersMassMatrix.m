function [ M ] = BurgersMassMatrix(Approx)
    degs = Approx.degs;
    %Determine outer border
    [~,out_border,~,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    M1 = eye(prod(degs));
    M2 = M1;
    M1(out_border,:) = zeros(sum(out_border),prod(degs));
    M2(out_border,:) = zeros(sum(out_border),prod(degs));
    M = blkdiag(M1,M2);
end


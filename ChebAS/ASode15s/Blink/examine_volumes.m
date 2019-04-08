load('blink_00.mat');

vols = zeros(size(U,1),1);
min_H = zeros(size(U,1),1);

for i=1:size(U,1)
    H_sol = U(i,1:length(H))';
    H.sample(H_sol);
    vols(i) = BlinkVolume(Blinks,H,t(i));
    min_H(i) = min(H_sol);
end



function [ cutoff ] = simplify(values,tol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Chebfun's StandardChop algorithm needs at least 17
[n,m] = size(values);
N = max(17, round(n*1.25 + 5));

%We use the coefficients to 
coeffs = chebtech2.vals2coeffs(values);

%
nDiff = N - n;

%Pad the coefficients with zeros
coeffs = [coeffs ; zeros(nDiff, m) ];

coeffs = chebtech2.vals2coeffs(chebtech2.coeffs2vals(coeffs));

tol = max(tol)*ones(1, m);

% Loop through columns to compute CUTOFF.
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff, standardChop(coeffs(:,k), tol(k)));
end

% Take the minimum of CUTOFF and LENGTH(F). This is necessary when padding was
% required.
if cutoff == n-1
    cutoff = n;
end

cutoff = min(cutoff, n);

end


function [WKS] = calcWks(V, F, wks_size, wks_variance, numEV)
% CALC_WKS Computes wks_size wave kernel signatures on M with n eigenfunctions
% using wks_variance. 
%
% written by Emanuele Rodolà
M.VERT = V; 
M.TRIV = F; 
M.n = length(V);
M.m = length(F);

M.S_tri = calc_tri_areas(M);
[M.evecs, M.evals, M.S, ~] = cotan_LB(M.VERT, M.TRIV, numEV);

WKS = zeros(M.n, wks_size);

log_E = log(max(M.evals(1:numEV), 1e-6))';
e = linspace(log_E(2), (max(log_E))/1.02, wks_size);
sigma = (e(2)-e(1)) * wks_variance;

C = zeros(1, wks_size);
for i = 1:wks_size
    WKS(:,i) = sum(...
        (M.evecs(:,1:numEV)).^2 .* repmat( exp((-(e(i) - log_E).^2) ./ (2*sigma.^2)),M.n,1), 2);
    C(i) = sum(exp((-(e(i)-log_E).^2)/(2*sigma.^2)));
end

WKS(:,:) = WKS(:,:)./repmat(C,M.n,1);

end


function S_tri = calc_tri_areas(M)
%CALC_TRI_AREAS Calculates the area of each triangle in M.
%   M needs to have fields VERT, TRIV, n and m.
% 
%   copyright (c) 2016 Emanuele Rodolà 

S_tri = zeros(M.n,1);

for k=1:M.m
    e1 = M.VERT(M.TRIV(k,3),:) - M.VERT(M.TRIV(k,1),:);
    e2 = M.VERT(M.TRIV(k,2),:) - M.VERT(M.TRIV(k,1),:);
    S_tri(k) = 0.5*norm(cross(e1,e2));
end

end





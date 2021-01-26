function ci = bca(estimate,boots,jacks,thresh)
% bca: Bootstrap bias-corrected accelerated (BCa) intervals, ripped from bootci
%
% This is necessary to maintain bootstrap and jackknife balanced sampling for
%    the separate hyperspheres.

if ~exist('thresh','var') || isempty(thresh), thresh = 0.05; end

% this is the bias correction
z_0 = fz0(boots,estimate);

N = size(jacks,1);
weights = repmat(1/N,N,1);

% acceleration finding, see DiCiccio and Efron (1996)
mjstat = mean(jacks,1);
score = bsxfun(@minus,mjstat,jacks); % score function at stat; ignore (N-1) factor because it cancels out in the skew
iszer = all(score==0,1);
skew = mean(score.^3,1) ./ (mean(score.^2,1).^1.5) /sqrt(N); % skewness of the score function
skew(iszer) = 0;
acc = skew/6;  % acceleration

% transform back with bias corrected and acceleration
z_alpha1 = norminv(thresh/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z_0 +(z_0+z_alpha1)./(1-acc.*(z_0+z_alpha1)));
pct1(z_0==Inf) = 100;
pct1(z_0==-Inf) = 0;
pct2 = 100*normcdf(z_0 +(z_0+z_alpha2)./(1-acc.*(z_0+z_alpha2)));
pct2(z_0==Inf) = 100;
pct2(z_0==-Inf) = 0;

% inverse of ECDF
m = numel(estimate);
lower = zeros(1,m);
upper = zeros(1,m);
for i=1:m
    lower(i) = prctile(boots(:,i),pct2(i),1);
    upper(i) = prctile(boots(:,i),pct1(i),1);
end

% return
ci = sort([lower;upper],1);

end % bootbca()


function z0=fz0(bstat,stat)
% Compute bias-correction constant z0
z0 = norminv(mean(bsxfun(@lt,bstat,stat),1) + mean(bsxfun(@eq,bstat,stat),1)/2);
end   % fz0()


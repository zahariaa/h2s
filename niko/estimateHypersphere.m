function [loc,locCI,rad,radCI] = estimateHypersphere(points,nBootstrapSamples)

% uses expected distance between pairs of points to estimate the radius of
% the hypersphere.
% uses the mean (or median) to estimate the location of the center.
% uses bootstrap resampling of the points to estimate 95% confidence
% intervals for both.

%% preparation
[n,d] = size(points);

%% Bootstrap via recursion
if exist('nBootstrapSamples','var') && nBootstrapSamples > 1
   [loc,~,rad] = estimateHypersphere(points); % All-data estimate
   % Bootstrap
   locs = NaN(nBootstrapSamples,d);
   rads = NaN(nBootstrapSamples,1);
   for bootstrapI = 1:nBootstrapSamples
      bsIs = ceil(rand(n,1)*n); % boostrap sample indices
      [locs(bootstrapI,:),~,rads(bootstrapI)] = estimateHypersphere(points(bsIs,:));
   end
   locCI = prctile(locs,[2.5 97.5]);
   radCI = prctile(rads,[2.5 97.5]);
   return
else  locCI = [];   radCI = [];
end

%% estimate location
% optimises L2 cost, corresponds to L1 force equilibrium
loc = mean(points,1);

%% estimate radius via skew-based MVUEs
points = points - repmat(loc,[n 1]);
radii = sqrt(sum(points.^2,2));
dvarrad = sqrt(d)*var(radii/median(radii));
if dvarrad > 0.2 % Assume Gaussian
   rad = std(points(:))*sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
else             % Assume Uniform
   rad = stdPer(d)*std(radii) + median(radii);
end

return
end


%% stdPer
function v = stdPer(d)

expectedStds = [1.2733    1.0115    0.8796    0.8107    0.8384    0.8638    0.9579    1.0403    1.1938  1.4268    1.8384    2.4485];
v = interp1(2.^(1:12),expectedStds,d);
return
end



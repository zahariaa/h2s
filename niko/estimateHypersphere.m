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
   locs = NaN(nBoostrapSamples,d);
   rads = NaN(nBoostrapSamples,1);
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

%% estimate the expected-distance-to-radius ratio
% by interpolation among estimates precomputed in QT_expectedDistBetweenTwoPointsWithinHypersphere
ds= [1 2 3 5 8 15 30 60 100 1000];
expectedDists = [0.6582    0.9073    1.0261    1.1490    1.2343    1.3181    1.3651    1.3877    1.3993  1.4127];
if d<1000
    expDistPerRad = interp1(ds,expectedDists,d);
else 
    expDistPerRad = sqrt(2);
end

%% estimate radius
dists = pdist(points,'Euclidean');
meanDist = mean(dists);

rad = meanDist/expDistPerRad;

return


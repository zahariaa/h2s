function [loc,locCI,rad,radCI] = estimateHypersphere(points,nBootstrapSamples)

% uses expected distance between pairs of points to estimate the radius of
% the hypersphere.
% uses the mean (or median) to estimate the location of the center.
% uses bootstrap resampling of the points to estimate 95% confidence
% intervals for both.


%% preparation
[n,d] = size(points);

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

%% bootstrap confidence intervals
dists_sq = squareform(dists);

loc_bs = nan(nBootstrapSamples,d);
rad_bs = nan(nBootstrapSamples,1);

for bootstrapI = 1:nBootstrapSamples
    
    % boostrap sample indices
    bsIs = ceil(rand(n,1)*n);
    
    % location estimate
    points_bs = points(bsIs,:);
    loc1 = mean(points_bs,1);   % optimises L2 cost, corresponds to L1 force equilibrium
    loc2 = median(points_bs,1); % optimises L1 cost, corresponds to L0 force equilibrium
    loc_bs(bootstrapI) = mean([loc1; loc2],1);
    
    % radius estimate
    dists_sq_bs = dists_sq(bsIs,bsIs); showRDMs(dists_sq_bs);
    rad_bs(bootstrapI) = mean(squareform(dists_sq_bs))/expDistPerRad;
    
end

locCI = prctile(loc_bs,[2.5 97.5]);
radCI = prctile(loc_bs,[2.5 97.5]);


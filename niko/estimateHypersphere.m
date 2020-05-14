function hyp = estimateHypersphere(points,varargin)
% hyp = estimateHypersphere(points,nBootstrapSamples) % output is a Hypersphere obj
% hyps = estimateHypersphere(points,categories) % output is an array of Hypersphere objs
%
% IF POINTS IS 3D:
% hyp = estimateHypersphere(points,categories) % assumes 3rd dim are bootstraps
% ASSUMES YOU WANT TO COMBINE the 2D points across frames in 3rd dim,
%    using the same categories in categories, but to keep all frames in same space
%    e.g. for a movie:
% hyp = estimateHypersphere(points,categories,1)
%

%% preparation
[n,d,f] = size(points);
nBootstraps= 1;
STRATIFIED = true;     % stratified bootstrap (default)
NORMALIZE  = true;     % pre-normalize points
ESTIMATOR  = 'meandist';
DBUGPLOT   = false;
VERBOSE    = false;

for v = 1:numel(varargin)
   if isa(varargin{v},'Categories')
      categories  = varargin{v};
   elseif isnumeric(varargin{v}) && numel(varargin{v})==1
      nBootstraps = varargin{v};
   elseif ischar(   varargin{v})
      switch( lower(varargin{v}) )
         case 'permute';     STRATIFIED = false;
         case 'stratified';  STRATIFIED = true;
         case 'normalize';   NORMALIZE  = true;
         case 'raw';         NORMALIZE  = false;
         case {'meandist','distance','mcmc','jointml','gaussian','uniformball','uniformcube'}
            ESTIMATOR = lower(varargin{v});
         otherwise warning('bad string input option: %s', varargin{v})
      end
   end
end

% Auto-detect distance matrix (do easy/fast checks first)
if f==1 && ( (n==d && all(diag(points)==0) && issymmetric(points)) ...
   || strcmp(ESTIMATOR,'distance') || d==1 )% this line assumes points are de-squareform'd

   DISTMATRIX = true;
   normtype   = [];    % if normalizing, make full matrix min be 0 and max be 1
   ESTIMATOR  = 'distance'; % override estimator choice
   STRATIFIED = false; % bootstrap must be w/o replacement; stratified doesn't make sense

   if n==d, points = squareform(points); end
else
   DISTMATRIX = false;
   normtype   = 'col'; % if normalizing, adjust each column by overall mean and std
end

if f==1 && NORMALIZE, points = normalize(points,normtype); end

%% recurse to slice/bootstrap points
if exist('categories','var') && numel(categories.labels)>1
   nCats   = numel(categories.labels);

   %% EVALUATE ALL FRAMES IN SAME SPACE
   if f > 1
      points  = reshape(permute(points,[1 3 2]),n*f,d);
      if NORMALIZE
         points = normalize(points,normtype);
      end

      % EVALUATE ALL TOGETHER!
      % Recursive call
      hypEst  = estimateHypersphere(points,categories.internalrepmat(f),'raw',ESTIMATOR);
      % SEPARATE COMBINED HYP
      centers = mat2cell(hypEst.centers,nCats*ones(f,1));
      radii   = num2cell(reshape(hypEst.radii,[nCats f]),1)';
      hyp     = Hypersphere(centers,radii,categories);
      return
   end
   %% permutation or stratified bootstrap test
   if nBootstraps > 1 && ~strcmp(ESTIMATOR,'mcmc')
      catperm = [categories; categories.permute(nBootstraps,STRATIFIED)];
      % Recursive call
      hyp = arrayfun(@(c) estimateHypersphere(points,c,'raw',ESTIMATOR),catperm);
   else % mcmc or single sample (of bootstrap or not)
      if DISTMATRIX
         % build inter-category mean distance matrix
         catix = nchoosek_ix(nCats);
         distsAcross = zeros(1,size(catix,2));
         for i = 1:size(catix,2)
            distsAcross(i) = mean(points( ...
                              categories.select(catix(:,i)).vectorsForDistanceMatrix ));
         end
         % mds on the category-level distance matrix to get centers
         % TODO: what should the dimensionality of the high-dimensional embedding be?
         nDimsEmbedding = max(2,nCats-1);
         centers = mdscale(distsAcross,nDimsEmbedding,'Criterion','metricstress');
         % not doing this anymore: compute centers from all unique distances

         % collect within-category distances, so later recursive calls compute radii
         for i = 1:nCats
            p{i} = points( categories.select(i).vectorsForDistanceMatrix );
         end
      else
         for i = 1:numel(categories.labels)
            ix = categories.select(i).vectors;
            if categories.ispermuted
               ix = undeal(2,@() unique(ix)); ix = ix(2:end);
            end
            if isnumeric(ix) && ~any(ix>1), ix = ~~ix; end
            p{i} = points(ix,:);
         end
      end
      % Recursive call
      hyp = cellfun(@(x) estimateHypersphere(x,'raw',ESTIMATOR,nBootstrapSamples),...
               p,'UniformOutput',false);
      switch(ESTIMATOR)
         case 'mcmc'
            % do some fancy combining
            hyp = cellfun(@unconcat, hyp,'UniformOutput',false); % n-cell of nSamples
            hyp = num2cell(cell2mat_concat(hyp),1); % nSamples-cell of n-Hyperspheres
            hyp = cellfun(@(x) x.concat(categories),hyp); % nSamples n-Hyperspheres
         case 'distance'
            hyp = cell2mat_concat(hyp).concat(categories);
            % replace true centers in all hyp's
            hyp.centers = centers;
         otherwise
            hyp = cell2mat_concat(hyp).concat(categories);
      end
   end
   return
end

%% DO THE ACTUAL ESTIMATION
%% estimate location, re-center points
switch(ESTIMATOR)
   case 'distance' % do nothing
   case 'mcmc'     % do nothing
   otherwise       loc    = mean(points,1); % optimises L2 cost (same as L1 force equilibrium)
                   points = points - repmat(loc,[n 1]);
end

switch(ESTIMATOR)
   case 'distance'
      loc = []; % consistent center points are in a lower call
      % points here are actually distances
      rad = mean(points)/expDistPerRad(d);
   case 'mcmc'
      [loc,rad] = inferHyperspherePosterior(points,nBootstraps,DBUGPLOT,VERBOSE);
   case 'meandist'
      % dists = pdist(points,'Euclidean');
      % cv = cvindex(n,10);
      % loc = mean(cell2mat_concat(cv.crossvalidate(@mean,points)));

      %% estimate radius via skew-based MVUEs
      radii = sqrt(sum(points.^2,2));
      dvarrad = log2(std(radii/median(radii)))+log2(d)/1.5;
      if dvarrad > -0.5 % Assume Gaussian
         rad = mvue_gaussian(points,d,n);
      else              % Assume Uniform
         rad = stdPer(d)*std(radii) + median(radii);
      end
end

%% OUTPUT FORMATTING
if ~exist('categories','var')
   categories = Categories(numel(rad));
end

hyp = Hypersphere(loc,rad,categories);

return
end


%% HELPER FUNCTIONS
function v = stdPer(d)
expectedStds = [1.2733    1.0115    0.8796    0.8107    0.8384    0.8638    0.9579    1.0403    1.1938  1.4268    1.8384    2.4485];
v = interp1(2.^(1:12),expectedStds,d);
return
end


function r = expDistPerRad(d)
% by interpolation among estimates precomputed in QT_expectedDistBetweenTwoPointsWithinHypersphere
expectedDists = [0.9067    1.1037    1.2402    1.3215    1.3660    1.3902    1.4018    1.4081    1.4112  1.4126    1.4134    1.4138];

if d<4097,   r = interp1(2.^(1:12),expectedDists,d);
else         r = sqrt(2);
end
return
end

function rad = mvue_gaussian(points,d,n)
   target    = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
   % rad       = std(points)*target;
   rad       = sqrt(sum(var(points,[],2))/(n-1))*target;
   return
end


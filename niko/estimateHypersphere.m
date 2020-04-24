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
STRATIFIED = true;     % stratified bootstrap (default)
NORMALIZE  = true;     % pre-normalize points
ESTIMATOR  = '';
DBUGPLOT   = false;
VERBOSE    = false;

% Auto-detect distance matrix (do easy/fast checks first)
if f==1 && n==d && all(diag(points)==0) && issymmetric(points)
   DISTMATRIX = true;
   normtype   = [];    % if normalizing, make full matrix min be 0 and max be 1
else
   DISTMATRIX = false;
   normtype   = 'col'; % if normalizing, adjust each column by overall mean and std
end

for v = 2:nargin
   if isa(varargin{v-1},'Categories'), categories        = varargin{v-1};
   elseif  isnumeric(varargin{v-1}),   nBootstrapSamples = varargin{v-1};
   elseif  ischar(   varargin{v-1})
      switch(lower(varargin{v-1}))
         case 'permute';     STRATIFIED = false;
         case 'stratified';  STRATIFIED = true;
         case 'normalize';   NORMALIZE  = true;
         case 'raw';         NORMALIZE  = false;
         otherwise,          ESTIMATOR  = varargin{v-1};
      end
   end
end
if ~exist('nBootstrapSamples','var'),  nBootstrapSamples = 1; end

if f==1 && NORMALIZE, points = normalize(points,normtype); end

%% recurse cell array of points
if exist('categories','var') && numel(categories.labels)>1
%% EVALUATE ALL FRAMES IN SAME SPACE
   if f > 1
      points  = reshape(permute(points,[1 3 2]),n*f,d);
      if NORMALIZE
         points = normalize(points,normtype);
      end
      nCats   = numel(categories.labels);

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
   if nBootstrapSamples > 1 && ~strcmpi(ESTIMATOR,'mcmc')
      catperm = [categories; categories.permute(nBootstrapSamples,STRATIFIED)];
      % Recursive call
      hyp = arrayfun(@(c) estimateHypersphere(points,c,'raw',ESTIMATOR),catperm);
   else
      for i = 1:numel(categories.labels)
         ix = categories.select(i).vectors;
         if categories.ispermuted
            ix = undeal(2,@() unique(ix)); ix = ix(2:end);
         end
         p{i} = points(ix,:);
      end
      % Recursive call
      hyp = cellfun(@(x) estimateHypersphere(x,'raw',ESTIMATOR,nBootstrapSamples),...
               p,'UniformOutput',false);
      if strcmpi(ESTIMATOR,'mcmc')
         % do some fancy combining
         hyp = cellfun(@unconcat, hyp,'UniformOutput',false); % n cell of nSamples
         hyp = num2cell(cell2mat_concat(hyp),1); % nSamples cell of n Hyperspheres
         hyp = cellfun(@(x) x.concat(categories),hyp); % nSamples of Hyperspheres
      else
         hyp = cell2mat_concat(hyp).concat(categories);
      end
   end
   return
end

switch(lower(ESTIMATOR))
   case 'mcmc'
      [loc,rad,llh] = inferHyperspherePosterior(points,nBootstrapSamples,DBUGPLOT,VERBOSE);
   otherwise
      %% estimate location
      loc = mean(points,1); % optimises L2 cost, corresponds to L1 force equilibrium
      % dists = pdist(points,'Euclidean');
      % cv = cvindex(n,10);
      % loc = mean(cell2mat_concat(cv.crossvalidate(@mean,points)));

      %% estimate radius via skew-based MVUEs
      points = points - repmat(loc,[n 1]);
      radii = sqrt(sum(points.^2,2));
      dvarrad = log2(std(radii/median(radii)))+log2(d)/1.5;
      if dvarrad > -0.5 % Assume Gaussian
         rad = sqrt(sum(var(points,[],2))/(n-1))*sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
      else             % Assume Uniform
         rad = stdPer(d)*std(radii) + median(radii);
      end
end

if ~exist('categories','var')
   categories = Categories(numel(rad));
end

hyp = Hypersphere(loc,rad,categories);

return
end



function v = stdPer(d)
expectedStds = [1.2733    1.0115    0.8796    0.8107    0.8384    0.8638    0.9579    1.0403    1.1938  1.4268    1.8384    2.4485];
v = interp1(2.^(1:12),expectedStds,d);
return
end


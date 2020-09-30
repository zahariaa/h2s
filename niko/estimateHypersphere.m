function [hyp,loc_cv] = estimateHypersphere(points,varargin)
% hyp = estimateHypersphere(points,nBootstraps) % output is a Hypersphere obj
% hyps = estimateHypersphere(points,categories) % output is an array of Hypersphere objs
% [hyp,loc_cv] = estimateHypersphere(points) % output is a Hypersphere obj and
%    its 2-fold cross-validated center
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
NORMALIZE  = false;    % pre-normalize points
CVDISTS    = true;     % distances are cross-validated by default
cvdists    = false;    % placeholder if don't have cross-validated distance values
ESTIMATOR  = 'meandist';
DBUGPLOT   = false;
VERBOSE    = false;
INDEPENDENT= false;

for v = 1:numel(varargin)
   if isa(varargin{v},'Categories')
      categories  = varargin{v};
      varargin{v} = [];
   elseif isa(varargin{v},'cvindex')
      cv = varargin{v};
   elseif isnumeric(varargin{v}) && numel(varargin{v})==1
      nBootstraps = varargin{v};
      varargin{v} = [];
   elseif ischar(   varargin{v})
      switch( lower(varargin{v}) )
         case 'permute';     STRATIFIED = false;
         case {'bootstrap','stratified'};  STRATIFIED = true;
         case 'calcstats';   STRATIFIED = [true false];
         case 'normalize';   NORMALIZE  = true;  varargin{v} = [];
         case 'raw';         NORMALIZE  = false; varargin{v} = [];
         case 'nocvdists';     CVDISTS  = false;
         case 'cvdists';       CVDISTS  = true;
         case 'independent'; INDEPENDENT= true;  varargin{v} = [];
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
      if INDEPENDENT % normalize must be false (think about this some more)
         for i = 1:f
            [hyp(i),loc_cv{i}] = estimateHypersphere(points(:,:,i),categories,'raw',varargin{:});
         end
         loc_cv = reshape(cat(1,loc_cv{:}),[nCats f])';
         return
      end
      points  = reshape(permute(points,[1 3 2]),n*f,d);
      if NORMALIZE
         points = normalize(points,normtype);
      end

      % EVALUATE ALL TOGETHER!
      % Recursive call
      [hypEst,loc_cv] = estimateHypersphere(points,categories.internalrepmat(f),...
                                            'raw',nBootstraps,'nocvdists',varargin{:});
      % SEPARATE COMBINED HYP
      centers = mat2cell(hypEst.centers,nCats*ones(f,1));
      radii   = num2cell(reshape(hypEst.radii,[nCats f]),1)';
      loc_cv  = reshape(loc_cv,[nCats f])';
      if CVDISTS, cvdists = cvCenters2cvSqDists(loc_cv); end

      hyp     = Hypersphere(centers,radii,categories,'cvdists',cvdists);
      return
   end

   %% permutation or stratified bootstrap test
   if nBootstraps > 1 && ~strcmp(ESTIMATOR,'mcmc')
      catperm = [categories; categories.permute(nBootstraps,STRATIFIED)];
      % Recursive call
      % [hyp,loc_cv] = arrayfun(@(c) estimateHypersphere(points,c,'raw',varargin{:}),...
      %                              catperm,'UniformOutput',false);
      % hyp = cell2mat_concat(hyp);
      loc_cv = cell(numel(catperm),1);
      if CVDISTS
         parfor iboot = 1:numel(catperm)
            [hyp(iboot),loc_cv{iboot}] = estimateHypersphere(points,catperm(iboot),'raw',varargin{:});
         end
      else
         parfor iboot = 1:numel(catperm)
            hyp(iboot) = estimateHypersphere(points,catperm(iboot),'raw',varargin{:});
         end
      end
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
         
         if CVDISTS
            % Generate cross-validation objects for center bootstraps that makes sure
            %    bootstrap resamplings keep points in each fold independent
            cvs = arrayfun(@(i) cvindex(categories.select(i).vectorsForDistanceMatrix,2),1:nCats);
         end
      else
         [p,ix] = categories.slice(points);
         if CVDISTS, cvs = cellfun(@(i) cvindex(i,2),ix); end
      end

      % Recursive call
      if ~CVDISTS
         hyp = cellfun(@(x) estimateHypersphere(x,'raw',nBootstraps,varargin{:}),...
                            p,'UniformOutput',false);
      else
         [hyp,loc_cv] = cellfun(@(x,cv) estimateHypersphere(x,'raw',nBootstraps,cv,varargin{:}),...
                                     p,num2cell(cvs),'UniformOutput',false);
         % Compute cross-validated distances
         cvdists = cvCenters2cvSqDists(loc_cv);
      end

      switch(ESTIMATOR)
         case 'mcmc'
            % do some fancy combining
            hyp = cellfun(@unconcat, hyp,'UniformOutput',false); % n-cell of nSamples
            hyp = num2cell(cell2mat_concat(hyp),1); % nSamples-cell of n-Hyperspheres
            hyp = cellfun(@(x) x.concat(categories,cvdists),hyp); % nSamples n-Hyperspheres
         case 'distance'
            hyp = cell2mat_concat(hyp).concat(categories,cvdists);
            % replace true centers in all hyp's
            hyp.centers = centers;
         otherwise
            hyp = cell2mat_concat(hyp).concat(categories,cvdists);
      end
   end
   return
end

%% DO THE ACTUAL ESTIMATION
if CVDISTS && ~exist('cv','var')
   cv = cvindex(n,2); % for cross-validated centers estimation
end

%% estimate location, re-center points
switch(ESTIMATOR)
   case {'distance','mcmc','jointml'} % do nothing
   otherwise
      if CVDISTS
         % loc_cv = [mean(points(cv.train(1),:));
         %           mean(points(cv.train(2),:)) ];
         loc_cv = cell2mat_concat(cv.crossvalidate(@mean,points));
      end
      loc    =  mean(points,1); % optimises L2 cost (same as L1 force equilibrium)
      points = points - repmat(loc,[n 1]);
end

switch(ESTIMATOR)
   case 'distance'
      loc = []; loc_cv = [];% consistent center points are in a lower call
      % points here are actually distances
      rad = mean(points)/expDistPerRad(d);
   case 'mcmc'
      [loc,rad,llh] = inferHyperspherePosterior(points,nBootstraps,DBUGPLOT,VERBOSE);

      % just take the best estimate
      ix     = maxix(llh);
      loc    = loc(ix,:);
      rad    = rad(ix);
      
      if CVDISTS
         loc_cv    = cv.crossvalidate(@inferHyperspherePosterior,...
                                               points,nBootstraps,DBUGPLOT,VERBOSE);
         loc_cv    = squeeze(mat2cell(permute(cell2mat_concat(loc_cv,3),[2 3 1]),...
                                      d,2,ones(numel(rad),1)));
         loc_cv = loc_cv{ix};
      end
   case 'jointml'
      tol  = 10^(-7-log2(d));
      opts = optimoptions('fminunc','TolX',tol,'TolFun',tol,'Algorithm','quasi-newton',...
                          'Display','off','GradObj','on');%,'DerivativeCheck','on');
      % search for the center that minimizes the maximum radius, starting at centroid
      loc0      =  mean(points,1);
      objfcn    = @(points) fminunc(@(m) maxRadiusGivenCenter(m,points),loc0,opts);
      [loc,rad] = objfcn(points);
      if CVDISTS, loc_cv = cell2mat_concat(cv.crossvalidate(objfcn,points)); end
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
expectedDists = [0.6673 0.9039 1.1043 1.2407 1.3230 1.3657 1.3898 1.4020 1.4081 1.4111 1.4127 1.4134 1.4138];

if d<4097,   r = interp1(2.^(0:12),expectedDists,d);
else         r = sqrt(2);
end
return
end

%% Objective function to optimize, with gradient
function [fval,grad] = maxRadiusGivenCenter(m,X)

grad = X - repmat(m(:)',[size(X,1) 1]);
fval = sqrt(sum(grad.^2,2));

% Apply max function
[fval,ix] = max(fval);
grad = -grad(ix,:)/fval;
return
end

function rad = mvue_gaussian(points,d,n)
   target    = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
   % rad       = std(points)*target;
   rad       = sqrt(sum(var(points,[],2))/(n-1))*target;
   return
end

function cvdists = cvCenters2cvSqDists(loc_cv)
% Compute SQUARED cross-validated distances from a cell array of 2-fold cross-
%    validated centers (a standard output of estimateHypersphere).
   [f,n] = size(loc_cv);
   if f > 1 % recurse
      for i = 1:f
         cvdists{i} = cvCenters2cvSqDists(loc_cv(i,:));
      end
      return
   end

   cvdists = NaN(nchoosek(n,2),1);
   i=0;
   for a = 1:n-1
      for b = a+1:n, i=i+1;
         cvdists(i) =  (loc_cv{a}(1,:)-loc_cv{b}(1,:))...
                      *(loc_cv{a}(2,:)-loc_cv{b}(2,:))';
      end
   end
   % cvdists = sign(cvdists).*sqrt(abs(cvdists));
   return
end

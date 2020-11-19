function [hyp,loc_cv] = estimateHypersphere(points,varargin)
% estimateHypersphere: the first step in the hypersphere2sphere algorithm,
%    this function takes points as input and outputs a Hypersphere object.
% 
% e.g.
% hyp = estimateHypersphere(points,nBootstraps) % output is a Hypersphere obj
% hyps = estimateHypersphere(points,categories) % output is an array of Hypersphere objs
% [hyp,loc_cv] = estimateHypersphere(points) % output is a Hypersphere obj and
%    % its 2-fold cross-validated center
% hyp = estimateHypersphere(points,'nocvdists') % cross-validated dists not computed
% hyp = estimateHypersphere(points,'calcStats',nBootstraps) % compute both
%    % permutation- and bootstrap-based statistical tests
% 
% Required input:
%    points (2D or 3D array): [n x d x f] array, where n is the number of points/
%       samples, d is the dimensionality of the points, and f (DEFAULT=1) is the
%       numnber of 'frames' in a common space (as with temporal data in a movie)
%       or independent sets of hyperspheres that can be computed in parallel.
%
%% [NOT FULLY IMPLEMENTED] Points can also be a distance matrix--this case will be
%       automatically detected, and a distance-based estimator will be applied.
% 
% Optional inputs:
%    categories (DEFAULT = none): A Categories object. This is perhaps the most
%       critical optional input; it can be used to specify the subsets of the
%       points that denote distinct hyperspheres to estimate, which
%       estimateHypersphere will recursively compute. The output will be a single
%       Hypersphere object, with the Categories object embedded within it. When
%       multiple hyperspheres are given for estimation, the cross-validated
%       distance(s) between them will be computed and also embedded into the
%       Hypersphere object. The cross-validated center coordinates (used to
%       compute the cross-validated inter-center distances) are an optional output.
% 
%    nBootstraps (DEFAULT = 1): A scalar that determines how many bootstraps (or
%       permutations, or both) to compute via recursive calls. A value of 0 or 1
%       means no bootstrapping or permuting is done.
% 
%    'permute'/'stratified'/'bootstrap'/'calcStats' (DEFAULT = 'bootstrap'):
%       string inputs that determine what type(s) of bootstrapping or permutation
%       to do. 'bootstrap'/'stratified' do the same thing: they resample the
%       points with replacement and compute new Hyperspheres on those resamplings.
%       If multiple hyperspheres are given by the Categories object, the
%       bootstrap is 'stratified' because it disallows points belonging to one
%       hypersphere to be resampled in another. 'permute' resamples the category
%       labels of each point, whichout replacement. 'calcStats' does both: it
%       computes nBootstraps bootstraps and nBootstraps permutations.
% 
%       With any of these options, the final output of estimateHypersphere is an
%       array of Hypersphere objects, each with their categories objects which
%       specifies whether and how the points were bootstrapped or permuted. The
%       first Hypersphere will always be the estimate based on the full (non-
%       resampled) set of points. These can all be consolidated into a single
%       SetOfHyps object (hypset) using [Hyersphere/SetOfHyps].meanAndMerge. The
%       full (non-resampled) estimate is made the primary set of center and radius
%       parameters, and the bootstrapped/permuted Hyperspheres are stored in
%       hypset.ci.bootstraps and hypset.ci.permutations.
% 
%       NOTE: if points have been bootstrapped, the 2-fold cross-validation of the
%       centers will not allow the same point to occur in both folds.
% 
%    'meandist'/'distance'/'mcmc'/'jointml'/'gaussian'/'uniformball'/'uniformcube'
%       (DEFAULT = 'meandist'): The estimator used to estimate the Hypersphere(s).
%       Options include:
%          'meandist': An empirically-derived estimator that automatically selects
%             either the Gaussian- or Uniform-based Minimum Variance Unbiased
%             Estimators (MVUEs) based on the skewness of the radius estimates.
%          'distance': An empirically-derived estimator based on the distance
%             matrix.
%          'mcmc': An MCMC-based estimator that estimates the joint posterior of
%             the centers and radii under the assumption of a uniform ball
%             distribution.
%          'jointml': A maximum likelihood estimator that optimizes the centers and
%             radii jointly under the assumption of a uniform ball distribution.
% [NOT FULLY IMPLEMENTED/TESTED YET:]
%          'gaussian': The Gaussian MVUE.
%          'uniformball': The uniform ball MVUE.
%          'uniformcube': The uniform cube MVUE.
% 
%    'nocvdists' (DEFAULT = false): Does not compute cross-validated centers and
%       distances. Very important when using the 'joint' option.
% 
%    'independent'/'joint' (DEFAULT = 'joint'): These only apply when the points
%       input is a 3D array. With the 'joint' option, the third dimension is
%       treated as an index of, for example, frames in the time dimension, where
%       the same points are changing within a common space. Hyperspheres are
%       therefore estimated at once in the whole 3D array, while ignoring distances
%       of Hyperspheres across frames.
% 
%       The 'independent' option estimates each 'frame' along the 3rd dimension
%       independently. This is useful if points were pre-bootstrapped/permuted, but
%       in that case, it is up to the user to keep track of this fact.
% 
%    'normalize'/'raw' (DEFAULT = 'raw'): 'normalize' standardizes each dimension
%       in points to have zero mean and unit variance. 'raw' just ensures this
%       does not happen.
% 
%    cv (DEFAULT = cvindex(n,2)): A cvindex object used for center/distance cross-
%       validation. This is important for bootstrapping: estimateHypersphere will 
%       use this to ensure that resampled duplicate points don't show up in both
%       folds.
%
% SEE ALSO HYPERSPHERE, HYPERSPHERE.ESTIMATE, HYPERSPEHERE.CONCAT,
%    HYPERSPHERE.MOVIE, SETOFHYPS, SETOFHYPS.MEANANDMERGE,, SETOFHYPS.H2S,
%    CATEGORIES, CATEGORIES.PERMUTE, CATEGORIES.INTERNALREPMAT, CATEGORIES.SELECT,
%    CATEGORIES.SLICE

%% preparation
[n,d,f] = size(points);
nBootstraps= 1;
SAMPLING   = 'calcStats'; % stratified bootstrap (default)
NORMALIZE  = false;       % pre-normalize points
CVDISTS    = true;        % distances are cross-validated by default
cvdists    = false;       % placeholder if don't have cross-validated distance values
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
         case 'normalize';   NORMALIZE  = true;  varargin{v} = [];
         case 'raw';         NORMALIZE  = false; varargin{v} = [];
         case 'nocvdists';     CVDISTS  = false;
         case 'cvdists';       CVDISTS  = true;
         case 'independent'; INDEPENDENT= true;  varargin{v} = [];
         case 'joint';       INDEPENDENT= false; varargin{v} = [];
         case {'bootstrap','stratified','permute','jackknife','calcstats'};
               SAMPLING = varargin{v};
         case {'meandist','distance','mcmc','maxradius','jointml',...
               'gaussian','uniformball','uniformcube'}
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
   SAMPLING   = 'permute';  % bootstrap must be w/o replacement; stratified doesn't make sense

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
         if CVDISTS
            parfor i = 1:f
               [hyp(i),loc_cv{i}] = estimateHypersphere(points(:,:,i),categories,'raw',varargin{:});
            end
            loc_cv = reshape(cat(1,loc_cv{:}),[nCats f])';
         else
            parfor i = 1:f
               hyp(i) = estimateHypersphere(points(:,:,i),categories,'raw',varargin{:});
            end
         end
         return
      end
      points  = reshape(permute(points,[1 3 2]),n*f,d);
      if NORMALIZE
         points = normalize(points,normtype);
      end

      % EVALUATE ALL TOGETHER!
      % Recursive call
      hypEst = estimateHypersphere(points,categories.internalrepmat(f),...
                                   'raw',nBootstraps,'nocvdists',varargin{:});
      % SEPARATE COMBINED HYP
      centers = mat2cell(hypEst.centers,nCats*ones(f,1));
      radii   = num2cell(reshape(hypEst.radii,[nCats f]),1)';
      %loc_cv  = reshape(loc_cv,[nCats f])';
      %if CVDISTS, cvdists = cvCenters2cvSqDists(loc_cv); end

      hyp = Hypersphere(centers,radii,categories,'commonSpace',true);%,'cvdists',cvdists);
      return
   end

   %% permutation or stratified bootstrap test
   if nBootstraps > 1 && ~strcmp(ESTIMATOR,'mcmc')
      catperm = [categories; categories.permute(nBootstraps,SAMPLING)];
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
         cvdists = false;
      else
         [hyp,loc_cv] = cellfun(@(x,cv) estimateHypersphere(x,'raw',nBootstraps,cv,varargin{:}),...
                                     p,num2cell(cvs),'UniformOutput',false);
         % Compute cross-validated distances
         cvdists = cvCenters2cvSqDists(loc_cv);
      end

      switch(ESTIMATOR)
         case 'mcmc'
            % do some fancy combining
            hyp = cat(2,hyp{:}); % n-cell of nSamples to [nSamples x n] array
            % nSamples n-Hyperspheres
            if nBootstraps > 1
               catperm = categories.permute(nBootstraps,true);
               hyp = cellfun(@(x,c) x.concat(c,cvdists),num2cell(hyp,2),num2cell([categories;catperm]));
            else
               hyp = hyp.concat(categories,cvdists);
            end
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
   case {'distance','mcmc'} % do nothing
   case {'jointml','maxradius'}
      loc =  mean(points,1);
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

      % Store samples as bootstraps
      if any(strcmpi(varargin,'bootstrap') | strcmpi(varargin,'permute') | strcmpi(varargin,'calcstats')) 
         hypBoots = Hypersphere(num2cell(loc,2),num2cell(rad));
      end

      % just take the best estimate
      ix     = maxix(llh);
      loc    = loc(ix,:);
      rad    = rad(ix);
      
      if CVDISTS
         loc_cv    = cv.crossvalidate(@inferHyperspherePosterior,...
                                               points,nBootstraps,DBUGPLOT,VERBOSE);
         loc_cv    = squeeze(mat2cell(permute(cell2mat_concat(loc_cv,3),[2 3 1]),...
                                      d,2,ones(numel(llh),1)));
         loc_cv = loc_cv{ix};
      end
   case 'maxradius'
      rad = maxRadiusGivenCenter(loc,points);
      if CVDISTS, loc_cv = cell2mat_concat(cv.crossvalidate(@(p) ...
                                           maxRadiusGivenCenter(loc,p),points));
      end
   case 'jointml'
      tol  = 10^(-7-log2(d));
      opts = optimoptions('fminunc','TolX',tol,'TolFun',tol,'Algorithm','quasi-newton',...
                          'Display','off','GradObj','on');%,'DerivativeCheck','on');
      % search for the center that minimizes the maximum radius, starting at centroid
      objfcn    = @(points) fminunc(@(c) maxRadiusGivenCenter(c,points),loc,opts);
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
if exist('hypBoots','var')
   hyp = [hyp;hypBoots(:)];
end

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

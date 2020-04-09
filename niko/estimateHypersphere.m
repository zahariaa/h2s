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

for v = 2:nargin
   if isa(varargin{v-1},'Categories'), categories = varargin{v-1};
   elseif  isnumeric(varargin{v-1}),   nBootstrapSamples = varargin{v-1};
   elseif  ischar(   varargin{v-1})
      switch(lower(varargin{v-1}))
         case 'permute';     STRATIFIED = false;
         case 'stratified';  STRATIFIED = true;
      end
   end
end
if ~exist('nBootstrapSamples','var'),  nBootstrapSamples = 1; end
%% recurse cell array of points
if exist('categories','var') && numel(categories.labels)>1

%% EVALUATE ALL FRAMES IN SAME SPACE
   if f > 1
      points  = reshape(permute(points,[1 3 2]),n*f,d);
      nCats   = numel(categories.labels);

      % EVALUATE ALL TOGETHER!
      hypEst  = estimateHypersphere(points,categories.internalrepmat(f));
      % SEPARATE COMBINED HYP
      centers = mat2cell(hypEst.centers,nCats*ones(f,1));
      radii   = num2cell(reshape(hypEst.radii,[nCats f]),1)';
      hyp     = Hypersphere(centers,radii,categories);
      return
   end
%% permutation or stratified bootstrap test
   if nBootstrapSamples > 1
      catperm = [categories; categories.permute(nBootstrapSamples,STRATIFIED)];
      for i = 1:nBootstrapSamples+1
         hyp(i) = estimateHypersphere(points,catperm(i));
      end
   else
      for i = 1:numel(categories.labels)
         p{i} = points(~~categories.select(i).vectors,:);
      end
      hyp = cell2mat_concat(cellfun(@estimateHypersphere,p,'UniformOutput',false));
      hyp = hyp.concat(categories);
   end
   return
end

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


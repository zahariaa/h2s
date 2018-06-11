function varargout = estimateHypersphere(points,varargin)
% hyp = estimateHypersphere(points,nBootstrapSamples) % output is a Hypersphere obj
% hyps = estimateHypersphere(points,categories) % output is an array of Hypersphere objs
% [loc,locCI,rad,radCI] = estimateHypersphere(points,nBootstrapSamples)
% uses expected distance between pairs of points to estimate the radius of
% the hypersphere.
% uses the mean (or median) to estimate the location of the center.
% uses bootstrap resampling of the points to estimate 95% confidence
% intervals for both.

for v = 2:nargin
   if isa(varargin{v-1},'Categories'), categories = varargin{v-1};
   elseif  isnumeric(varargin{v-1}),   nBootstrapSamples = varargin{v-1};
   end
end
if ~exist('nBootstrapSamples','var'),  nBootstrapSamples = 1; end
%% recurse cell array of points
if exist('categories','var')
   for i = 1:numel(categories.labels)
      p{i} = points(~~categories.vectors(:,i),:);
   end
   hyp = cellfun(@(p) estimateHypersphere(p,nBootstrapSamples), ...
                      p,'UniformOutput',false);
   varargout{1} = [hyp{:}];
   return
end

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
loc = mean(points,1); % optimises L2 cost, corresponds to L1 force equilibrium

%% estimate radius via skew-based MVUEs
points = points - repmat(loc,[n 1]);
radii = sqrt(sum(points.^2,2));
dvarrad = log2(std(radii/median(radii)))+log2(d)/1.5;
if dvarrad > -0.5 % Assume Gaussian
   rad = stdEst(points)*sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
else             % Assume Uniform
   rad = stdPer(d)*std(radii) + median(radii);
end

if nargout==1, varargout = {Hypersphere(loc,rad)};
else           varargout = {loc,locCI,rad,radCI};
               varargout = varargout(1:nargout);
end

return
end


%% auxiliary functions
function s = stdEst(X)
% Estimate standard deviation by calculating variance along shorter matrix dimension
% (for numerical accuracy)
[n,d] = size(X);
s = sqrt(sum(var(X,[],maxix([d n])))/(max([n d])-1));
return
end

function v = stdPer(d)
expectedStds = [1.2733    1.0115    0.8796    0.8107    0.8384    0.8638    0.9579    1.0403    1.1938  1.4268    1.8384    2.4485];
v = interp1(2.^(1:12),expectedStds,d);
return
end



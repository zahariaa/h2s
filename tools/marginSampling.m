function hyps = marginSampling(points,varargin)
% marginSampling: Does either jackknife- or bootstrap-based sampling to
%    assess margin significance.
% 
% NOTE: ONLY JACKKNIFE-BASED TEST IS IMPLEMENTED RIGHT NOW.
% 
% Required inputs:
%    points: [n x d] numeric array
%    categories: Categories object specifying point-category labels
% 
% Optional inputs:
%    samplingMethod: bootstrap, jackknife (DEFAULT)
% NOT YET IMPLEMENTED{
%    bootstrap + maxradius estimator
%    alternative radius percentile choices (either 100 percentile (max), or other)
% }
% 
% Output:
%    hyps: Hyperspheres of resampled data
% 
% SEE ALSO SETOFHYPS.SIGNIFICANCE, HYPERSPHERE.MEANANDMERGE, CATEGORIES.LEAVEONEOUT

%% parse inputs
for v = 1:numel(varargin)
   if           isa(varargin{v},'Categories'),             cats = varargin{v};
   elseif isnumeric(varargin{v}) && numel(varargin{v})==1, N    = varargin{v};
   elseif    ischar(varargin{v})
      switch  lower(varargin{v})
         case 'bootstrap',                       samplingMethod = varargin{v};
         otherwise                               samplingMethod = 'jackknife';
      end
   end
end
[n,d,f] = size(points);

%% recurse on each category
if ~exist('N'  ,'var') || isempty(N), N = 1000; end
if exist('cats','var') %&& ~isempty(cats)
   nCats = size(cats.vectors,2);
   hyps = cellfun(@(p) marginSampling(p,N,samplingMethod),cats.slice(points),...
                  'UniformOutput',false);
   hyps = cell2mat_concat(hyps,2);
   switch lower(samplingMethod)
      case 'bootstrap' % NOT YET IMPLEMENTED
         % hyps = cell2mat_concat(cellfun(@(h) h.concat(cats),num2cell(hyps,2),...
         %          'UniformOutput',false))';
%% OR ACTUALLY JUST CALL ESTIMATEHYPERSPHERE STRAIGHT UP HERE
         hyps = estimateHypersphere(points,cats,'maxradius','bootstrap',N);
      case 'jackknife'
         % actually have to combine hyps for leave one out over *ALL* points
         [~,thiscat] = find(cats.vectors);
         notthiscat = cell2mat_concat(arrayfun(@(c) setdiff(1:nCats,c),thiscat,'UniformOutput',false));
         hyps = num2cell([vectify(hyps(2:end,:)) reshape(hyps(1,notthiscat),size(notthiscat))],2);
         hyps = cell2mat_concat(cellfun(@(h,c) concat(h,c),...
                          hyps,num2cell(cats.leaveoneout),'UniformOutput',false))';
   end
   return
end

% %% TODO: different percentiles

%% continue on to centroid/radius estimation
loc = mean(points,1); % full data centroid

switch lower(samplingMethod)
   case 'bootstrap' % NOT YET IMPLEMENTED
      % bootstrap this:         maxRadiusGivenCenter(loc,points)
      % don't forget to output bootstrapped cats
      keyboard
   case 'jackknife'
      % leave one out (but first one is full data estimate)
      locs = repmat(loc,[n+1 1]) - [zeros(1,d); points/n];
      locs(2:end,:) = (n-1)*locs(2:end,:)/n;
end

rads = cellfun(@(c) maxRadiusGivenCenter(c,points),num2cell(locs,2));

hyps = Hypersphere(num2cell(locs,2),num2cell(rads));

return
end


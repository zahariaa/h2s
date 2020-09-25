classdef SetOfHyps < Hypersphere
   % Inherited from Hypersphere:
   % properties (SetObservable, GetObservable, AbortSet)
   %    centers    % [n x d] matrix: n centers of d-dimensional hyperspheres
   %    radii      % [1 x n] vector: radii of each of n hyperspheres
   %    categories(1,1) % a Categories object (see help Categories)
   % end
   properties (SetAccess = protected)
      volume     % volume of hypersphere(s) (auto-computed)
      dists      % distances between centers of hyperspheres (auto-computed)
      margins    % margins between hyperspheres along center-connection lines (auto-computed)
      ci         = [];  % confidence intervals, with raw bootstrap data
      error      = NaN; % [n-choose-2 x 1] vector of errors for each summary stat
      sig        = [];  % significance tests results
      sigp       = [];  % significance tests p-values
      sigdiff    = [];  % significant differences tests results
      sigdiffp   = [];  % significant differences tests p-values
      msflips    = [];  % margin sign flips
   end
   properties (Dependent, SetAccess = protected)
      overlap    % overlaps (negative margins) between hyperspheres along center-connection lines
   end
   properties (GetAccess = 'private', SetAccess = 'private')
      boxPos     = [0.65 0.65];
      boxSize    = 0.3;
      nSigLevs   = 3; % 3 means [0.95 0.99 0.999] levels
   end

   methods
      function obj = SetOfHyps(centers,radii,varargin)
      % Constructor for SetOfHyps object for hypersphere2sphere, Hypersphere
      % e.g.:
      % hypset = SetOfHyps(centers,radii,categories)                     % (1)
      % hypset = SetOfHyps('estimate',points,categories,extraargs)       % (2)
      % hypset = SetOfHyps(hyp,ci)                                       % (3)
      % hypset = SetOfHyps(hypstruct)                                    % (4)
      % hypset = SetOfHyps({hyps},radii,extraargs)                       % (5)
      % hypset = SetOfHyps(centers,radii,categories,'cvdists',cvdists)   % (6)
      % hypset = SetOfHyps(hypset,hypsetTarget)                          % (7)
      % 
      % Constructor input options:
      %    (1-2) Same as Hypersphere input options (1-2), but yields a SetOfHyps.
      %    (3) Similar to Hypersphere input option (3), except this SetOfHyps
      %       constructor: (a) converts a Hypersphere object to a SetOfHyps
      %       object (instead of vice-versa), and (b) can take in an optional
      %       ci argument, which is a structure containing the confidence
      %       intervals and bootstrapped Hypersphere/SetOfHyps objects (computed
      %       using SetOfHyps.significance), and saves it in the SetOfHyps object.
      %    (4-6) Same as Hypersphere input options (4-6), but yields a SetOfHyps.
      %    (7) Regardless of the preceding inputs, if a SetOfHyps object is
      %       included as an additional (i.e., not the first) input argument, it
      %       is treated as the "target" SetOfHyps from which errors of the
      %       existing/first input/newly constructed SetOfHyps are computed.
      %       This is equivalent to calling hypset.stressUpdate(hypsetTarget).
      % 
      % N.B.: unlike Hypersphere, volume, dists, margins, and overlap are all
      %    automatically computed once and stored as properties. These are
      %    computed upon SetOfHyps object construction and re-computed if a
      %    change was made to either centers or radii. The methods for
      %    computing these quantities are inherited from Hypersphere.
      % 
      % Inherited Properties (from Hypersphere):
      %    centers    - [n x d] matrix: n centers of d-dimensional hyperspheres
      %    radii      - [1 x n] vector: radii of each of n hyperspheres
      %    categories - a Categories object (see help Categories)
      %    distsCV - (protected) false or vector of squared cross-validated distances
      % Properties:
      %    volume     - volume of hypersphere(s) (auto-computed)
      %    dists      - distances between centers of hyperspheres (auto-computed)
      %    margins    - margins between hyperspheres along center-connection
      %                    lines (auto-computed)
      %    ci         - confidence intervals, with raw bootstrap data
      %    error      - [n-choose-2 x 1] vector of errors for each summary stat
      %    sig        - significance tests results
      %    sigdiff    - significant differences tests results
      %    msflips    - margin sign flips
      %    overlap    - overlaps (negative margins) between hyperspheres along
      %                    center-connection lines (auto-computed)
      % 
      % Inherited Methods (from Hypersphere):
      %    SetOfHyps.resetRandStream
      %    SetOfHyps.getRandStream
      %    SetOfHyps.isempty
      %    SetOfHyps.select
      %    SetOfHyps.snip
      %    SetOfHyps.concat
      %    SetOfHyps.unconcat
      %    SetOfHyps.meanAndMerge
      %    SetOfHyps.sample
      %    SetOfHyps.plotSamples
      %    SetOfHyps.show
      %    SetOfHyps.camera
      %    SetOfHyps.movie
      % Methods:
      %    SetOfHyps.significance
      %    SetOfHyps.h2s
      %    SetOfHyps.plotOverlapErrors
      %    SetOfHyps.plotDynamics
      %    SetOfHyps.showSig
      %    SetOfHyps.showValues
      %    SetOfHyps.showSigLegend
      %    SetOfHyps.toJSON
      %    SetOfHyps.stressUpdate
      %    stress (Static method: h2s objective function and gradients)
      % 
      % 2018-06-07 AZ Created
      % 
      % See also HYPERSPHERE, CATEGORIES, ESTIMATEHYPERSPHERE, SHOWMODEL
         if nargin==0; return; end

         if isa(centers,'SetOfHyps'), obj = centers;    % input option  (7)
         elseif isa(centers,'Hypersphere')              % input option  (3)
            if numel(centers) == 1
               for fn = fieldnames(centers)'
                  obj.(fn{1}) = centers.(fn{1});
               end
               if exist('radii','var') && isstruct(radii)
                  obj.ci = radii; % assumes merged Hypersphere object
               end
            else
               obj = SetOfHyps(centers.meanAndMerge,varargin{:});
            end
         else                                           % input options (1-2,4-6)
            % Use Hypersphere constructor to handle other possible inputs
            if ~exist('radii','var'), radii = []; end
            obj = Hypersphere(centers,radii,varargin{:});

            % Deal with multiple Hyperspheres, convert to SetOfHyps
            MCMC = any(strcmpi(varargin,'mcmc'));
            if obj(end).categories.ispermuted || MCMC
               obj = obj.meanAndMerge(MCMC);
            else
               % assume multiple objects are not bootstraps, but different
               % frames etc, leave them as multiple objects
               obj = arrayfun(@(x) SetOfHyps(x,[],varargin{:}),obj);
               return
            end
         end

         for v = 1:numel(varargin)                      % input option  (7)
            if     isa(varargin{v},'SetOfHyps'),  hypsTarget     = varargin{v};
            elseif isa(varargin{v},'Categories'), obj.categories = varargin{v};
            end
         end

         if ~isa(obj.categories,'Categories')
            obj.categories = Categories(numel(obj.radii));
         end         

         % Update dependent properties, but save them to reduce on-line compute
         obj.volume = obj.volume();
         obj = setPostCenters(obj,[],[],true);

         % Set up listeners to update dependent properties only on centers/radii set
         addlistener(obj,'centers','PostSet',@obj.setPostCenters);
         addlistener(obj,'radii'  ,'PostSet',@obj.setPostRadii);

         % Update error if another SetOfHyps given
         if exist('hypsTarget','var'), obj = obj.stressUpdate(hypsTarget); end
      end

      %% SET FUNCTIONS TO COMPUTE AND STORE VOLUME, DISTS, MARGINS, OVERLAPS
      function self = set.volume(self,~)   % Compute volume of a single hypersphere
         self.volume = Hypersphere.calcVolume(self.centers,self.radii);
      end
      function self = set.dists(self,~)    % Compute distance between two hyperspheres
         self.dists = Hypersphere.calcDists(self.centers);
      end
      function self = set.margins(self,~)
         self.margins = Hypersphere.calcMargins(self.radii,self.dists);
      end
      function ov = get.overlap(self,~)    % Compute overlap distance of two hyperspheres
         ov = -self.margins;
      end


      %% HIGHER ORDER FUNCTIONS
      function [self,loc_cv] = significance(self,points,N,varargin)
      % SetOfHyps.significance: tests 'statistics of interest' (distances,
      %    margins/overlaps, radii), and their pairwise differences, for
      %    significance. It uses bootstrapped and permuted SetOfHyps or
      %    Hypersphere objects found in self.ci.bootstraps and
      %    self.ci.permutations or if they don't exist, does the bootstrapping
      %    /permutation on inputted points, and saves the bootstraps/permutations
      %    (and 95% confidence intervals) in self.ci. Alternatively, the ci
      %    struct can be passed as an input and embedded into the SetOfHyps
      %    object. One can choose to run only permutation- or bootstrap-based
      %    tests with the 'permute' or 'bootstrap'/'stratified' input options.
      % 
      %    Significance testing is based on the percentile confidence interval
      %    that contains 0 (for primary significance tests)*. Significant
      %    difference testing involves computing the percentile confidence
      %    interval (two-tailed) of the difference of two of the same type of
      %    statistics of interest from different hyperspheres. This is true for
      %    all tests except (first-order) distance significance, which is a
      %    permutation-based test asking if the distance computed on the
      %    unpermuted data is within the 95% confidnece interval of the
      %    permuted distances.
      %
      %    p-values are stored in self.sigp and self.sigdiffp, and logicals
      %    indicating whether the null hypothesis was rejected (true) in self.sig
      %    and self.sigdiff.
      %
      %    *Radii are always significant by definition.
      % 
      % e.g.:
      % hyps = hyps.significance(points)
      % hyps = hyps.significance(points,N)
      % hyps = hyps.significance(points,N,'permute') % only runs permutation tests
      % hyps = hyps.significance(ci)
      % hyps = hyps.significance           % works only if hyps.ci is populated
      % 
      % SEE ALSO ESTIMATEHYPERSPHERE, HYPERSPHERE.MEANANDMERGE,
      %    SETOFHYPS.NULLHYPOTHESISREJECTED
         loc_cv = [];
         testtype = 'calcStats';
         for v = 1:numel(varargin)
            if ischar(varargin{v})
               switch lower(varargin{v})
                  case{'calcstats','permute','stratified','bootstrap'}
                     testtype = varargin{v};
                     varargin{v} = [];
               end
            end
         end
         if ~exist('N','var') || isempty(N), N=100; end
         if ~exist('points','var') && isstruct(self.ci) && ~isempty(self.ci.bootstraps)
            % then use self.ci.bootstraps
         elseif isstruct(points) && ~isempty(points)
            self.ci = points;
         elseif ~size(self.categories.vectors,1)==size(points,1)
            error('self.categories.vectors needs to index points');
         else
            % Compute confidence intervals on bootstrapped overlaps & radii for significance
            [bootsNperms,loc_cv] = Hypersphere.estimate(points,self.categories,N,...
                                                        testtype,varargin{:});
            bootsNperms = bootsNperms.meanAndMerge;
            self.ci = bootsNperms.ci;
            if islogical(self.distsCV)
               self.distsCV = bootsNperms.distsCV;
            end
         end
         %% set up significance measures
         n   = numel(self.radii);
         nc2 = nchoosek(n,2);
         % compute indices of radii corresponding to overlaps
         ix = nchoosek_ix(n);
         if n < 3, nc2c2 = 0;
         else      nc2c2 = nchoosek(nc2,2);
                   ixc2  = nchoosek_ix(nc2);
         end

         self.sigp = struct('ra',[],'ma',NaN(1,nc2),'ov',NaN(1,nc2),'di',NaN(1,nc2));
         self.sigdiffp = struct('ra',NaN(1,nc2),...
                                'ma',NaN(1,nc2c2),'ov',NaN(1,nc2c2),'di',NaN(1,nc2c2));

         %% Collect relevant permutations/bootstraps
         if numel(self.ci.permutations)>0
            distCV_perm  = cat(1,self.ci.permutations.distsCV);
            %dist_perm    = cat(1,self.ci.permutations.dists);
         end
         if numel(self.ci.bootstraps)>0
            radii_boot   = cat(1,self.ci.bootstraps.radii);
            margin_boot  =       self.ci.bootstraps.margins;
            dist_boot    = cat(1,self.ci.bootstraps.dists);
         end
         %% COMPUTE SIGNIFICANCE (smaller values more significant)
         % helper functions
         smallertail     = @(x) min(cat(3,1-x,x),[],3);
         ciprctile       = @(x,t) sum([-Inf(1,size(x,2));x;Inf(1,size(x,2))]<=t)/(size(x,1)+2);
         ciprctileSmTail = @(x,t) smallertail(ciprctile(x,t));
         % What percentile confidence interval of bootstrapped margins contains 0?

         if numel(self.ci.bootstraps)>0
            % What percentile confidence interval of bootstrapped overlap/margins contains 0?
            % Overlap/margin is significant if 0 does not exceed the lowest 5% range
            for i = 1:nc2
               self.sigp.ma(i) = ciprctile(margin_boot(:,i),0);
            end
            self.sigp.ov = 1-self.sigp.ma;
         end
         if numel(self.ci.permutations)>0
            % At what percentile confidence interval of permuted distances does
            % the unpermuted distance estimate occur?
            for i = 1:nc2
               self.sigp.di(i) = ciprctileSmTail(distCV_perm(:,i),self.distsCV(i));
            end
         end

         %% SECOND-ORDER COMPARISONS
         if numel(self.ci.bootstraps)>0
            for i = 1:nc2
               self.sigdiffp.ra(i) = ciprctileSmTail(diff(  radii_boot(:,  ix(:,i)),[],2),0);
            end
            for i = 1:nc2c2
               self.sigdiffp.ov(i) = ciprctileSmTail(diff(-margin_boot(:,ixc2(:,i)),[],2),0);
               self.sigdiffp.di(i) = ciprctileSmTail(diff(   dist_boot(:,ixc2(:,i)),[],2),0);
            end
            self.sigdiffp.ma = self.sigdiffp.ov; % Margin/overlap sigdiff is same bc smaller tail
         end

         self.sig     = SetOfHyps.nullHypothesisRejected(self.sigp,0.05);
         self.sigdiff = SetOfHyps.nullHypothesisRejected(self.sigdiffp,0.05);

         % %% DEBUG distances
         % fh = newfigure([nc2 1]);
         % for i=1:nc2
         %    axtivate(i)
         %    hist(dist_boot(:,i));
         %    plot((self.dists(i)^2)*[1 1],[0 N/nc2],'r-','LineWidth',2);
         % end
         % matchy(fh.a.h,{'XLim','YLim'})
      end

      function self = stressUpdate(self,hypsTarget)
      % SetOfHyps.stressUpdate: populates error and msflips properties, computed
      %    by evaluating the error between the self SetOfHyps object and the
      %    hypsTarget SetOfHyps input.
      % 
      % e.g.:
      % hypset = hypset.stressUpdate(hypsTarget)
      % 
      % SEE ALSO SETOFHYPS.STRESS
         for i = 1:numel(self)
            [~,~,self(i).error,self(i).msflips] = self(i).stress(self(i),hypsTarget(i));
         end
      end
      
      function model = h2s(self,varargin)
      % SetOfHyps.h2s: optimizes the hypersphere2sphere algorithm, and outputs
      %    a SetOfHyps object that has dimensionality reduced to dimLow (2 or 3).
      %    If self is an array of SetOfHyps objects, self.h2s will be recursively
      %    run on each element.
      % 
      % e.g.
      % hypslo = hypshi.h2s
      % hypslo = hypslo.h2s(hypsTarget)
      % hypsjointlo = severalhyps.h2s(hypsTargets,'joint') % jointly optimizes all 
      %    SetOfHyps objects in severalhyps (useful for movies)
      % hypslo = hypshi.h2s([dimLow ninits]).  % (optional numeric inputs)
      % hypslo = hypshi.h2s({alpha})           % alpha values must be in a cell
      % hypslo = hypshi.h2s('mdsinit','debugplot','verbose')  % (string flags)
      %    ... or any combination of the above, with the restriction that if ninits
      %    is specified, dimLow must also be specified first in the array.
      % hypslo = hypshi.h2s('mdsinit','debug') % equivalent to previous, if all
      %                                          inputs in previous are provided
      % hypslo = hypshi.h2s('randinit')        % (as opposed to 'mdsinit')
      % 
      % Optional inputs:
      %    hypsTarget (DEFAULT = self): a SetOfHyps object, which contains the
      %       high-dimensional summary statistics of interest against which the
      %       h2s stress function will measure errors and be optimized to match.
      %    dimLow (DEFAULT = 3, or the number of hyperspheres, or the
      %       dimensionality of the high-dimensional space, whichever is
      %       smallest): a scalar numeric, either 2 or 3, that determines the
      %       dimensionality of the visualization space.
      %       N.B. providing a negative value tells h2s to use the default value,
      %       so an ninits value can be provided.
      %    ninits (DEFAULT = 0): the number of additional (random)
      %       initializations (plus the 1 MDS/random initialization) for stress
      %       function optimization, keeping the best result. Must be provided
      %       as a second element in a numeric vector.
      %       e.g.:
      %          self.h2s([3 5])  sets dimLow to 3 and sets ninits to 5
      %          self.h2s(3)      sets dimLow to 3
      %          self.h2s([-1 5]) uses dimLow default and sets ninits to 5
      %    {alpha} (DEFAULT = {1 1 1}): a numeric 3-element cell (or 3-vector
      %       in a cell) specifying the regularization hyperparameters on the
      %       3 components of the objective function: the elements weight the
      %       distances, margins, and radii, in that order, in the error
      %       (stress) function.
      %       e.g. the following disables margins in error:
      %          self.h2s({[1 0 1]}) or self.h2s({1 0 1})
      %          self.h2s({1 1 0}) and self.h2s({1 0 0}) keep radii fixed, and
      %                              are both equivalent to MDS on the centers
      %    String input options:
      %       'joint'/'separate' (DEFAULT = 'joint'): if an array of multiple
      %          SetOfHyps objects are provided, h2s jointly optimizes them. If
      %          'separate' option provided, recursive calls to h2s optimize them
      %          separately.
      %       'mdsinit'/'randinit' (DEFAULT = 'mdsinit'): use MDS for low centers
      %          initialization (as opposed to random initialization(s)).
      %       'debugplot' (DEFAULT = false): display the h2s optimization
      %          debugging plot, which updates on each iteration. Must be second
      %          element in logical vector.
      %       'verbose' (DEFAULT = false): prints fmincon optimization messages
      %          to the command prompt on every iteration. Must be third element
      %          in logical vector.
      %       'debug': equivalent to 'debugplot' and 'verbose'
      %       'quiet' (DEFAULT): debug and verbose are off
      %          e.g.:
      %             self.h2s('debug') uses MDS for low centers
      %                initialization, displays the fit debugging plot, and
      %                prints optimization messages to the command prompt.
      % 
      % SEE ALSO SETOFHYPS.STRESS, SETOFHYPS.STRESSUPDATE, SETOFHYPS.STRESSPLOTFCN
         lo = self;
         for v = 2:nargin
            if isa(varargin{v-1},'SetOfHyps')
               hi = varargin{v-1};
            elseif isnumeric(varargin{v-1})
               [dimLow,ninits] = dealvec(varargin{v-1});
               if numel(varargin{v-1})==3
                  varargin{v-1}(3) = 1; % change to 1 for recursive call below
               end
            elseif iscell(varargin{v-1})
               alpha = [varargin{v-1}{:}];
            else
               switch lower(varargin{v-1})
                  case 'joint';     JOINT    = true;
                  case 'separate';  JOINT    = false;
                  case 'mdsinit';   MDS_INIT = true;
                  case 'randinit';  MDS_INIT = false;
                  case 'debugplot'; DBUGPLOT = true;
                  case 'verbose';   VERBOSE  = true;
                  case 'debug';     DBUGPLOT = true;   VERBOSE  = true;
                  case 'quiet';     DBUGPLOT = false;  VERBOSE  = false;
               end
            end
         end
         if ~exist('hi'      ,'var'),                 hi = self;          end
         n  = numel(hi);
         [nr,d] = size(cat(1,hi.centers));
         if ~exist('dimLow'  ,'var') || dimLow<0, dimLow = min([3 nr d]); end
         if ~exist('ninits'  ,'var'),             ninits = 0;             end
         if ~exist('MDS_INIT','var'),           MDS_INIT = true;          end
         if ~exist('DBUGPLOT','var'),           DBUGPLOT = false;         end
         if ~exist('VERBOSE' ,'var'),            VERBOSE = false;         end
         if ~exist('JOINT'   ,'var'),              JOINT = true;          end
         if ~exist('alpha'   ,'var'),              alpha = [1 1 1];       end

         % Recurse for multiple SetOfHyps objects
         if ~JOINT && n > 1
            for i = 1:n
               stationarycounter(i,n);
               model(i) = lo(i).h2s(varargin);
            end
            return
         end

         % Setup optimization and run
         if MDS_INIT
            multimds = arrayfun(@(h) mdscale(h.dists,dimLow,'Criterion','metricstress'),...
                                     hi,'UniformOutput',false);
            x0  = [cat(2,hi.radii)' cell2mat_concat(multimds)];
         else
            allcentervals = cat(1,hi.centers);
            cminmax = [min(allcentervals(:)) max(allcentervals(:))];
            %% initialization of centers: random noise added to projection onto first 2 dimensions
            x0  = [cat(2,hi.radii)' ...
                     -0.1*diff(cminmax)*rand(nr,dimLow) + allcentervals(:,1:dimLow)];
         end
         
         if (all(alpha(1:2)==1) && alpha(3)==0)
            fit = x0;
         else
            opts = optimoptions(@fmincon,'TolFun',1e-4,'TolX',1e-4,'Display','off',...
                        'SpecifyObjectiveGradient',true);%,'DerivativeCheck','on');
            if DBUGPLOT, opts.OutputFcn = @self.stressPlotFcn; end
            if VERBOSE,  opts.Display   = 'iter';              end
            lb = [zeros(nr,1) -Inf(nr,dimLow)];
            ub = Inf(nr,1+dimLow);
            fit = fmincon(@(x) SetOfHyps.stress(x,hi,alpha),x0,...
                          [],[],[],[],lb,ub,[],opts);
         end
         % Output reduced SetOfHyps model
         nc = size(fit,1)/n;
         model = SetOfHyps(mat2cell(fit(:,2:end),repmat(nc,[n 1])),...
                           mat2cell(fit(:,1    ),repmat(nc,[n 1])),self.categories);
         model.stressUpdate(hi);

         % Recurse for multiple random initializations
         if ninits > 0
            % Make sure other initializations are random
            if MDS_INIT
               for v = 2:nargin
                  if islogical(varargin{v-1})
                     MDS_INIT = false;
                     varargin{v-1}(2) = MDS_INIT;
                  end
               end
            end
            % force randinit if mdsinit is specified
            tochange = strcmpi(varargin,'mdsinit');
            if any(tochange)
               varargin{tochange} = 'randinit';
            end
            % Loop through recursive function calls
            for iinit = 1:ninits
               tmp = self.h2s(varargin);
               if max(cat(2,tmp.error)) < max(cat(2,model.error))
                  model = tmp;   % keep lower error solutions
               end
            end
         end
      end


      function errlines = plotOverlapErrors(self)
      % SetOfHyps.plotOverlapErrors: plots a line where a true overlap is
      %    visualized as a margin, or a true margin is visualized as an
      %    overlap, but *ONLY IF THE TRUE OVERLAP/MARGIN is SIGNIFICANT*. Uses
      %    self.sig and self.msflips properties, populated by self.significance.
      %    If no information about significance is available, then all overlap/
      %    margins are considered significant and any flip is plotted.
      % 
      % e.g.:
      % hyps.plotOverlapErrors
      % 
      % SEE ALSO HYPERSPHERE.SHOW, SETOFHYPS.SIGNIFICANCE
         if isempty(self.sig), sigov = true(1,numel(self.margins));
         else                  sigov = sum([self.sig.ov;self.sig.ma]);
         end
         if isempty(self.msflips), self.msflips = zeros(size(sigov));
         end
         sigov = sigov .* self.msflips;

         if isempty(sigov) || ~any(sigov)
            errlines = [];
            return
         end
         n = numel(self.radii);
         ix = nchoosek_ix(n);
         for i = find(sigov)
            err      = self.error(n+i); % assumes order of errors is radii, margins, distances
            centers  = self.centers(ix(:,i),:);
            dxy      = diff(centers);
            dxy      = dxy/norm(dxy);
            midpoint = mean(centers);
            % Compute overlap midpoint
            % need to assign +/-1 to larger/smaller center coordinate?
            ov_edges = centers + [1;-1].*self.radii(ix(:,i))'*dxy;
            ov_midpoint = mean(ov_edges);
            % % debug
            % keyboard
            % plot(centers(:,1),centers(:,2),'ko')
            % plot(midpoint(:,1),midpoint(:,2),'ro')
            % plot(ov_edges(:,1),ov_edges(:,2),'go')
            % plot(ov_midpoint(:,1),ov_midpoint(:,2),'rx')

            if size(self.centers,2)==2
               errlines(i) = plot( ov_midpoint(1) + [-1 1]*dxy(1)*err/2,...
                                   ov_midpoint(2) + [-1 1]*dxy(2)*err/2,'k-','LineWidth',2);
            else
               errlines(i) = plot3(ov_midpoint(1) + [-1 1]*dxy(1)*err/2,...
                                   ov_midpoint(2) + [-1 1]*dxy(2)*err/2,...
                                   ov_midpoint(3) + [-1 1]*dxy(3)*err/2,'k-','LineWidth',2);
            end
         end
      end

      function plotDynamics(self,prop,varargin)
      % SetOfHyps.plotDynamics: Plots the progression of specified SetOfHyps
      %    object properties ('prop' string input argument) through the array
      %    of SetOfHyps objects this is called on. Can optionally add a legend,
      %    plot multiple properties, and specify the x-axis values each object
      %    corresponds to (e.g., the times of the frames each object represents).
      % 
      % e.g.:
      % bunchohyps.plotDynamics('dists')     % plots all distance comparisons
      % bunchohyps.plotDynamics({'radii','margins','dists'})
      % bunchohyps.plotDynamics('all')       % equivalent to line above
      % bunchohyps.plotDynamics('all',times) % specifies x-axis values
      % 
      % Required input argument:
      %    prop (string or cell of strings): The summary stat(s) to plot. Can be:
      %       'radii'
      %       'margins'
      %       'dists'
      %       'volume'
      %       'error'
      %       'all': equivalent to {'radii','margins','dists'}
      %       Or a cell of any combination of the above
      % 
      % Optional inputs:
      %    ax (DEFAULT = gca): Axes handle(s) where the plots will be placed.
      %    times (DEFAULT = 1:nframes): x-axis values for each frame, e.g., the
      %       time of that frame.
      %    LEGEND (DEFAULT = false): If true, places a categories.legend at the
      %       top left.
      % 
      % SEE ALSO CATEGORIES.LEGEND
         for v = 3:nargin
            if      ishandle(varargin{v-2}), ax     = varargin{v-2};
            elseif islogical(varargin{v-2}), LEGEND = varargin{v-2};
            elseif isnumeric(varargin{v-2}), times  = varargin{v-2};
            end
         end
         if ~exist('LEGEND','var') || isempty(LEGEND), LEGEND = false; end
         if ~exist('ax'    ,'var') || isempty(ax),     ax     = gca;  end
         if ~exist('times' ,'var') || isempty(times),  times  = 1:numel(self); end
         if ischar(prop) && strcmpi(prop,'all')
            prop = {'radii','margins','dists'};
         end
         if iscell(prop)
            for i = 1:numel(prop)
               self.plotDynamics(prop{i},times,ax(i));
            end
            matchy(ax,'XLim')
            if LEGEND, self(1).categories.legend; end
            return
         end

         n  = numel(self(1).radii);
         nk = nchoosek(n,2);
         ix = nchoosek_ix(n)';

         axtivate(ax);
         set(ax,'ColorOrder',self(1).categories.colors(ix(:),:))
         if strcmpi(prop,'radii') || strcmpi(prop,'volume')
            plot(times,reshape([self.(prop)],n ,[])','LineWidth',1)
         else
            plot(times,reshape([self.(prop)],nk,[])','LineWidth',1)
            plot(times,reshape([self.(prop)],nk,[])','--','LineWidth',1)
         end
         switch prop
            case 'dists',   prop = 'distances';
            case 'margins' % do nothing
            case 'radii',   prop = 'sizes';
            case 'volume',  prop = 'volumes';
            case 'error',   prop = 'errors';
         end
         ylabel(prop)
         xlim([min(times) max(times)]);

         if LEGEND, self(1).categories.legend; end
      end

      function showSig(self,varargin)
      % SetOfHyps.showSig: Visualizes SIGNIFICANT summary statistics (radii,
      %    distances, and overlaps/margins), and/or the significances of their
      %    DIFFERENCES, in an elaborated Hinton diagram.
      % 
      %    As in SetOfHyps.showValues, the diagonal corresponds to radii, the
      %    upper triangle the overlaps/margins, and the lower triangle the
      %    distances. If differences are being plotted, these respective pairwise
      %    differences. Circle means positive value, square means negative value
      %    (relevant for overlaps/margins).
      % 
      %    N.B. As opposed to SetOfHyps.showValues, each shape's grayscale value
      %    corresponds to the level of significance. 
      % 
      % e.g.:
      % hyps.showSig
      % hyps.showSig('legend')           % shows significances with legend
      % hyps.showSig(true)               % shows significances with legend
      % hyps.showSig('sigdiff',axHandle) % significant differences in axHandle
      % hyps.showSig(twoAxHandles)       % sig & sigdiff in each axes handles
      % 
      % Optional inputs:
      %    ax (DEFAULT = gca): Axes handle where the visualization(s) should be
      %       plotted. If two axes are given, showSig automatically plots both
      %       significances and significant differences.
      %    'sig' (DEFAULT): plots significances.
      %    'sigdiff'/'diff': plots significant differences.
      %    'legend': plots a legend explaining the grayscale values and their
      %       corresponding significance levels.
      % 
      % SEE ALSO SETOFHYPS.SHOWVALUES, SETOFHYPS.SIGNIFICANCE, STATSMAT,
      %    SETOFHYPS.SHOWSIGLEGEND, SETOFHYPS.SHAPESINABOX
         for v = 1:nargin-1
            if   ishandle(varargin{v}),      ax  = varargin{v};
            elseif ischar(varargin{v})
               switch lower(varargin{v})
                  case 'legend';            LGND = varargin{v};
                  case 'sig';               sig  = self.sigp;
                                            DIFF = lower(varargin{v});
                  case {'sigdiff', 'diff'}; sig  = self.sigdiffp;
                                            DIFF = lower(varargin{v});
               end
            end
         end
         if ~exist('ax'  ,'var') || isempty(ax)  , ax   = gca;        end
         if ~exist('sig' ,'var') || isempty(sig) , sig  = self.sigp;
                                                   DIFF = 'sig';      end
         if ~exist('LGND','var') || isempty(LGND), LGND = 'nolegend'; end

         if isempty(sig)
            error(['need to populate obj.sigp/sigdiffp first,'...
                   ' e.g. with:\nobj = obj.significance(points);']);
         end

         if numel(ax) == 2
            % if two axes given, plot sig in first, sigdiff in second
            self.showSig(ax(1),LGND);   % if legend requested, put here
            self.showSig(ax(2),'diff');
            return
         end

         nSigLevs = self.nSigLevs; % 3 means [0.05 0.01 0.001] levels

         % Convert to sigmas significance, maxed to nSig(=3) and floored
         % Use False Discovery Rate correction for significance computation
         sigThresh = SetOfHyps.nullHypothesisRejected(sig,[0.05 10.^(-2:-1:-nSigLevs)]);
         sigThresh = structfun(@(x) max(0,x/nSigLevs-1e-5),sigThresh(1),'UniformOutput',false);

         n = numel(self.radii);
         % Time to get a-plottin'
         self.elaborateHinton(statsmat(sigThresh.ov-sigThresh.ma,sigThresh.di,sigThresh.ra),...
                              'color',ax,DIFF)
         if strcmpi(DIFF,'sig'), title('Significant values')
         else                    title('Significant differences')
         end
         if strcmpi(LGND,'legend'), self.showSigLegend(ax); end
      end

      function showValues(self,varargin)
      % SetOfHyps.showValues: Visualizes significant summary statistic (radii,
      %    distances, and overlaps/margins) VALUES in an elaborated Hinton
      %    diagram. A reference SetOfHyps object can be additionally provided to
      %    compare estimated vs visualized values.
      % 
      %    As in SetOfHyps.showSig, the diagonal corresponds to radii, the upper
      %    triangle the overlaps/margins, and the lower triangle the distances.
      %    Circle means positive value, square means negative value (relevant for
      %    overlaps/margins).
      % 
      %    N.B. As opposed to SetOfHyps.showSig, the radii of each shape
      %    represents the corresponding summary statistic's value. These radii
      %    are normalized to the maximum of the absolute value of ALL summary
      %    statistics provided.
      % 
      %    Furthermore, each shape can be split into left/right sides, their
      %    radii corresponding to the high/low, estimated/visualized values. This
      %    happens only when an additional SetOfHyps object is provided as input.
      % 
      % e.g.:
      % hyps.showValues
      % hyps.showValues(hypsTarget) % splits shapes left/right hypsTarget/hyps
      % hyps.showValues(axHandle)   % puts diagram in axHandle
      % 
      % Optional inputs:
      %    hypsTarget (DEFAULT = self): reference (high-dimensional) SetOfHyps
      %       object containing summary statistics to be rendered on left half of
      %       each shape.
      %    axHandle (DEFAULT = gca): Axes handle where the visualization should
      %       be plotted.
      % 
      % SEE ALSO SETOFHYPS.SHOWSIG, SETOFHYPS.SIGNIFICANCE, STATSMAT,
      %    SETOFHYPS.SHOWSIGLEGEND, SETOFHYPS.SHAPESINABOX
         % Parse inputs
         for v = 1:nargin-1
            if isa(varargin{v},'SetOfHyps'), hi = varargin{v};
            elseif ishandle(varargin{v}),    ax = varargin{v};
            end
         end
         if ~exist('ax','var') || isempty(ax), ax = gca; axis ij off; end
         if ~exist('hi','var') || isempty(hi), hi = self; self = [];  end

         himat = statsmat( -hi.margins, hi.dists, hi.radii);

         if ~isempty(self)
            lomat = statsmat(-self.margins,self.dists,self.radii);
            % Normalize both to same scale
            maxval = max(abs([himat(:);lomat(:)]));

              hi.elaborateHinton(himat/maxval,'left', ax)
            self.elaborateHinton(lomat/maxval,'right',ax)
            title('Values comparison')
         else
              hi.elaborateHinton(himat/max(abs(himat(:))),'size',ax)
            title('Values')
         end
      end

      function showSigLegend(self,ax)
      % SetOfHyps.showSigLegend: Renders legend for a significance/significant
      %    differences elaborated Hinton diagram. Indicates the correspondence
      %    between color and significance level.
      % 
      % e.g.:
      % hyps.showSigLegend
      % hyps.showSigLegend(axHandle)
      % 
      % Optional input:
      %    axHandle (DEFAULT = gca): Axes handle where the legend should be
      %       rendered. Should be the same one where SetOfHyps.showSig was/will
      %       be.
      % 
      % SEE ALSO SETOFHYPS.SHOWSIG
         n        = numel(self.radii);
         boxPos   = self.boxPos;
         sqSz     = self.boxSize/n;
         nSigLevs = self.nSigLevs; % 3 means [0.95 0.99 0.999] levels
         
         if ~exist('ax','var') || isempty(ax), ax = gca; axis ij off; end
         
         % Time to get a-plottin'
         set(0,'CurrentFigure',get(ax,'Parent'))
         set(get(ax,'Parent'),'CurrentAxes',ax,'Units','normalized')
         
         %% Draw Legend
         for i = 1:nSigLevs
            rectangle('Position',[boxPos+[0.01+sqSz*(n+1/3) sqSz*(i-1)/3] sqSz/3 sqSz/3],...
                      'FaceColor',[1 1 1]*(1-i/nSigLevs),'EdgeColor','none')
         end
         % Probability labels
            text('Position',[boxPos+[0.01+sqSz*(n+0.7) sqSz/6]],'String','p<0.05')
         for i = 2:nSigLevs
            text('Position',[boxPos+[0.01+sqSz*(n+0.7) sqSz*(2*i-1)/6]],'String',...
                 ['p<' num2str(10^-i)])
         end
         axis equal ij off
      end

      function varargout = toJSON(self, saveFile)
      % SetOfHyps.toJSON(<saveFile>): exports SetOfHyps object as an array of
      %    Hypersphere objects in JSON format. Optional saveFile input argument
      %    specifies where to save the JSON.
      % e.g.:
      % hyps.toJSON             % prints JSON in command line
      % hyps.toJSON('hypsjson') % saves JSON in file named hypsjson.json
      % 
      % SEE ALSO HYPERSPHERE.UNCONCAT
         for i = 1:numel(self)
            self.categories.vectors = [];
         end
         json = jsonencode(Hypersphere(self).unconcat);
         % Saves JSON to file, if requested
         if exist('saveFile','var') && ~isempty(saveFile) && any(saveFile)
            fileID = fopen([saveFile '.json'],'w');
            fprintf(fileID,'%s',json);
            fclose(fileID);
            fprintf('JSON saved to %s.json\n',saveFile);
         end
         if nargout || ~exist('fileID','var')
            varargout = {json};
         end
      end
   end % methods (public)


   methods(Access = 'private')
      %% UPDATE DISTS, MARGINS IF CENTERS CHANGED
      function self = setPostCenters(self,src,event,FIRSTRUN)
         if ~exist('FIRSTRUN','var') || isempty(FIRSTRUN)
            FIRSTRUN = false;
         end
         self.dists   = dists(self);
         self.margins = margins(self);
         if ~FIRSTRUN
            if ~any(isnan(self.error))
               warning(['SetOfHyps object updated, but obj.errors/ci/sig/sigdiff '...
                        'are out of sync']);
            elseif self.distsCV
               warning(['SetOfHyps object updated, but self.distsCV, '...
                        'and thus self.dists, are out of sync'])
            end
         end
      end
      function self = setPostRadii(self,src,event,FIRSTRUN)
         if ~exist('FIRSTRUN','var') || isempty(FIRSTRUN)
            FIRSTRUN = false;
         end
         self.volume  = volume(self);
         self.margins = margins(self);
         if ~FIRSTRUN && ~any(isnan(self.error))
            warning('SetOfHyps object updated, but obj.errors/ci/sig/sigdiff are out of sync');
         end
      end
      function stop = stressPlotFcn(self,varargin)%x,optimValues,state)
      % stressPlotFcn: private function that plots errors and summary statistics
      %    as a function of optimization iteration. This is only relevant as a
      %    debugging plot for h2s.
      % 
      % e.g.
      % hypslo = hypshi.h2s('debugplot')
      % 
      % SEE ALSO SETOFHYPS.H2S, SETOFHYPS.STRESS, SETOFHYPS.STRESSUPDATE
         % unpack inputs
         x           = varargin{1};
         optimValues = varargin{2};
         state       = varargin{3};
         % preliminaries
         n           = size(x,1);
         cols        = colormap('lines');
         % copy hi radii to lo SetOfHyps
         lo          = SetOfHyps(x,self.radii);
         % Recalculate error from high dimensional SetOfHyps self
         lo          = lo.stressUpdate(self);
         loerr       = lo.error.^2;

         % Plot
         figure(99);
         if strcmpi(state,'init')
            clf; subplot(3,3,1); hold on; lo.show; title('initial');
         end
         subplot(3,3,2); hold on; lo.show; title('fit')

         subplot(3,3,3); hold on % Plot error value
         plot(optimValues.iteration,optimValues.fval,'ko')
         ylabel('Error')
         xlabel('Iteration')

         for i = 1:(n^2-n)/2
            subplot(3,3,9); hold on;title('dists'  )
            plot(optimValues.iteration,lo.dists(i)  ,'o','Color',cols(i,:));
            subplot(3,3,8); hold on;title('margins')
            plot(optimValues.iteration,lo.margins(i),'o','Color',cols(i,:));
            subplot(3,3,7); hold on;title('overlaps')
            plot(optimValues.iteration,lo.overlap(i),'o','Color',cols(i,:));
            if size(loerr,1)>2
            subplot(3,3,6); hold on;title('dists'  )
            plot(optimValues.iteration,loerr(3,i),'o','Color',cols(i,:));
            end
            if size(loerr,1)>1
            subplot(3,3,5); hold on;title('margins')
            plot(optimValues.iteration,loerr(2,i),'o','Color',cols(i,:));
            end
            subplot(3,3,4); hold on;title('overlaps')
            plot(optimValues.iteration,loerr(1,i),'o','Color',cols(i,:));
         end
         ylabel('Relative error')

         stop = false;
      end
   end

   methods(Static)
      function [err2mean,grad,err,msflips] = stress(radii_and_centers,hi,alpha)
      % SetOfHyps.stress: Stress/objective function that is optimized for h2s.
      %    Computes analytic gradients, margin/overlap sign flips between high
      %    and low-dimensional summary stats, and can output individual errors.
      % 
      % e.g.:
      % [err2mean,grad,err,msflips] = stress([radii(:) centers],hi)
      % err2mean = stress(lo,hi)
      % err2mean = stress(lo,hi,alpha)
      % [~,~,hyps.error,hyps.msflips] = hyps.stress(hyps,hypsTarget); %stressUpdate
      % fit = fmincon(@(x) SetOfHyps.stress(x,hi,alpha),x0,[],[],[],[],lb,ub,[],opts);
      % 
      % Requires at least two inputs:
      %    radii_and_centers: Can be either: (1) a SetOfHyps object, or
      %       (2) a numeric matrix with the form [radii(:) centers].
      %    hi: The 'target' high dimensional SetOfHyps object, relative to which
      %       the stress function computes its errors.
      % Optional input:
      %    alpha (DEFAULT = [1 1 1]): a numeric 3-vector that contains the
      %       regularization hyperparameters, such that the first element is the
      %       weight on the distance errors, the second is on the overlap errors,
      %       and the third is on the radius errors.
      % 
      % Outputs:
      %    err2mean: Total of all error terms, with regularization (alpha)
      %       hyperparameters applied. If alpha=[1 1 1], err2mean = mean(err.^2)
      %    grad: Computed analytic gradients wrt radii_and_centers.
      %    err: Vector containing each individual error term, with regularization
      %       (alpha) hyperparameters applied.
      %    msflips: Logical vector representing whether a high-dimensional margin
      %       flipped sign and became an overlap, or vice-versa.
      % 
      % SEE ALSO SETOFHYPS.H2S, SETOFHYPS.STRESSUPDATE, HYPERSPHERE.SHOW
         n = numel(hi);
         if ~exist('alpha','var') || isempty(alpha)
            alpha  = [1 1 1]; % for testing gradients with hyperparameters
         end
         if isa(radii_and_centers,'SetOfHyps'), lo = radii_and_centers;
         else % recurse and combine
            nc = size(radii_and_centers,1)/n;
            lo = SetOfHyps(mat2cell(radii_and_centers(:,2:end),repmat(nc,[n 1])),...
                           mat2cell(radii_and_centers(:,1    ),repmat(nc,[n 1])));
            % separately optimize each SetOfHyps (recursively)...
            [err2mean,grad,err,msflips] = arrayfun(@(l,h) ...
                  SetOfHyps.stress(l,h,alpha),lo(:),hi(:),'UniformOutput',false);
            % ...then combine the outputs
            % (to avoid additional/incorrect pairwise comparisons being included)
            err2mean = mean([err2mean{:}]);
            grad     = cat(1,grad{:});
            err      = [err{:}];
            msflips  = [msflips{:}];
            return
         end

         [nc,dimLow] = size(lo.centers);
         errd        = hi.dists   - lo.dists;
         erro        = hi.margins - lo.margins;
         erra        = hi.radii   - lo.radii;
         err         = abs([alpha(1)*errd alpha(2)*erro alpha(3)*erra]); % not exactly correct
         err2mean    = mean([alpha(1)*errd.^2 alpha(2)*erro.^2 alpha(3)*erra.^2]);

         msflips = ~(sign(hi.margins) == sign(lo.margins));

         % Gradient: partial centers derivative of distances (should be size of centers)
         % Gradient: partial derivatives of error, relative to distances and overlaps
         ix = nchoosek_ix(nc);
         nc2 = size(ix,2);
         
         % center gradients
         dddc = zeros(nc2,nc,dimLow);% 1./lo.dists;
         centersdiff = lo.centers(ix(2,:),:) - lo.centers(ix(1,:),:);
         for i = 1:nc2
            dddc(i,ix(:,i),:) = dddc(i,ix(:,i),:) ...
                                + shiftdim([1;-1]*centersdiff(i,:),-1);
         end
         dddc = dddc./repmat(lo.dists',[1 nc dimLow]);
         dEdd = 2*(alpha(1)*errd + alpha(2)*erro );
         dEdc = reshape(dEdd*reshape(dddc,nc2,[]),[nc dimLow]);

         % radius gradients
         dEdr = -2*alpha(3)*erra';
         for i = 1:nc
            dEdr(i) = dEdr(i) + 2*alpha(2)*sum(erro(~~sum(ix==i)));
         end
         % Gradient: put it all together and normalize properly
         grad = [dEdr dEdc]/(2*nc2+nc);
      end
   end

   methods(Hidden = true)
      function elaborateHinton(self,values,varargin)
      % SetOfHyps.elaborateHinton: Generates elaborated Hinton diagrams as
      %    alternative visualizations of the summary statistic values,
      %    significances, and significant differences. Includes color keys
      %    indicating which diagram entry corresponds to which category.
      % 
      % e.g.:
      % hyps.elaborateHinton(hyps.overlaps,'uppert',ax);
      % hyps.elaborateHinton(hyps.dists,'lowert',ax);
      % hyps.elaborateHinton(hyps.radii,'diagonal',ax);
      % hyps.elaborateHinton(statsmat(overlaps,distances,radii)) % ~equivalent to 3 lines above
      % hi.elaborateHinton(hivals,'left',ax);lo.elaborateHinton(lovals,'right',ax);
      % hyps.elaborateHinton(values,'size')  % values-based shape sizing 
      % hyps.elaborateHinton(values,'color') % values-based shape coloring (grayscale)
      % hyps.elaborateHinton(values,'diff')  % forces second-order comparison diagram
      % 
      % Required input:
      %    values: Must be numeric. Can be a vector or matrix. If it's a matrix,
      %       it's assumed to be arranged in the same way as the diagram itself.
      %       If it's a vector, it should be accompanied with a 'boxpart' string
      %       input.
      % 
      % Optional inputs:
      %    colors (DEFAULT: mat2cell(zeros(N,3),ones(1,N)) => all black): A cell
      %       containing the colors of the shapes. The colors can either be
      %       based on the value for the shape, or the identity the shape
      %       corresponds to (as in the case of the radii in a second-order
      %       comparison diagram).
      %    boxpart ('uppert'/'diagonal'/'lowert'): The part of the matrix the
      %       values are meant to populate (upper triangle / diagonal / lower
      %       triangle). DEFAULT is 'uppert' or none (the whole matrix),
      %       depending on the size of the values provided.
      %    side ('left'/'right'): The side of the shapes the values correspond to
      %       (both in terms of square vs circle, the radius, and the color [this
      %       particular scenario is not a default one]). DEFAULT is none (which 
      %       really means both).
      %    FIRSTORDER ('sig'/sigdiff'/'diff'): Forces the diagram to be either
      %       first-order (values or significance) or second-order (significant
      %       differences). The color key is either on the inside (diagonal) or
      %       the outside. The only time this is necessary is when nCategories<4.
      %       DEFAULT is 'sig'.
      % 
      % SEE ALSO SETOFHYPS.SHOWVALUES, SETOFHYPS.SHOWSIG,
      %    SETOFHYPS.SHOWSIGLEGEND, STATSMAT
         % Parse inputs
         if ~exist('values','var'), values = []; end
         nVals   = numel(values);
         nCats   = numel(self.radii);
         side    = '';
         boxpart = '';
         for v = 1:nargin-2
            if        ischar(varargin{v})
               switch lower(varargin{v})
                  case {'size','left','right'}; side = varargin{v}; sizes = abs(values);
                  case 'sig';             FIRSTORDER = true;
                  case {'sigdiff','diff'};FIRSTORDER = false;
                  case 'color';               colors = mat2cell(1-abs(values(:))*[1 1 1],ones(1,nVals));
                  otherwise                  boxpart = varargin{v};
               end
            elseif isnumeric(varargin{v}), colors = varargin(v);
            elseif    iscell(varargin{v}), colors = varargin{v};
            elseif  ishandle(varargin{v}), ax     = varargin{v};
            end
         end
         if ~exist('ax','var') || isempty(ax), ax = gca; end
         if ~any(values(:)), axtivate(ax); axis off; return; end

         % If whole matrix of values passed, recursively call on each matrix
         % part, after normalizing by abs-max
         if all(size(values)==nCats)
            self.elaborateHinton(values(triu(true(nCats), 1)),'uppert'  ,varargin{:})
            self.elaborateHinton(diag(values)                ,'diagonal',varargin{:})
            self.elaborateHinton(values(tril(true(nCats),-1)),'lowert'  ,varargin{:})
            return
         end
         % Continue input parsing
         if ~exist('sizes' ,'var'),  sizes = ones(nVals,1); end
         if ~exist('colors','var') || isempty(colors)
            colors = mat2cell(zeros(nVals,3),ones(1,nVals));
         end
         % separate values (color) and value signs (circle or square)
         circleNotSquare = double(values>-eps); %values = abs(values); % would set values>=0, but not needed
         if numel(colors)==1, colors = repmat(colors,[1 nVals]); end
         edgewidth = 1;
         % Set up axes
         axtivate(ax); axis ij off;
         set(get(ax,'Parent'),'CurrentAxes',ax,'Units','normalized')

         % Determine size of matrix, save it in figure
         if nCats < 3, nCatsc2   = 0; nCatsc2c2 = 0;
         else          nCatsc2   = nchoosek(nCats,2);
                       nCatsc2c2 = nchoosek(nCatsc2,2);
         end

         % Determine nDiag (number of entries along DIAGonal in the DIAGram)
         SETUP = isempty(ax.UserData);
         if exist('FIRSTORDER','var')
            if FIRSTORDER, nDiag = nCats;
            else           nDiag = nCatsc2;
            end
            ax.UserData = [nDiag FIRSTORDER];
         elseif ~SETUP
            [nDiag,FIRSTORDER] = dealvec(ax.UserData);
         elseif ~isempty(boxpart) % don't use side for this
            switch boxpart
               case {'uppert','lowert'}; nDiag = ceil(sqrt(nVals*2));
               case 'diagonal';          nDiag = nVals;
            end
            FIRSTORDER  = nDiag==nCats; % ambiguous if nCats<4; effectively defaults to true
            ax.UserData = [nDiag FIRSTORDER];
         else
            error('is this a ''sig'' or ''sigdiff''/''diff'' diagram?')
         end
         n  = nDiag;              % shorthand for code readability below
         ix = nchoosek_ix(n,2);

         % Box parameters
         b     = self.boxPos;     % starting coordinates for box containing matrix
         z     = self.boxSize/n;  % size of individual squares in matrix box
         sqScl = 0.97;%0.8;       % scale factor for squares
         csep  = -2*z/3;          % separation between matrix box and circle key

         if SETUP
            %% TEXT LABELS
            if FIRSTORDER
               text('Position',[b+[z*(n-1)/2 z*(n+0.5)]],'HorizontalAlignment','center','String',...
                                                         'Separation','FontSize',12)
               text('Position',[b+[z*(n+0.5) z*(n-1)/2]],'HorizontalAlignment','center','String',...
                                                         'Overlap','Rotation',90,'FontSize',12)
               text('Position',[b+z*n],'String','Radius','Rotation',-45,'FontSize',12)
            else
               text('Position',[b+[z*(n-1)/2 z*(n+0.5)]],'HorizontalAlignment','center','String',...
                                                         'Separation difference','FontSize',12)
               text('Position',[b+[z*(n+0.5) z*(n-1)/2]],'HorizontalAlignment','center','String',...
                                                         'Overlap difference','Rotation',90,'FontSize',12)
               text('Position',[b+z*n],'String','Radius difference','Rotation',-45,'FontSize',12)
            end

            %% COLOR KEY: circles on sides (only for second-order/difference comparisons)
            if ~FIRSTORDER
               cix = nchoosek_ix(nCats,2);
               for i = 1:nCatsc2
                  selfcolors = self.categories.colors(cix(:,i),:);
                  rectangle('Position',[b+[(i-0.5)*z csep] 0.4*[z z]],...
                            'FaceColor',selfcolors(1,:),'Curvature',1,'EdgeColor',selfcolors(1,:))
                  rectangle('Position',[b+[csep (i-0.5)*z] 0.4*[z z]],...
                            'FaceColor',selfcolors(1,:),'Curvature',1,'EdgeColor',selfcolors(1,:))
                  rectangle('Position',[b+[(i-0.9)*z csep] 0.4*[z z]],...
                            'FaceColor',selfcolors(2,:),'Curvature',1,'EdgeColor',selfcolors(2,:))
                  rectangle('Position',[b+[csep (i-0.9)*z] 0.4*[z z]],...
                            'FaceColor',selfcolors(2,:),'Curvature',1,'EdgeColor',selfcolors(2,:))
               end
            end
         end

         %% DRAW ACTUAL CIRCLES/SQUARES
         tt = linspace(-pi/2,pi/2,50)+pi;
         if     strcmpi(side,'left'),   lside =  1;
         elseif strcmpi(side,'right'),  lside = -1;
         end
         if     strcmpi(boxpart,'lowert'), ix = flip(ix);
         elseif strcmpi(boxpart,'diagonal') % plot as circles!
            ix = repmat(1:max(2,nchoosek(n,2)),[2 1]);
            if FIRSTORDER, colors = mat2cell(self.categories.colors,ones(n,1)); end
         end
         for i = 1:nVals
            if ~all(colors{i}==1)
               switch lower(side)
                  case {'left';'right'}
                     if circleNotSquare(i)
                        patchx = -0.5 + ix(1,i) + lside*cos(tt)*sizes(i)*sqScl/2;
                        patchy = -0.5 + ix(2,i) + lside*sin(tt)*sizes(i)*sqScl/2;
                     else
                        patchx = -0.5 + ix(1,i) + [0 -lside -lside 0]*sizes(i)*sqScl/2;
                        patchy = -0.5 + ix(2,i) + [1 1 -1 -1]*sizes(i)*sqScl/2;
                     end
                     patch(b(1)+z*patchx,b(2)+z*patchy,colors{i},...
                           'EdgeColor','None','LineWidth',edgewidth);
                  otherwise % both sides
                     rectangle('Position',[b+z*(ix(:,i)'-(1+sqScl*sizes(i))/2) ...
                                           sizes(i)*sqScl*[z z]],'Curvature',circleNotSquare(i),...
                               'FaceColor',colors{i},'EdgeColor','None','LineWidth',edgewidth)
               end
            end
         end
         axis equal ij off
      end
   end % hidden methods

   methods(Hidden = true, Static = true)
      function sigThresh = nullHypothesisRejected(sig,threshes)
      % SEE ALSO SETOFHYPS.SIGNIFICANCE
         % Convert to sigmas significance, maxed to nSig(=3) and floored
         % Use False Discovery Rate correction for significance computation
         if ~exist('threshes','var') || isempty(threshes), threshes = 0.05; end

         fn = fieldnames(sig)';
         DIFF = ~isempty(sig.ra);

         sigThresh = structfun(@(x) false(1,numel(x)),sig,'UniformOutput',false);
         % Compute increasing levels of significance
         for this_thresh = threshes
            for f = fn; f=f{1};
               if DIFF || strcmpi(f,'di'), t = this_thresh/2;
               else                        t = this_thresh;
               end
               if numel(threshes)>1
                  sigThresh.(f) = sigThresh.(f) + double(fdr_bh(sig.(f),t));
               else
                  sigThresh.(f) = fdr_bh(sig.(f),t);
               end
            end
         end
      end
   end
end


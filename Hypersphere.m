classdef Hypersphere < handle
   properties (SetObservable, GetObservable, AbortSet)
      centers    % [n x d] matrix: n centers of d-dimensional hyperspheres
      radii      % [1 x n] vector: radii of each of n hyperspheres
      categories(1,1) % a Categories object (see help Categories)
   end

   methods
      function obj = Hypersphere(centers,radii,varargin)
      % Contructor for Hypersphere object for hypersphere2sphere, SetOfHyps
      % e.g.:
      % hyp = Hypersphere(centers,radii,<categories>)               % (1)
      % hyp = Hypersphere('estimate',points,categories,<extraargs>) % (2)
      % hyp = Hypersphere(hypset)                                   % (3)
      % hyp = Hypersphere(hypstruct)                                % (4)
      % hyp = Hypersphere({hyps},radii,<extraargs>)                 % (5)
      % 
      % Constructor input options:
      %    (1) centers is an [n x d] numeric matrix, and radii is a [1 x n]
      %       vector, where n is the number of (hyper-)spheres, d is the
      %       number of dimensions. If a Categories object is provided, it
      %       is copied into the Hypersphere object. Otherwise, a dummy
      %       Categories object is automatically generated.
      %       This is the native format of Hypersphere.
      %    (2) If the first argument is 'estimate', then the following inputs
      %       (including the optional <extraargs>) are passed to
      %       estimateHypersphere. That function, as called here, outputs
      %       Hyperpshere object(s). points and categories are required
      %       inputs; other optional valid inputs are the number of bootstraps
      %       and the 'stratified' or 'permute' inputs to determine whether
      %       to sample with or without replacement during bootstrapping.
      %    (3) A SetOfHyps object input is stripped down to a Hypersphere
      %       object.
      %    (4) Mainly for compatibility with old code, this takes a struct
      %       with the same fields as a Hypersphere object and converts it to
      %       Hypersphere object, creating a dummy Categories object if none
      %       exists.
      %    (5) centers is a cell of centers matrices, and radii is a vector
      %       of radii. This recursively calls Hypersphere on the individual
      %       centers cell elements and the corresponding radii vector
      %       elements. This makes sense for bootstrapped center positions.
      % 
      % Properties:
      %    centers    - [n x d] matrix: n centers of d-dimensional hyperspheres
      %    radii      - [1 x n] vector: radii of each of n hyperspheres
      %    categories - a Categories object (see help Categories)
      % 
      % Methods:
      %   Hypersphere.select
      %   Hypersphere.concat
      %   Hypersphere.unconcat
      %   Hypersphere.meanAndMerge
      %   Hypersphere.volume
      %   Hypersphere.dists
      %   Hypersphere.margins
      %   Hypersphere.overlap
      % 
      % 2018-06-07 AZ Created
      % 
      % See also SETOFHYPS, CATEGORIES, ESTIMATEHYPERSPHERE

         % create dummy object (for initializing multiple objects, e.g. option (5))
         if nargin==0; return; end

         if isa(centers,'SetOfHyps')                           % input option (3)
            obj = Hypersphere(centers.centers,centers.radii,...
                              centers.categories);
            return;
         elseif isa(centers,'Hypersphere'), obj = centers; return;
         elseif ischar(centers) && strcmpi(centers,'estimate') % input option (2)
            obj = estimateHypersphere(radii,varargin{:});
            return
         elseif isstruct(centers)                              % input option (4)
            % Helper for older struct-based code
            if isfield('categories',centers)
               obj = Hypersphere(centers.centers,centers.radii,...
                                 centers.categories);
            else
               obj = Hypersphere(centers.centers,centers.radii);
            end
            return
         elseif iscell(centers)                                % input option (5)
            % Recurse if multiple centers provided (e.g., bootstrapping)
            obj = repmat(Hypersphere(),[numel(centers) 1]);
            for i = 1:numel(centers)
               obj(i) = Hypersphere(centers{i},radii(i),varargin{:});
            end
            return
         end

         % input option (1)
         obj.centers = centers;
         obj.radii = radii;

         for v = 1:numel(varargin)
            if isa(varargin{v},'Categories')
               obj.categories = varargin{v};
            end
         end
         % Auto-generate dummy Categories if none exists by now
         if ~isa(obj.categories,'Categories')
            obj.categories = Categories(numel(obj.radii));
         end
      end

      function newobj = select(self,i)
      % Hypersphere.select: outputs a Hypersphere object that has been
      %    subsampled to have one or more hyperspheres, indexed by input i
      % e.g.:
      % fewerHyps = allHyps.select(i)
      %    where i can be a logical vector or list of indices
         if isa(self,'SetOfHyps')
            newobj =   SetOfHyps(self.centers(i,:),...
                                 self.radii(i),    ...
                                 self.categories.select(i));
         else
            newobj = Hypersphere(self.centers(i,:),...
                                 self.radii(i),    ...
                                 self.categories.select(i));
         end
      end

      function self = concat(self)
      % Hypersphere.concat: collapses an array of Hypersphere objects into
      %    a single Hypersphere object, with centers and radii concatenated,
      %    i.e., an [n x 1] array of Hypersphere objects with one hypersphere
      %.   each turns into one Hypersphere object that has n hyperspheres.
      % e.g.:
      % oneBigHyp = aFewHypsArray.concat
         if isa(self,'SetOfHyps')
            self =   SetOfHyps(cat(1,self.centers),[self.radii]);
         else
            self = Hypersphere(cat(1,self.centers),[self.radii]);
         end
      end

      function hyps = unconcat(self)
      % Hypersphere.unconcat: takes a single Hypersphere object with multiple
      %    centers/radii and breaks it out into an array of Hypersphere objects
      %    with a single hypersphere in each. Useful for JSON export
      % e.g.:
      % arrayOfIndividualHyps = hypsallinone.unconcat;
      % 
      % SEE ALSO HYPERSPHERE.MERGE, SETOFHYPS.TOJSON
         n = numel(self.radii);
         for i = 1:n
            hyps(i) = self.select(i);
         end
      end

      function self = meanAndMerge(self)
      % Hypersphere.meanAndMerge: takes an array of Hypersphere (or SetOfHyps)
      %    objects as input. Outputs a single SetOfHyps object, where the
      %    centers and radii are the means of those in the object array. 95%
      %    confidence intervals are computed and their range values are stored
      %    in self.ci.centers and self.ci.radii. The object array originally
      %    input to this function is stored for posterity in self.ci.bootstraps,
      %    and h2s significance tests (from SetOfHyps.significance) are run,
      %    with the results embedded in self.sig and self.sigdiff.
      % e.g.:
      % hypset = arrayOfHyps.meanAndMerge
      % 
      % SEE ALSO SETOFHYPS.SIGNIFICANCE
         ci = [];
         % uses 1st object's categories
         if numel(self)>1 % assumes Hypersphere was bootstrapped
            cperms = ispermuted([self.categories]);
            ci.bootstraps = self(cperms);
            ci.centers = prctile(cat(3,self(cperms).centers),[2.5 97.5],3);
            ci.radii   = prctile(vertcat(self.radii),[2.5 97.5])';
            if all(cperms) % Do we always want this?
               self = Hypersphere(mean(cat(3,self(cperms).centers),3),...
                                  mean(cat(1,self(cperms).radii)),...
                                  self(1).categories);
            else % preserve best (unpermuted) estimate if exists (should be self(1))
               self = self(find(~cperms,1,'first'));
            end
         end
         self = SetOfHyps(self,ci);

         % might as well populate significance fields if we have the bootstraps
         if ~isempty(ci) && isempty(self.sig)
            self = self.significance(ci);
         end
      end

      %% FUNCTIONS TO COMPUTE VOLUME, DISTS, MARGINS, OVERLAPS
      % (they call the hidden static functions below)
      function V = volume(self)
      % Hypersphere.volume: Compute volume(s) of hypersphere(s)
      % e.g.:
      % V = hyp.volume
         if numel(self) > 1
            V = cell2mat_concat(arrayfun(@volume,self,'UniformOutput',false));
            return
         end
         V = self.calcVolume(self.centers,self.radii);
      end

      function D = dists(self)
      % Hypersphere.dists: Compute distance(s) between pair(s) of hyperspheres
      %    along their center-connection lines.
      % e.g.:
      % D = hyp.dists
         if numel(self) > 1
            D = cell2mat_concat(arrayfun(@dists,self,'UniformOutput',false));
            return
         end
         D = self.calcDists(self.centers,self.radii);
      end

      function M = margins(self)
      % Hypersphere.margins: Compute margin(s) between pair(s) of hyperspheres.
      %    Margins are defined as the distance along the center-connection line
      %    between two hyperspheres, minus the radii. Margins are exactly the 
      %    same as negative overlaps.
      % e.g.:
      % M = hyp.margins
         if numel(self) > 1
            M = cell2mat_concat(arrayfun(@margins,self,'UniformOutput',false));
            return
         end
         M = self.calcMargins(self.radii,self.dists);
      end

      function ov = overlap(self)
      % Hypersphere.overlap: Compute overlap(s) between pair(s) of hyperspheres.
      %    Overlaps are defined as the distance along the center-connection line
      %    between two hyperspheres, minus the radii. Overlaps are exactly the 
      %    same as negative margins.
      % e.g.:
      % OV = hyp.overlap
         ov = -self.margins;
      end


      %% SET FUNCTION TO VALIDATE CENTERS (they must be numeric, real, and finite)
      function self = set.centers(self,centers)
         arguments
            self
            centers {mustBeNumeric, mustBeFinite, mustBeReal}
         end
         % Ensure consistent formatting
         if isempty(self.radii)
            self.centers = centers;
         else
            if     size(centers,1) == numel(self.radii)
               self.centers = centers;
            elseif size(centers,2) == numel(self.radii)
               self.centers = centers';
            else
               error('size of radii and centers don''t match')
            end
         end
      end
      %% SET FUNCTION TO VALIDATE RADII (they must be numeric, nonnegative, and
      %%    finite)
      function self = set.radii(self,radii)
         arguments
            self
            radii {mustBeNumeric, mustBeFinite, mustBeNonnegative}
         end
         if ~isempty(radii) && ~isempty(self.centers) ...
            && numel(radii) ~= size(self.centers,1)
            error('size of radii and centers don''t match')
         end
         self.radii = radii(:)';
      end
      %% SET FUNCTION TO VALIDATE CATEGORIES (they must be a Categories class
      %%    object)
      function self = set.categories(self,categories)
         arguments
            self
            categories {mustBeCategoriesClass}
         end
         self.categories = categories;
      end
   end

   methods(Static, Hidden = true)
      % Made these functions static so SetOfHyps can use them
      function V = calcVolume(centers,radii)
      % Compute volumes of each hypersphere
         d = size(centers,2);
         V = (radii.^d).*(pi.^(d/2))./gamma((d/2)+1);
      end

      function D = calcDists(centers,radii)
      % Compute distances between pairs of hyperspheres (along center-connection
      %    lines)
         n = numel(radii);
         D = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               D(i) = sqrt(sum(diff(centers([a b],:)).^2));
            end
         end
      end

      function M = calcMargins(radii,dists)
      % Compute margins between pairs of hyperspheres
         n = numel(radii);
         M = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               M(i) = dists(i) - radii(a) - radii(b);
            end
         end
      end
   end
end

% Helper (validation) function
function mustBeCategoriesClass(c)
   if ~isa(c,'Categories')
      error('categories must be Categories class');
   end
end

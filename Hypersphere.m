classdef Hypersphere < handle
   properties (SetObservable, GetObservable, AbortSet)
      centers    % [n x d] matrix: n centers of d-dimensional hyperspheres
      radii      % [1 x n] vector: radii of each of n hyperspheres
      categories(1,1) % Categories object
   end

   methods
      function obj = Hypersphere(centers,radii,varargin)
         % Contructor for Hypersphere object for hypersphere2sphere, SetOfHyps
         % e.g.:
         % hyp = Hypersphere(centers,radii)
         % hyp = Hypersphere(centers,radii,categories)
         % hyp = Hypersphere('estimate',points,categories)
         % hyp = Hypersphere({hyps})
         % 

         % WRITE ACTUAL DOCUMENTATION HERE

         % Methods:
         % Hypersphere.select
         % Hypersphere.merge
         % Hypersphere.concat
         % 
         % 2018-06-07 AZ Created
         % 
         % See also SETOFHYPS, CATEGORIES

         if nargin==0; return; end

         if isa(centers,'Hypersphere'), obj = centers; return;
         elseif ischar(centers) && strcmpi(centers,'estimate')
            obj = estimateHypersphere(radii,varargin{:});
            return
         elseif isstruct(centers) % Helper for older struct-based code
            obj = Hypersphere(centers.centers,centers.radii);
            return
         end

         obj.centers = centers;
         obj.radii = radii;

         for v = 1:numel(varargin)
            if isa(varargin{v},'Categories')
               obj.categories = varargin{v};
            end
         end
         if ~isa(obj.categories,'Categories')
            % keyboard
            obj.categories = Categories(numel(obj.radii));
         end
      end

      function obj = select(obj,i)
         % Hypersphere.select: outputs a Hypersphere object that has been
         %    subsampled to have one or more hyperspheres, indexed by input i
         % e.g.:
         % fewerhyps = allhyps.select(i)
         %    where i can be a logical vector or list of indices
         % 
         obj.centers    = obj.centers(i,:);
         obj.radii      = obj.radii(i);
         obj.categories = obj.categories.select(i);
      end

      function self = concat(self)
         self = Hypersphere(cat(1,self.centers),[self.radii]);
      end

      %% SET FUNCTIONS TO COMPUTE AND STORE VOLUME, DISTS, MARGINS, OVERLAPS
      function V = volume(self)   % Compute volume of a single hypersphere
         V = self.calcVolume(self.centers,self.radii);
      end

      function D = dists(self)    % Compute distance between two hyperspheres
         if numel(self) > 1
            D = cell2mat_concat(arrayfun(@dists,self,'UniformOutput',false));
            return
         end
         D = self.calcDists(self.centers,self.radii);
      end

      function M = margins(self)
         if numel(self) > 1
            M = cell2mat_concat(arrayfun(@margins,self,'UniformOutput',false));
            return
         end
         M = self.calcMargins(self.radii,self.dists);
      end

      function ov = overlap(self)      % Compute overlap distance of two hyperspheres
         ov = -self.margins;
      end

      function hyps = merge(self)
         hyps = self.mergeToSetOfHyps(self);
      end


      %% SET FUNCTION TO VALIDATE CENTERS
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
      %% SET FUNCTION TO VALIDATE RADII
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
      %% SET FUNCTION TO VALIDATE CATEGORIES
      function self = set.categories(self,categories)
         arguments
            self
            categories {mustBeCategoriesClass}
         end
         self.categories = categories;
      end
   end

   methods(Static, Hidden = true)
      % Made these functions static so SetOfHyps can use them (not sure why they won't 
      %    inherit directly when SetOfHyps < Hypersphere) 
      function V = calcVolume(centers,radii)   % Compute volume of a single hypersphere
         d = size(centers,2);
         V = (radii.^d).*(pi.^(d/2))./gamma((d/2)+1);
      end

      function D = calcDists(centers,radii)    % Compute distance between two hyperspheres
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
         n = numel(radii);
         d = dists;
         M = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               M(i) = d(i) - radii(a) - radii(b);
            end
         end
      end

      function hyps = mergeToSetOfHyps(hyps)
         ci = [];
         % uses 1st object's categories
         if numel(hyps)>1 % assumes Hypersphere was bootstrapped
            cperms = ispermuted([hyps.categories]);
            ci.bootstraps = hyps(cperms);
            ci.centers = prctile(cat(3,hyps(cperms).centers),[2.5 97.5],3);
            ci.radii   = prctile(vertcat(hyps.radii),[2.5 97.5])';
            if all(cperms) % Do we always want this?
               hyps = Hypersphere(mean(cat(3,hyps(cperms).centers),3),...
                                  mean(cat(1,hyps(cperms).radii)),...
                                  hyps(1).categories);
            else % preserve best (unpermuted) estimate if exists (should be hyps(1))
               hyps = hyps(find(~cperms,1,'first'));
            end
         end
         hyps = SetOfHyps(hyps,ci);

         % might as well populate significance fields if we have the bootstraps
         if ~isempty(ci) && isempty(hyps.sig)
            hyps = hyps.significance(ci);
         end
      end
   end
end
function mustBeCategoriesClass(c)
   if ~isa(c,'Categories')
      error('categories must be Categories class');
   end
end


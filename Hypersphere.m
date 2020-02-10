classdef Hypersphere < handle
   properties
      centers    % [n x d] matrix: n centers of d-dimensional hyperspheres
      radii      % [1 x n] vector: radii of each of n hyperspheres
      categories(1,1) % Categories object
   end

   methods
      function obj = Hypersphere(centers,radii,varargin)
         % Contructor for Hypersphere object for hypersphere2sphere, SetOfHyps
         % e.g.:
         % hyp = Hypersphere(centers,radii)
         % 

         % WRITE ACTUAL DOCUMENTATION HERE

         % Methods:
         % Hypersphere.select
         % Hypersphere.merge
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

         % Recurse if multiple centers provided (e.g., bootstrapping)
         if iscell(centers) 
            obj = repmat(Hypersphere(),[numel(centers) 1]);
            for i = 1:numel(centers)
               obj(i) = Hypersphere(centers{i},radii(i));
            end
            return
         end
         % Ensure consistent formatting
         if size(centers,1) ~= numel(radii), obj.centers = centers';
         else                                obj.centers = centers;
         end
         obj.radii = radii(:)';

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

      function obj = concat(self)
         obj = Hypersphere(cat(1,self.centers),[self.radii]);
      end

      %% SET FUNCTIONS TO COMPUTE AND STORE VOLUME, DISTS, MARGINS, OVERLAPS
      function V = volume(self)   % Compute volume of a single hypersphere
         d = size(self.centers,2);
         V = (self.radii.^d).*(pi.^(d/2))./gamma((d/2)+1);
      end
      
      function D = dists(self)    % Compute distance between two hyperspheres
         if numel(self) > 1
            D = cell2mat_concat(arrayfun(@dists,self,'UniformOutput',false));
            return
         end
         n = numel(self.radii);
         D = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               D(i) = sqrt(sum(diff(self.centers([a b],:)).^2));
            end
         end
      end

      function M = margins(self)
         if numel(self) > 1
            M = cell2mat_concat(arrayfun(@margins,self,'UniformOutput',false));
            return
         end
         n = numel(self.radii);
         d = self.dists;
         M = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               M(i) = d(i) - self.radii(a) - self.radii(b);
            end
         end
      end

      function ov = overlap(self)      % Compute overlap distance of two hyperspheres
         ov = -self.margins;
      end

      function hyps = mergeToSetOfHyps(self)
         ci = [];
         % uses 1st object's categories
         if numel(self)>1 % assumes Hypersphere was bootstrapped
            ci.bootstraps = self;
            ci.centers = prctile(cat(3,self.centers),[2.5 97.5],3);
            ci.radii   = prctile(vertcat(self.radii),[2.5 97.5])';
            self       = Hypersphere(mean(cat(3,self.centers),3),...
                                     mean(cat(1,self.radii)),...
                                     self(1).categories);
         end
         hyps = SetOfHyps(self,ci);
      end
   end
end


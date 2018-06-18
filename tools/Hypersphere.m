classdef Hypersphere
   properties
      centers    = [0 0];   % [1 x d] vector: center of hypersphere
      radii      = 1;       % scalar:         radius of hypersphere
   end
   methods
      function obj = Hypersphere(centers,radii)     % Constructor
         if nargin==0; return; end
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
      end
      function obj = select(obj,i)
         obj.centers = obj.centers(i,:);
         obj.radii   = obj.radii(i);
      end
      function obj = merge(obj,ix)
         if ~exist('ix','var'), ix = 1:numel(obj); end
         if isa(ix,'Hypersphere')
            nself = numel(obj);
            n2mrg = numel(ix);
            if nself==n2mrg && n2mrg > 1 % merge parallel arrays into one array
               for i = 1:n2mrg
                  obj(i) = obj(i).merge(ix(i));
               end
            elseif n2mrg > 1             % merge second array into first object
               for i = 1:n2mrg
                  obj = obj.merge(ix(i));
               end
            else                         % merge two objects
               obj.centers = [obj.centers;ix.centers];
               obj.radii   = [obj.radii   ix.radii  ];
            end
            return
         else                            % merge (indexed) array of objects into one object
            obj(1).centers = vertcat(obj(ix).centers);
            obj(1).radii   = horzcat(obj(ix).radii  );
            obj = obj(1);
         end
      end
      function V = volume(obj)          % Compute volume of a single hypersphere
         r = obj.radii;
         d = size(obj.centers,1);
         V = (r.^d).*(pi.^(d/2))./gamma((d/2)+1);
      end
      function D = dists(obj)           % Compute distance between two hyperspheres
         n = numel(obj.radii);
         D = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               D(i) = sqrt(sum(diff(obj.centers([a b],:)).^2));
            end
         end
      end
      function M = margins(obj)
         n = numel(obj.radii);
         d = obj.dists;
         M = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               M(i) = d(i) - obj.radii(a) - obj.radii(b);
            end
         end
      end
      function V = overlap(obj)        % Compute overlap distance of two hyperspheres
         n = numel(obj.radii);
         V = zeros(1,(n^2-n)/2);
         % Recurse if more than two objs passed
         if n > 2, i=0;
            for a = 1:n-1
               for b = a+1:n, i=i+1;
                  V(i) = obj.select([a b]).overlap;
               end
            end
            return
         end
         % Collect dimensionality and radii, compute V if hyperspheres are
         % separate or enclosed
         r = obj.radii;
         V = sum(r)-obj.dists;
         V = min(V,2*min(r));
      end
   end
end


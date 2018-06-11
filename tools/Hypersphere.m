classdef Hypersphere
   properties
      centers    = [0 0];   % [1 x d] vector: center of hypersphere
      radii      = 1;       % scalar:         radius of hypersphere
   end
   methods
      function obj = Hypersphere(centers,radii)     % Constructor
         if nargin==0; return; end
         if size(centers,1) ~= numel(radii), obj.centers = centers';
         else                                obj.centers = centers;
         end
         obj.radii = radii(:)';
      end
      function obj = select(obj,i)
         obj.centers    = obj.centers(i,:);
         obj.radii      = obj.radii(i);
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
      function V = overlap(obj)        % Compute overlap volume of two hyperspheres
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
         d = obj.dists;
         r = obj.radii;
         if     d >= sum(r), V = 0;   return
         elseif d < 1e-6 || d-abs(diff(r)) < 1e-6
            V = obj.select(minix(r)).volume;  return
         else % continue
         end
         % Assume partial overlap of hyperspheres. Compute sector cap heights.
         h = [(r(2)-r(1)+d)*(r(2)+r(1)-d) (r(1)-r(2)+d)*(r(1)+r(2)-d)]/(2*d);
         n = numel(r);
         % Compute overlap (lens) volume by computing each sector cap volume
         % see http://mathworld.wolfram.com/Sphere-SphereIntersection.html
         % and https://en.wikipedia.org/wiki/Spherical_cap#Hyperspherical_cap
         V = 0.5*obj.volume ...
               .*betainc((2*r.*h-h.^2)./(r.^2),(n+1)/2,1/2) ...
               ./beta((n+1)/2,1/2);
         V = sum(V);
      end
   end
end


classdef Hypersphere
   properties
      cntr       % [1 x d] vector: center of hypersphere
      r          % scalar:         radius of hypersphere
      d          % scalar: number of dimensions
   end
   methods
      function obj = Hypersphere(center,radius)     % Constructor
         obj.cntr = center(:)';
         obj.r    = radius;
         obj.d    = numel(center);
      end
      function V = volume(obj)          % Compute volume of a single hypersphere
         V = ([obj.r].^[obj.d]).*(pi.^([obj.d]/2))./gamma(([obj.d]/2)+1);
      end
      function D = dist(objs)           % Compute distance between two hyperspheres
         n = numel(objs);
         D = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               D(i) = sqrt(sum([objs([a b]).cntr].^2));
            end
         end
      end
      function V = overlap(objs)        % Compute overlap volume of two hyperspheres
         n = numel(objs);
         V = zeros(1,(n^2-n)/2);
         % Recurse if more than two objs passed
         if n > 2, i=0;
            for a = 1:n-1
               for b = a+1:n, i=i+1;
                  V(i) = objs([a b]).overlap;
               end
            end
            return
         end
         d = objs.dist;
         if     d >= sum([objs.r]),       V = 0;   return
         elseif d <= abs(diff([objs.r])), V = objs(minix([objs.r])).volume;  return
         else % continue
         end
         % Compute sector cap heights
         h = [(objs(2).r - objs(1).r + d)*(objs(2).r + objs(1).r + d)/(2*d);
              (objs(1).r - objs(2).r + d)*(objs(1).r + objs(2).r + d)/(2*d) ];
         % Compute overlap (lens) volume by computing each sector cap volume
         % see http://mathworld.wolfram.com/Sphere-SphereIntersection.html
         % and https://en.wikipedia.org/wiki/Spherical_cap#Hyperspherical_cap
         V = 0.5*objs.volume ...
               .*betainc((2*[objs.r].*h(:)'-h(:)'.^2)/([objs.r].^2),([objs.d]+1)/2,1/2) ...
               ./beta(([objs.d]+1)/2,1/2);
         V = sum(V);
      end
   end
end


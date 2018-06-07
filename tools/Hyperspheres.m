classdef Hyperspheres < Hypersphere
   properties
      categories = true;
      error      = NaN;
      dists      = NaN;
      margins    = NaN;
   end
   methods
      function obj = Hyperspheres(h,varargin)  % Contructor
         if ~isa(h,'Hypersphere')
            h = Hypersphere(h,varargin{1});
            varargin = varargin(2:end);
         end
         obj.centers = h.centers;
         obj.radii   = h.radii;
         obj.dists   = h.dists;
         obj.margins = h.margins;
         if nargin>3, obj.categories = varargin{1};
         else         obj.categories = Categories(numel(obj.radii));
         end
      end
      function obj = select(obj,i)
         n = numel(obj.radii);
         obj = select@Hypersphere(obj,i);
         obj.categories = obj.categories.select(i);
         % select pairs
         j=0; sel=[]; 
         for a = 1:n-1
            for b = a+1:n, j=j+1;
               if sum(ismember(i,[a b]))==2, sel = [sel j]; end
            end
         end
         obj.dists   = obj.dists(sel);
         obj.margins = obj.margins(sel);
      end
      function err = stress(centers_and_radii,hi)
         if isa(centers_and_radii,'Hyperspheres'), lo = centers_and_radii;
         else lo = Hyperspheres(centers_and_radii(:,2:end),centers_and_radii(:,1));
         end
         err = sum( (hi.dists - lo.dists).^2 + (hi.overlap - lo.overlap).^2 );
      end
      function model = h2s(obj,dimLow)
         if ~exist('dimLow','var'), dimLow = 2; end
         x0  = [obj.radii(:) obj.centers(:,1:dimLow)];
         [fit,err] = fminunc(@(x) stress(x,obj),x0);
         model = Hyperspheres(fit(:,2:end),fit(:,1));
         model.error = err;
      end
   end
end


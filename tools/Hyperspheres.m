classdef Hyperspheres < Hypersphere
   properties
      categories = NaN;
      error      = NaN;
      dists      = NaN;
      margins    = NaN;
   end
   methods
      function obj = Hyperspheres(h,varargin)  % Contructor
         if isstruct(h) % Helper for older struct-based code
            h = Hypersphere(h.centers,h.radii);
         elseif isstr(h) && strcmpi(h,'estimate')
            obj = Hyperspheres(estimateHypersphere(varargin{:}));
            return
         elseif isa(h,'Hypersphere') && numel(h)>1 % merge Hypersphere array
            h(1).centers = vertcat(h.centers);
            h(1).radii   = horzcat(h.radii  );
            h = h(1);
         elseif ~isa(h,'Hypersphere')
            h = Hypersphere(h,varargin{1});
            varargin = varargin(2:end);
         end
         obj.centers    = h.centers;
         obj.radii      = h.radii;
         obj.dists      = h.dists;
         obj.margins    = h.margins;
         obj.categories = Categories(numel(obj.radii));
         for v = 1:numel(varargin)
            if isa(varargin{v},'Categories')
               obj.categories = varargin{v};
            elseif isa(varargin{v},'Hyperspheres')
               obj.error = obj.stress(varargin{v});
            end
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
      function err = stress(centers,hi)
         if isa(centers,'Hyperspheres'), lo = centers;
         else lo = Hyperspheres(centers,hi.radii); % recalculates dists, margins
         end
         err = sum( (hi.dists - lo.dists).^2 + (hi.overlap(true) - lo.overlap(true)).^2 );% + (hi.margins - lo.margins).^2 + (hi.overlap(true) - lo.overlap(true)).^2 );
      end
      function model = h2s(obj,varargin)
      % Optimizes stress of self.centers relative to hi
      % Takes as input: self.h2s(dimLow), self.h2s(hi), self.h2s(hi,dimLow)
      % dimLow defaluts to 2, hi defaults to self
         for v = 2:nargin
            if isa(varargin{v-1},'Hyperspheres')
               hi = varargin{v-1};
            elseif isnumeric(varargin{v-1}) && numel(varargin{v-1})==1
               dimLow = varargin{v-1};
            end
         end
         if ~exist('dimLow','var'), dimLow = 2;   end
         if ~exist('hi'    ,'var'), hi     = obj; end
         % Setup optimization and run
         x0  = mdscale(obj.dists,dimLow);
         %x0  = obj.centers(:,1:dimLow);
         opts = optimoptions('fminunc','Display','off');
         [model.centers,err] = fminunc(@(x) stress(x,hi),x0,opts);
         % Output reduced Hyperspheres model
         model.radii = hi.radii;
         model = Hyperspheres(model);
         model.error = err;
      end
      function show(obj,varargin)
         switch size(obj.centers,2)
            case 2; showCirclesModel(obj,varargin{:});
            case 3; showSpheresModel(obj,varargin{:});
         end
      end
   end
end


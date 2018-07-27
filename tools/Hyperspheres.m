classdef Hyperspheres < Hypersphere
   properties
      categories = NaN;
      error      = NaN;
      dists      = NaN;
      margins    = NaN;
      ci         = [];  % confidence intervals, with raw bootstrap data
   end
   methods
      function obj = Hyperspheres(h,varargin)  % Contructor
         if isstruct(h) % Helper for older struct-based code
            h = Hypersphere(h.centers,h.radii);
         elseif isstr(h) && strcmpi(h,'estimate')
            h = estimateHypersphere(varargin{:});
            obj = Hyperspheres(h(1));
            if numel(h)>1 % assumes Hypersphere was bootstrapped
               obj.ci.bootstraps = h(2:end);
               obj.ci.centers = prctile(cat(3,h(2:end).centers),[2.5 97.5],3);
               obj.ci.radii   = prctile(vertcat(h(2:end).radii),[2.5 97.5])';
            end
            return
         elseif isa(h,'Hypersphere') && numel(h)>1 % convert all to Hyperspheres objects
            % note: use Hypersphere.merge to merge Hypersphere objects
            for i = 1:numel(h)
               obj(i) = Hyperspheres(h(i),varargin{:});
            end
            return
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
      function [errtotal,grad,err] = stress(centers,hi)
         if isa(centers,'Hyperspheres'), lo = centers;
         else lo = Hyperspheres(centers,hi.radii(:));
         end
         fudge     = 1e-4;
         cfactor   = 1;%(numel(lo.radii)-1)/2;
         hoverlap  = max(fudge,hi.overlap);
         err_denom = [hoverlap hi.dists];
         err       = abs(err_denom - [lo.overlap lo.dists]).^8;
         err_denom(abs(err_denom)<fudge) = fudge;
         err       = err./abs(err_denom);
         errtotal  = sum(err);

         % Gradient: partial centers derivative of distances
         [n,d]= size(lo.centers);
         dddc = zeros([n d (n^2-n)/2]);
         i = 0;
         for a = 1:n-1
            for b = a+1:n, i = i+1;
               dddc(a,:,i) = (lo.centers(a,:) - lo.centers(b,:) ) / lo.dists(i);
               dddc(b,:,i) = -dddc(a,:,i);
            end
         end

         % Gradient: partial derivatives of error, relative to distances and overlaps
         dEdd = -8*((hi.dists - lo.dists  ).^7)./hi.dists;
         dEdo = -8*((hoverlap - lo.overlap).^7)./hoverlap;
         % Gradient: put it all together
         grad = sum(dddc.*permute(repmat(dEdd' - dEdo',[1 n d]),[2 3 1]),3);
         %grad = sum(dddc.*permute(repmat(1./hoverlap' - 1./hi.dists',[1 n d]),[2 3 1]),3);

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
         x0   = mdscale(hi.dists,dimLow);
         opts = optimoptions(@fmincon,'Display','iter','OutputFcn',@obj.stressPlotFcn,'TolFun',1e-4,'TolX',1e-4,'SpecifyObjectiveGradient',true);%,'DerivativeCheck','on');
         fit = fmincon(@(x) stress(x,hi),x0,[],[],[],[],[],[],[],opts);
         % Output reduced Hyperspheres model
         model = Hyperspheres(fit,hi.radii(:));
         model.error = model.stress(hi);
      end
      function show(obj,varargin)
         obj.error = mean(obj.error(:));
         switch size(obj.centers,2)
            case 2; showCirclesModel(obj,varargin{:});
            case 3; showSpheresModel(obj,varargin{:});
         end
      end
      function stop = stressPlotFcn(obj,varargin)%x,optimValues,state)
         % unpack inputs
         x           = varargin{1};
         optimValues = varargin{2};
         state       = varargin{3};
         % preliminaries
         n           = size(x,1);
         cols        = colormap('lines');
         % copy hi radii to lo Hyperspheres
         lo          = Hyperspheres(x,obj.radii);
         % Recalculate error from high dimensional Hyperspheres obj
         lo.error    = reshape(undeal(3,@() lo.stress(obj)),(n^2-n)/2,[])';

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
            if size(lo.error,1)>2
            subplot(3,3,6); hold on;title('dists'  )
            plot(optimValues.iteration,lo.error(3,i),'o','Color',cols(i,:));
            end
            if size(lo.error,1)>1
            subplot(3,3,5); hold on;title('margins')
            plot(optimValues.iteration,lo.error(2,i),'o','Color',cols(i,:));
            end
            subplot(3,3,4); hold on;title('overlaps')
            plot(optimValues.iteration,lo.error(1,i),'o','Color',cols(i,:));
         end
         ylabel('Relative error')

         stop = false;
      end
   end
end


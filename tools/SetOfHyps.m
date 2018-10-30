classdef SetOfHyps < Hypersphere
   properties
      categories = NaN;
      error      = NaN;
      dists      = NaN;
      margins    = NaN;
      ci         = [];  % confidence intervals, with raw bootstrap data
   end
   methods
      function obj = SetOfHyps(h,varargin)  % Contructor
         if isa(h,'SetOfHyps'), obj=h; return;
         elseif isstruct(h) % Helper for older struct-based code
            h = Hypersphere(h.centers,h.radii);
         elseif ischar(h) && strcmpi(h,'estimate')
            obj = estimateHypersphere(varargin{:});
            obj = SetOfHyps(obj,varargin{:});%.merge;
            return
         elseif isa(h,'Hypersphere') && numel(h)>1
            % convert all to SetOfHyps objects
            % note: use Hypersphere.merge to merge Hypersphere objects
            for i = 1:numel(h)
               obj(i) = SetOfHyps(h(i),varargin{:});
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
            elseif isa(varargin{v},'SetOfHyps')
               obj.error = obj.stress(varargin{v});
            end
         end
      end
      function obj = merge(self)
         obj = SetOfHyps(self(1));
         if numel(self)>1 % assumes Hypersphere was bootstrapped
            obj.ci.bootstraps = self(2:end);
            obj.ci.centers = prctile(cat(3,self(2:end).centers),[2.5 97.5],3);
            obj.ci.radii   = prctile(vertcat(self(2:end).radii),[2.5 97.5])';
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
      function ovp = permuteoverlaps(obj,points,N)
         % Permutation tests for margin significance:
         %    null hypothesis = full overlap
         if ~exist('N','var') || isempty(N), N=100; end
         n   = numel(obj.radii);
         n2  = nchoosek(n,2);
         ix  = nchoosek_ix(n);
         ovp = NaN(N,n2);
         map = NaN(N,n2);
         for i = 1:n2
            catperm = obj.select(ix(:,i)).categories.permute(N);
            for j = 1:N
               tmpobj   = SetOfHyps('estimate',points,catperm(j));
               % Normalize permuted overlaps by maximum possible overlap
               ovp(j,i) = tmpobj.overlap*0.5/min(tmpobj.radii);
            end
         end
      end
      function sig = significance(obj,points,N)
         if ~exist('N','var') || isempty(N), N=100; end
         overlap_perm = obj.permuteoverlaps(points,N);
         % Compute confidence intervals on bootstrapped overlaps & radii for significance
         % TODO: do we need to do this PAIRWISE???
         boots        = SetOfHyps('estimate',points,obj.categories,N);
         radii_boot   = reshape([boots.ci.bootstraps.radii],[],N)';
         overlap_boot = cell2mat_concat(arrayfun(@overlap,boots.ci.bootstraps,'UniformOutput',false));
         % compute indices of radii corresponding to overlaps
         n  = numel(obj.radii);
         ix = nchoosek_ix(n);
         % normalize overlaps by smaller of two radii
         overlap_meas = obj.overlap./min(obj.radii(ix));
         % compute significance
         for i = 1:nchoosek(n,2)
            % Is negative overlap significantly less than permuted overlaps?
            [~,sig.ma(i)] = ttest(overlap_perm(:,i),-overlap_meas(i),'Tail','left' ,'Alpha',0.95);
         end
         for i = 1:n
            % Is the measured radius significantly greater than zero?
            [~,sig.ra(i)] = ttest(  radii_boot(:,i),0,               'Tail','right','Alpha',0.95);
         end
         % What percentile confidence interval of bootstrapped overlaps contains 0?
         sig.ov = abs(2*(mean(overlap_boot<0)-0.5));
      end
      function [errtotal,grad,err] = stress(centers,hi)
         if isa(centers,'SetOfHyps'), lo = centers;
         else lo = SetOfHyps(centers,hi.radii(:));
         end
         fudge     = 1e-8;
         erro      = ( (hi.overlap - lo.overlap)./hi.overlap ).^2;
         ix        = abs(hi.overlap)<fudge;
         erro(ix)  = 10*abs(hi.overlap(ix) - lo.overlap(ix));
         errd      = (hi.dists/mean(hi.dists) - lo.dists/mean(lo.dists)).^2;
         err       = [errd erro];
         errtotal  = sum(err);

         % Gradient: partial centers derivative of distances
         [n,d]= size(lo.centers);
         nn2  = (n^2-n)/2;
         dddc = zeros([n d nn2]);
         i = 0;
         for a = 1:n-1
            for b = a+1:n, i = i+1;
               dddc(a,:,i) = (lo.centers(a,:) - lo.centers(b,:) ) / lo.dists(i);
               dddc(b,:,i) = -dddc(a,:,i);
            end
         end
         % Gradient: partial derivatives of error, relative to distances and overlaps
         otherlo = sum(lo.dists) - lo.dists;
         dEdd = 2*(nn2^2)*(lo.dists.*otherlo - sum(lo.dists.^2)+lo.dists.^2) / (sum(lo.dists)^3);
         dEdo = -2*(hi.overlap - lo.overlap) ./ (hi.overlap.^2);
         dEdo(ix) = 10*(lo.overlap(ix)-hi.overlap(ix))./abs(lo.overlap(ix)-hi.overlap(ix));
         % Gradient: put it all together
         grad = sum(dddc.*permute(repmat(dEdd' - dEdo',[1 n d]),[2 3 1]),3);

      end
      function model = h2s(obj,varargin)
      % Optimizes stress of self.centers relative to hi
      % Takes as input: self.h2s(dimLow), self.h2s(hi), self.h2s(hi,dimLow)
      % dimLow defaluts to 2, hi defaults to self
         for v = 2:nargin
            if isa(varargin{v-1},'SetOfHyps')
               hi = varargin{v-1};
            elseif isnumeric(varargin{v-1}) && numel(varargin{v-1})==1
               if     varargin{v-1} > 3, nboots = varargin{v-1};
               elseif varargin{v-1} > 0, dimLow = varargin{v-1};
               end
            end
         end
         if ~exist('dimLow','var'), dimLow = 2;   end
         if ~exist('nboots','var'), nboots = 0;   end
         if ~exist('hi'    ,'var'), hi     = obj; end
         % Setup optimization and run
         x0   = mdscale(hi.dists,dimLow);
         opts = optimoptions(@fmincon,'TolFun',1e-4,'TolX',1e-4,'SpecifyObjectiveGradient',true,'Display','off');%,'OutputFcn',@obj.stressPlotFcn);%,'DerivativeCheck','on');
         fit = fmincon(@(x) stress(x,hi),x0,[],[],[],[],[],[],[],opts);
         % Output reduced SetOfHyps model
         model = SetOfHyps(fit,hi.radii(:),hi.categories);
         model.error = model.stress(hi);
      end
      function show(obj,varargin)
         obj.error = mean(obj.error(:));
         switch size(obj.centers,2)
            case 2; showCirclesModel(obj,varargin{:});
            case 3; showSpheresModel(obj,varargin{:});
         end
      end
      function camSettings = cameraCalc(obj)
         if numel(obj)>1 % concatenate center & radii
            obj(1).centers = vertcat(obj.centers);
            obj(1).radii   = [obj.radii];
            obj = obj(1);
         end
         
         % set view angle orthogonal to the PC12 plane of the locations
         eigenvecs = pca(obj.centers);
         if size(eigenvecs,2)<2, normal = [0 0 0]';
         else                    normal = null(eigenvecs(:,1:2)');
         end
         
         if corr(obj.centers*normal,obj.radii')<0
            normal = -normal;
         end
         
         objectsCentroid = mean(obj.centers,1);
         camDist = 50*sqrt(sum(std(obj.centers).^2));
         camSettings.CameraPosition = objectsCentroid'-normal*camDist;
         camSettings.CameraTarget   = objectsCentroid;
      end
      function camera(obj,camSettings)
         if ~exist('camSettings','var') || isempty(camSettings)
            camSettings = obj.cameraCalc;
         end
         for f = fields(camSettings)'
            set(gca,f{1},camSettings.(f{1}));
         end
         % Force lighting to camlight headlight
         delete(findall(gca,'type','light'));
         camlight
         set(findobj(gca,'Type','Surface'),'AmbientStrength',  0,...
                                           'SpecularStrength', 1,...
                                           'DiffuseStrength',  1);
      end
      function stop = stressPlotFcn(obj,varargin)%x,optimValues,state)
         % unpack inputs
         x           = varargin{1};
         optimValues = varargin{2};
         state       = varargin{3};
         % preliminaries
         n           = size(x,1);
         cols        = colormap('lines');
         % copy hi radii to lo SetOfHyps
         lo          = SetOfHyps(x,obj.radii);
         % Recalculate error from high dimensional SetOfHyps obj
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
      function showSig(~,self,ax)
         % Parse inputs
         if isa(self,'SetOfHyps'), sig = self.sig;
         else                      sig = self;
         end
         if ~exist('ax','var') || isempty(ax), ax = gca; end

         % Reshape into matrix
         sig.ov = sig.ov <= 0.95;
         mat = statsmat(1-sig.ma,sig.ov,1-sig.ra);

         % Time to get a-plottin'
         set(0,'CurrentFigure',get(ax,'Parent'))
         set(get(ax,'Parent'),'CurrentAxes',ax,'Units','normalized')
         %imagesc(mat)

         n = size(mat,1); ix = nchoosek_ix(n,2); N = size(ix,2);
         % Box parameters
         boxSize = 0.3;
         boxPos  = [0.65 0.65];
         sqSz    = boxSize/n;

         col = [1 0.5 0.5];
         for i = 1:N
            rectangle('Position',[boxPos+sqSz*(ix(:,i)'-1) sqSz sqSz],'Curvature',1,'FaceColor',[1 0 0]*(1-sig.ma(i)),'EdgeColor','k')
            rectangle('Position',[boxPos+sqSz*(fliplr(ix(:,i)')-1) sqSz sqSz],'Curvature',1,'FaceColor',[0 0 1]*sig.ov(i),'EdgeColor','k')
         end
         for i = 1:n
            rectangle('Position',[boxPos+[sqSz sqSz]*(i-1) sqSz sqSz],'Curvature',1,'FaceColor',[1 1 1]*(1-sig.ra(i)),'EdgeColor','k')
         end
      end
   end
end


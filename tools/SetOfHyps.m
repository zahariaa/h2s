classdef SetOfHyps < Hypersphere
   properties
      categories = NaN;
      error      = NaN;
      dists      = NaN;
      margins    = NaN;
      ci         = [];  % confidence intervals, with raw bootstrap data
   end
   properties(GetAccess = 'private', SetAccess = 'private')
      boxPos     = [0.65 0.65];
      boxSize    = 0.3;
      nSigLevs   = 3; % 3 means [0.95 0.99 0.999] levels
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
            obj.ci.bootstraps = self;
            obj.ci.centers = prctile(cat(3,self.centers),[2.5 97.5],3);
            obj.ci.radii   = prctile(vertcat(self.radii),[2.5 97.5])';
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
         for i = 1:n2
            catperm = obj.select(ix(:,i)).categories.permute(N);
            for j = 1:N
               tmpobj   = SetOfHyps('estimate',points,catperm(j));
               % Normalize permuted overlaps by maximum possible overlap
               ovp(j,i) = tmpobj.overlap*0.5/min(tmpobj.radii);
            end
         end
      end
      function [sig,sec] = significance(obj,points,N)
         if ~exist('N','var') || isempty(N), N=100; end
         % Compute confidence intervals on bootstrapped overlaps & radii for significance
         boots        = SetOfHyps('estimate',points,obj.categories,N,'stratified').merge;
         radii_boot   = reshape([boots.ci.bootstraps.radii],[],N)';
         margin_boot  =  vertcat(boots.ci.bootstraps.margins);
         overlap_boot = cell2mat_concat(arrayfun(@overlap,boots.ci.bootstraps,'UniformOutput',false));
%% compute significance
         ciprctileLtail = @(x)        mean(x>0);
         ciprctile2tail = @(x) abs(2*(mean(x<0)-0.5));
         % What percentile confidence interval of bootstrapped margins contains 0?
         sig.ma = ciprctileLtail(margin_boot);
         % What percentile confidence interval of bootstrapped overlaps contains 0?
         sig.ov = ciprctileLtail(overlap_boot);
         sig.ra = [];

%% SECOND-ORDER COMPARISONS
         % compute indices of radii corresponding to overlaps
         n     = numel(obj.radii);
         nc2   = nchoosek(n,2);
         nc2c2 = nchoosek(nc2,2);
         ix    = nchoosek_ix(n);
         ixc2  = nchoosek_ix(nc2);
         %dists_boot = cell2mat_concat(arrayfun(@dists,boots.ci.bootstraps,'UniformOutput',false));

         sec = struct('ra',NaN(1,nc2),'ma',NaN(1,nc2c2),'ov',NaN(1,nc2c2));
         for i = 1:nc2
            sec.ra(i) = ciprctile2tail(diff(  radii_boot(:,  ix(:,i)),[],2));
         end
         for i = 1:nc2c2
            sec.ma(i) = ciprctile2tail(diff( margin_boot(:,ixc2(:,i)),[],2));
            sec.ov(i) = ciprctile2tail(diff(overlap_boot(:,ixc2(:,i)),[],2));
         end
      end
      function [errtotal,grad,err] = stress(centers_and_radii,hi)
         if isa(centers_and_radii,'SetOfHyps'), lo = centers_and_radii;
         else lo = SetOfHyps(centers_and_radii(:,2:end),centers_and_radii(:,1));
         end
         fudge     = 1e-8;
         erro      = ( hi.overlap - lo.overlap ).^2;
         errd      = ( hi.dists   - lo.dists   ).^2;
         erra      = ( hi.radii   - lo.radii   ).^2;
         err       = [errd erro erra];
         errtotal  = sum(err);


      end
      function model = h2s(obj,varargin)
      % Optimizes stress of self.centers relative to hi
      % Takes as input: self.h2s(dimLow), self.h2s(hi), self.h2s(hi,dimLow)
      % dimLow defaluts to 2, hi defaults to self
         for v = 2:nargin
            if isa(varargin{v-1},'SetOfHyps')
               hi = varargin{v-1};
               lo = obj;
            elseif isnumeric(varargin{v-1}) && numel(varargin{v-1})==1
               if     varargin{v-1} > 3, nboots = varargin{v-1};
               elseif varargin{v-1} > 0, dimLow = varargin{v-1};
               end
            end
         end
         if ~exist('dimLow','var'), dimLow = 2;   end
         if ~exist('nboots','var'), nboots = 0;   end
         if ~exist('hi'    ,'var'), hi = obj; lo = obj; end

         % Recurse
         n = numel(hi);
         if n > 1
            for i = 1:n
               stationarycounter(i,n);
               model(i) = lo(i).h2s(dimLow,nboots);
            end
            return
         end
         
         % Setup optimization and run
         x0  = [hi.radii(:) mdscale(hi.dists,dimLow)];
         opts = optimoptions(@fmincon,'TolFun',1e-4,'TolX',1e-4,'Display','iter');%,'SpecifyObjectiveGradient',true,'Display','off');%,'OutputFcn',@obj.stressPlotFcn);%,'DerivativeCheck','on');
         nonlcon = @(x) hi.constrain_pos_overlaps(x);
         fit = fmincon(@(x) stress(x,hi),x0,[],[],[],[],[],[],[],opts);
         % Output reduced SetOfHyps model
         model = SetOfHyps(fit(:,2:end),fit(:,1),hi.categories);
         model.error = model.stress(hi);
      end
      function [c,ceq] = constrain_pos_overlaps(self,x)
         new = SetOfHyps(struct('centers',x,'radii',self.radii));
         % Any negative margin (positive overlap) must stay a non-positive margin (non-negative overlap).
         c   = new.margins(self.margins<1000*eps);
         ceq = [];
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
      function movie(self,times,timeLabel,varargin)
      % function SetOfHyps.movie(times,timeLabel,<SAVEflag> or <StaticFrameAxes>)
      % SetOfHyps.movie(times,timeLabel,SAVEflag) generates a movie, optionally saves
      % SetOfHyps.movie(times,timeLabel,StaticFrameAxes) plots static frames in axes given
         nh     = numel(self);
         STATIC = false;
         if nargin > 3
            if        islogical(varargin{1}) , SAVE = varargin{1};
            elseif all(ishandle(varargin{1})), SAVE = false; STATIC = true;
                                               ax   = varargin{1};
            end
         else                                  SAVE = false;
         end
         if ~exist('ax'       ,'var'), ax        = repmat(gca,[1 nh]);   end
         if ~exist('timeLabel','var'), timeLabel = [];    end
         % Calculate frame times
         if exist('times','var') && ~isempty(times)
            if isempty(timeLabel), timeLabel = 'ms'; end
         else
            times = 1:nh;
         end
         dimLow = size(self(1).centers,2);
         % Calculate axis bounds
         axbounds = NaN(1,dimLow*2);
         for i = 1:nh
            for j = 1:dimLow
               axbounds(1+(j-1)*2) = min([axbounds(1+(j-1)*2);
                                    self(i).centers(:,j)-self(i).radii']);
               axbounds(   j   *2) = max([axbounds(   j   *2);
                                    self(i).centers(:,j)+self(i).radii']);
            end
         end
         % Calculate camera view
         if dimLow == 3, camSettings = self.cameraCalc; end
         % Legend text position
         txtY = [NaN;1.075;0.89];

         if SAVE
            fps = 20;
            vidObj           = VideoWriter([SAVE '.avi']);
            vidObj.FrameRate = fps;
            vidObj.Quality   = 100;
            open(vidObj);
         end
         %% Do stuff
         ann = [];
         for i = 1:nh
            axtivate(ax(i));
            % DRAW PLOT
            self(i).show(false,[]);
            if dimLow == 3,   self.camera(camSettings);   end
            axis(axbounds)
         
            %% Generate title string
            extratxt = sprintf('} %g %s',times(i),timeLabel);
            if STATIC
               title(extratxt(3:end));
               extratxt = [];
            end
            if i==1 || ~STATIC
               delete(ann);
               ann = self(1).categories.legend([0.01 txtY(dimLow) 1 0.1],extratxt);
            end
            if ~SAVE
               drawnow
               pause(0.1)
            else
               writeVideo(vidObj,getframe(gcf));
            end
         end
         if STATIC
            subplotResize([],[],0.01)
         elseif SAVE
            close(vidObj)
            fprintf('\nVideo written to %s\n',SAVE)
         end
      end
      function plotComparisons(self,type,times,ax)
         if ~exist('ax','var') || isempty(ax), ax = gca; end 
         n  = numel(self(1).radii);
         nk = nchoosek(n,2);
         ix = nchoosek_ix(n)';

         axtivate(ax);
         set(ax,'ColorOrder',self(1).categories.colors(ix(:),:))
         plot(times,reshape([self.(type)],nk,[])','LineWidth',1)
         plot(times,reshape([self.(type)],nk,[])','--','LineWidth',1)
         ylabel(type)
         xlim([min(times) max(times)]);
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
      function showSig(obj,self,ax)
         % Parse inputs
         if isa(self,'SetOfHyps'), sig = self.sig;
         else                      sig = self;
         end
         if ~exist('ax','var') || isempty(ax), ax = gca; axis ij off; end

         nSigLevs = obj.nSigLevs; % 3 means [0.95 0.99 0.999] levels

         % Convert to sigmas significance, maxed to nSig(=3) and floored
         %sig = structfun(@(x) erfinv(x)*sqrt(2),sig,'UniformOutput',false);
         sigThresh = structfun(@(x) x>0.95,sig,'UniformOutput',false);
         for i = 2:nSigLevs
            for f = fieldnames(sig)'; f=f{1};
               sigThresh.(f) = sigThresh.(f) + double(sig.(f)>(1-10^-i));
            end
         end
         sigThresh = structfun(@(x) x/nSigLevs,sigThresh,'UniformOutput',false);
         % Reshape into matrix
         mat = statsmat(sigThresh.ma,sigThresh.ov,sigThresh.ra);

         % Time to get a-plottin'
         set(0,'CurrentFigure',get(ax,'Parent'))
         set(get(ax,'Parent'),'CurrentAxes',ax,'Units','normalized')
         %imagesc(mat)

         n = size(mat,1); ix = nchoosek_ix(n,2); N = size(ix,2);
         % Box parameters
         boxPos = obj.boxPos;     % starting coordinates for box containing matrix
         sqSz   = obj.boxSize/n;  % size of individual squares in matrix box
         sqScl  = 0.8;            % scale factor for squares
         csep   = -2*sqSz/3;      % separation between matrix box and circle key

         %% DRAW UNDER-BOX OVERLAP/MARGIN SIGNIFIER AREAS
         upperX = [vectify(repmat(1:n,2,[]));1]*sqSz;
         plot(boxPos(1)+upperX,boxPos(2)+shift(upperX-sqSz,-1),'k-' )

         %% DRAW ACTUAL BOXES
         for i = 1:N
            rectangle('Position',[boxPos+sqSz*(ix(:,i)'-1+(1-sqScl)/2) sqScl*[sqSz sqSz]],...
                      'Curvature',0.2,'FaceColor',[1 1 1]*(1-sigThresh.ma(i)),'EdgeColor','k')
         end
         if isempty(sig.ra)
            plot(boxPos(1)+shift(upperX-sqSz,-1),boxPos(2)+upperX,'k--')
            for i = 1:N
            rectangle('Position',[boxPos+sqSz*(fliplr(ix(:,i)')-1+(1-sqScl)/2) sqScl*[sqSz sqSz]],...
                      'Curvature',0.2,'FaceColor',[1 1 1]*(1-sigThresh.ov(i)),'EdgeColor','k')
            end
            text('Position',[boxPos+[sqSz*(n-1)/2 sqSz*(n+0.1)]],'HorizontalAlignment','center','String','Overlap')
            text('Position',[boxPos+[sqSz*(n+0.1) sqSz*(n-1)/2]],'HorizontalAlignment','center','String','Margin','Rotation',90)
         else
            for i = 1:n
            rectangle('Position',[boxPos+[sqSz sqSz]*(i-1+(1-sqScl)/2) sqScl*[sqSz sqSz]],...
                      'Curvature',0.2,'FaceColor',[1 1 1]*(1-sigThresh.ra(i)),'EdgeColor','k')
            end
            text('Position',[boxPos+[sqSz*(n+0.1) sqSz*(n-1)/2]],'HorizontalAlignment','center','String','Margin/overlap difference','Rotation',90)
            text('Position',[boxPos+sqSz*n],'String','Radius difference','Rotation',-45)
         end

         %% COLOR KEY: circles
         nCats   = size(obj.categories.colors,1);
         nCatsc2 = nchoosek(nCats,2);
         if isempty(sig.ra)
            for i = 1:n
               rectangle('Position',[boxPos+[(i-0.75)*sqSz csep] 0.5*[sqSz sqSz]],...
                         'FaceColor',obj.categories.colors(i,:),'Curvature',1,...
                         'EdgeColor',obj.categories.colors(i,:))
               rectangle('Position',[boxPos+[csep (i-0.75)*sqSz] 0.5*[sqSz sqSz]],...
                         'FaceColor',obj.categories.colors(i,:),'Curvature',1,...
                         'EdgeColor',obj.categories.colors(i,:))
            end
         else
            ix = nchoosek_ix(nCats,2);
            for i = 1:nCatsc2
               rectangle('Position',[boxPos+[(i-0.5)*sqSz csep] 0.4*[sqSz sqSz]],...
                         'FaceColor',obj.categories.colors(ix(1,i),:),'Curvature',1,...
                         'EdgeColor',obj.categories.colors(ix(1,i),:))
               rectangle('Position',[boxPos+[csep (i-0.5)*sqSz] 0.4*[sqSz sqSz]],...
                         'FaceColor',obj.categories.colors(ix(1,i),:),'Curvature',1,...
                         'EdgeColor',obj.categories.colors(ix(1,i),:))
               rectangle('Position',[boxPos+[(i-0.9)*sqSz csep] 0.4*[sqSz sqSz]],...
                         'FaceColor',obj.categories.colors(ix(2,i),:),'Curvature',1,...
                         'EdgeColor',obj.categories.colors(ix(2,i),:))
               rectangle('Position',[boxPos+[csep (i-0.9)*sqSz] 0.4*[sqSz sqSz]],...
                         'FaceColor',obj.categories.colors(ix(2,i),:),'Curvature',1,...
                         'EdgeColor',obj.categories.colors(ix(2,i),:))
            end
         end
         axis equal ij off
      end
      function showSigLegend(self,ax)
         n        = numel(self.radii);
         boxPos   = self.boxPos;
         sqSz     = self.boxSize/n;
         nSigLevs = self.nSigLevs; % 3 means [0.95 0.99 0.999] levels
         if ~exist('ax','var') || isempty(ax), ax = gca; axis ij off; end
         % Time to get a-plottin'
         set(0,'CurrentFigure',get(ax,'Parent'))
         set(get(ax,'Parent'),'CurrentAxes',ax,'Units','normalized')
%% Draw Legend
         for i = 1:nSigLevs
            rectangle('Position',[boxPos+[0.01+sqSz*(n+1/3) sqSz*(i-1)/3] sqSz/3 sqSz/3],...
                      'FaceColor',[1 1 1]*(1-i/nSigLevs),'EdgeColor','none')
         end
         % Proability labels
            text('Position',[boxPos+[0.01+sqSz*(n+0.7) sqSz/6]],'String','p<0.05')
         for i = 2:nSigLevs
            text('Position',[boxPos+[0.01+sqSz*(n+0.7) sqSz*(2*i-1)/6]],'String',...
                 ['p<' num2str(10^-i)])
         end
         axis equal ij off
      end
   end
end


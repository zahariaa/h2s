classdef SetOfHyps < Hypersphere
   properties
      categories = NaN;
      error      = NaN;
      dists      = NaN;
      margins    = NaN;
      ci         = [];  % confidence intervals, with raw bootstrap data
      sig        = [];  % significance tests results
      msflips    = [];  % margin sign flips
   end
   properties(GetAccess = 'private', SetAccess = 'private')
      boxPos     = [0.65 0.65];
      boxSize    = 0.3;
      nSigLevs   = 3; % 3 means [0.95 0.99 0.999] levels
   end

   methods
      function obj = SetOfHyps(h,varargin)  % Contructor
         if ~exist('h','var'), return; % i don't understand why this line is necessary
         elseif isa(h,'SetOfHyps'), obj=h; return;
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
               [~,~,obj.error,obj.msflips] = obj.stress(varargin{v});
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
         dist_boot    =  vertcat(boots.ci.bootstraps.dists);
         radii_boot   = reshape([boots.ci.bootstraps.radii],[],N)';
         margin_boot  =  vertcat(boots.ci.bootstraps.margins);
         %overlap_boot = cell2mat_concat(arrayfun(@overlap,boots.ci.bootstraps,'UniformOutput',false));
%% compute significance
         ciprctileLtail = @(x)        mean(x>0);
         ciprctile2tail = @(x) abs(2*(mean(x<0)-0.5));
         %% TODO: DELETE margin_boot/sig.ma?
         % What percentile confidence interval of bootstrapped margins contains 0?
         sig.ma = ciprctileLtail(margin_boot);
         % What percentile confidence interval of bootstrapped overlap/margins contains 0?
         sig.ov = ciprctileLtail(-margin_boot);
         sig.ra = [];
         % What percentile confidence interval of bootstrapped distances contains 0?
         sig.di = ciprctileLtail(dist_boot);

%% SECOND-ORDER COMPARISONS
         % compute indices of radii corresponding to overlaps
         n     = numel(obj.radii);
         nc2   = nchoosek(n,2);
         nc2c2 = nchoosek(nc2,2);
         ix    = nchoosek_ix(n);
         ixc2  = nchoosek_ix(nc2);

         sec = struct('ra',NaN(1,nc2),'ma',NaN(1,nc2c2),'ov',NaN(1,nc2c2),'di',NaN(1,nc2c2));
         for i = 1:nc2
            sec.ra(i) = ciprctile2tail(diff(  radii_boot(:,  ix(:,i)),[],2));
         end
         for i = 1:nc2c2
            sec.ov(i) = ciprctile2tail(diff(-margin_boot(:,ixc2(:,i)),[],2));
            sec.di(i) = ciprctile2tail(diff(   dist_boot(:,ixc2(:,i)),[],2));
         end
      end

      function [errtotal,grad,err,msflips] = stress(centers_and_radii,hi)
         if isa(centers_and_radii,'SetOfHyps'), lo = centers_and_radii;
         else lo = SetOfHyps(centers_and_radii(:,2:end),centers_and_radii(:,1));
         end
         [nc,dimLow]= size(lo.centers);
         fudge     = 1e-8;
         errd      = hi.dists   - lo.dists;
         erro      = hi.margins - lo.margins;
         erra      = hi.radii   - lo.radii;
         err       = abs([errd erro erra]);
         errtotal  = mean( err.^2 );
         %errtotal  = mean( mean(errd.^2) + mean(erro.^2) + mean(erra.^2) );

         msflips = ~(sign(hi.margins) == sign(lo.margins));
         % Gradient: partial centers derivative of distances (should be size of centers)
         % Gradient: partial derivatives of error, relative to distances and overlaps
         ix = nchoosek_ix(nc);
         nc2 = size(ix,2);
         
         % center gradients
         dddc = zeros(nc2,nc,dimLow);% 1./lo.dists;
         centersdiff = lo.centers(ix(2,:),:) - lo.centers(ix(1,:),:);
         for i = 1:nc2
            dddc(i,ix(:,i),:) = dddc(i,ix(:,i),:) ...
                                + shiftdim([1;-1]*centersdiff(i,:),-1);
         end
         dddc = dddc./repmat(lo.dists',[1 nc dimLow]);
         dEdd = 2*( errd + erro );
         dEdc = reshape(dEdd*reshape(dddc,nc2,[]),[nc dimLow]);

         % radius gradients
         dEdr = -2*erra';
         for i = 1:nc
            dEdr(i) =  dEdr(i) + 2*sum(erro(~~sum(ix==i)));
         end
         % Gradient: put it all together
         grad = [dEdr dEdc]/(2*nc2+nc);
      end

      function model = h2s(obj,varargin)
      % Optimizes stress of self.centers relative to hi
      % Takes as input: self.h2s(dimLow), self.h2s(hi), self.h2s(hi,dimLow)
      % dimLow defaluts to 2, hi defaults to self
      % also can do self.h2s(true) to fix radii to hi dimensional ones
      % self.h2s([true true true]) fixes radii, constrains high positive
      %    overlaps to be positive, and uses MDS for low centers initialization
         for v = 2:nargin
            if isa(varargin{v-1},'SetOfHyps')
               hi = varargin{v-1};
               lo = obj;
            elseif isnumeric(varargin{v-1})
               switch numel(varargin{v-1})
                  case 1;    dimLow                =         varargin{v-1};
                  case 2;   [dimLow,nboots]        = dealvec(varargin{v-1});
                  case 3;   [dimLow,nboots,ninits] = dealvec(varargin{v-1});
                             varargin{v-1}(3) = 1; % change to 1 for recusive call below
               end
            elseif islogical(varargin{v-1})
               switch numel(varargin{v-1})
                  case 1;    FIXRADII           =         varargin{v-1};
                  case 2;   [FIXRADII,MDS_INIT] = dealvec(varargin{v-1});
               end
            end
         end
         if ~exist('hi'    ,'var'), hi = obj; lo = obj; end
         n  = numel(hi);
         [nr,d] = size(hi.centers);
         if ~exist('dimLow','var'), dimLow = min(3,min(nr,d)); end
         if ~exist('nboots','var'), nboots = 0;   end
         if ~exist('ninits','var'), ninits = 1;   end
         if ~exist('MDS_INIT','var'), MDS_INIT = true; end
         if ~exist('FIXRADII','var'), FIXRADII = true; end

         % Recurse for multiple SetOfHyps objects
         if n > 1
            for i = 1:n
               stationarycounter(i,n);
               model(i) = lo(i).h2s(varargin);
            end
            return
         end

         % Setup optimization and run
         if MDS_INIT
            x0  = [hi.radii(:) mdscale(hi.dists,dimLow,'Criterion','metricstress')];
         else
            cminmax = [min(hi.centers(:)) max(hi.centers(:))];
            %% initialization of centers: random noise added to projection onto first 2 dimensions
            x0  = [hi.radii(:) -0.1*diff(cminmax)*rand(nr,dimLow) + hi.centers(:,1:dimLow)];
         end
         if FIXRADII
            fit = x0;
         else
   %          x0  = x0 - repmat(mean(x0),numel(hi.radii),1);
   %          x0  = x0*2;
            opts = optimoptions(@fminunc,'TolFun',1e-4,'TolX',1e-4,'Display','iter','SpecifyObjectiveGradient',true);%,'DerivativeCheck','on');%'Display','off');%,'OutputFcn',@obj.stressPlotFcn);
            fit = fminunc(@(x) stress(x,hi),x0,opts);
         end
         % Output reduced SetOfHyps model
         model = SetOfHyps(fit(:,2:end),fit(:,1),hi.categories);
         [~,~,model.error,model.msflips] = model.stress(hi);

         % Recurse for multiple random initializations
         for iinit = 2:ninits
            tmp = obj.h2s(varargin);
            if max(tmp.error) < max(model.error)
               model = tmp;
            end
         end
      end

      function varargout = show(obj,varargin)
         maxerror = max(obj.error);
         switch size(obj.centers,2)
            case 2; [~,XY] = showCirclesModel(obj,varargin{:});
               %% Draw max error bar
               maxXY = max(XY);
               maxline = plot(maxXY([1 1]) - [0 maxerror],...
                              maxXY([2 2]),'k-','LineWidth',4);
               axis ij tight equal vis3d off % puts error bar at bottom right
            case 3; [~,XY] = showSpheresModel(obj,varargin{:});
               %% Draw max error bar
               axis ij tight equal vis3d off % puts error bar at bottom right
               ax = gca;
               ax.Units = 'normalized';
               maxline = annotation('line',[1-maxerror/diff(ax.XLim) 1]*ax.Position(3)+ax.Position(1),...
                                           ax.Position([2 2]),'LineWidth',4);
         end
         %% DRAW SIGN FLIP ERROR LINES
         ovlines = obj.plotOverlapErrors(obj.sig.ov<0.05);

         if nargout > 0,   varargout = {maxline,ovlines};   end
      end

      function errlines = plotOverlapErrors(obj,sigov)
         if ~exist('sigov','var')
            sigov = true(1,numel(obj.margins));
         end
         if ~isempty(obj.msflips)
            sigov = sigov .* obj.msflips;
         end
         if isempty(sigov) || ~any(sigov)
            errlines = [];
            return
         end
         n = numel(obj.radii);
         ix = nchoosek_ix(n);
         for i = find(sigov)
            err      = obj.error(n+i); % assumes order of errors is radii, margins, distances
            centers  = obj.centers(ix(:,i),:);
            dxy      = diff(centers);
            dxy      = dxy/norm(dxy);
            midpoint = mean(centers);
            % Compute overlap midpoint
            % need to assign +/-1 to larger/smaller center coordinate?
            ov_edges = centers + [1;-1].*obj.radii(ix(:,i))'*dxy;
            ov_midpoint = mean(ov_edges);
            % % debug
            % keyboard
            % plot(centers(:,1),centers(:,2),'ko')
            % plot(midpoint(:,1),midpoint(:,2),'ro')
            % plot(ov_edges(:,1),ov_edges(:,2),'go')
            % plot(ov_midpoint(:,1),ov_midpoint(:,2),'rx')

            if size(obj.centers,2)==2
               errlines(i) = plot( ov_midpoint(1) + [-1 1]*dxy(1)*err/2,...
                                   ov_midpoint(2) + [-1 1]*dxy(2)*err/2,'k-','LineWidth',2);
            else
               errlines(i) = plot3(ov_midpoint(1) + [-1 1]*dxy(1)*err/2,...
                                   ov_midpoint(2) + [-1 1]*dxy(2)*err/2,...
                                   ov_midpoint(3) + [-1 1]*dxy(3)*err/2,'k-','LineWidth',2);
            end
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
         lo.error    = lo.error.^2;

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

      function shapesInABox(obj,values,varargin)
         % Parse inputs
         side = '';
         for v = 1:nargin-2
            if        ischar(varargin{v})
               switch lower(varargin{v})
                  case {'size','left','right'}; side = varargin{v}; sizes = abs(values);
                  otherwise                  boxpart = varargin{v};
               end
            elseif isnumeric(varargin{v}), edgecolors= varargin(v);
            elseif    iscell(varargin{v}), edgecolors= varargin{v};
            elseif  ishandle(varargin{v}), ax        = varargin{v};
            end
         end
         if ~exist('values','var'), values = []; end
         N = numel(values);
         % If whole matrix of values passed, recursively call on each matrix part, after normalizing by abs-max
         if ~any(values), return
         elseif size(values,1)==sqrt(N) && size(values,2)==sqrt(N)
            obj.shapesInABox(values(triu(true(sqrt(N)), 1)),'uppert'  ,varargin{:})
            obj.shapesInABox(diag(values)                  ,'diagonal',varargin{:})
            obj.shapesInABox(values(tril(true(sqrt(N)),-1)),'lowert'  ,varargin{:})
            return
         end
         % Continue input parsing
         if ~exist('sizes' ,'var'),  sizes = ones(N,1); end
         if ~exist('ax','var') || isempty(ax), ax = gca; axis ij off; end
         if ~exist('edgecolors','var') || isempty(edgecolors)
            edgecolors = repmat({'None'},ones(1,N)); edgewidth = 1;
            colors = mat2cell(zeros(N,3),ones(1,N));
         else
            edgewidth = 2; colors = repmat({'None'},[1 N]);
         end
         % separate values (color) and value signs (circle or square)
         circleNotSquare = double(values>-eps); %values = abs(values); % would set values>=0, but not needed
         if numel(edgecolors)==1, edgecolors = repmat(edgecolors,[1 N]); end
         % Set up axes
         set(0,'CurrentFigure',get(ax,'Parent'))
         set(get(ax,'Parent'),'CurrentAxes',ax,'Units','normalized')

         % Determine size of matrix, save it in figure
         nCats   = size(obj.categories.colors,1);
         nCatsc2 = nchoosek(nCats,2);
         if ~isempty(ax.UserData)
            SETUP = false;
            [n,FIRSTORDER] = dealvec(ax.UserData);
         else
            SETUP = true;
            switch N
               case nCatsc2;             FIRSTORDER = true;  n = nCats;  % margins 1st order
               case nchoosek(nCatsc2,2); FIRSTORDER = false; n = nCatsc2;% ma/ov diff 2nd order
               case nCatsc2^2;           FIRSTORDER = false; n = nCatsc2;% full matrix 2nd order
               case 0;                   FIRSTORDER = true;  n = nCats;  % (empty) radii 1st order
               case nCats;               FIRSTORDER = true;  n = nCats;  % radii diff 2nd order
               case nCats^2;             FIRSTORDER = true;  n = nCats;  % full matrix 1st order
            end
            ax.UserData = [n FIRSTORDER];
         end
         ix = nchoosek_ix(n,2);

         % Box parameters
         boxPos = obj.boxPos;     % starting coordinates for box containing matrix
         sqSz   = obj.boxSize/n;  % size of individual squares in matrix box
         sqScl  = 0.97;%0.8;         % scale factor for squares
         csep   = -2*sqSz/3;      % separation between matrix box and circle key

         if SETUP
%          %% DRAW UNDER-BOX OVERLAP/MARGIN SIGNIFIER AREAS
%          upperX = [vectify(repmat(1:n,[2 1]));1]*sqSz;
%          plot(boxPos(1)+upperX,boxPos(2)+shift(upperX-sqSz,-1),'k-' )
            %% LABELS
            if FIRSTORDER
                  text('Position',[boxPos+[sqSz*(n-1)/2 sqSz*(n+0.5)]],'HorizontalAlignment','center','String','Separation','FontSize',12)
                  text('Position',[boxPos+[sqSz*(n+0.5) sqSz*(n-1)/2]],'HorizontalAlignment','center','String','Overlap','Rotation',90,'FontSize',12)
                  text('Position',[boxPos+sqSz*n],'String','Radius','Rotation',-45,'FontSize',12)
            else
                  text('Position',[boxPos+[sqSz*(n-1)/2 sqSz*(n+0.5)]],'HorizontalAlignment','center','String','Separation difference','FontSize',12)
                  text('Position',[boxPos+[sqSz*(n+0.5) sqSz*(n-1)/2]],'HorizontalAlignment','center','String','Overlap difference','Rotation',90,'FontSize',12)
                  text('Position',[boxPos+sqSz*n],'String','Radius difference','Rotation',-45,'FontSize',12)
            end

            %% COLOR KEY: circles (only for second-order/difference comparisons)
            if ~FIRSTORDER
               cix = nchoosek_ix(nCats,2);
               for i = 1:nCatsc2
                  rectangle('Position',[boxPos+[(i-0.5)*sqSz csep] 0.4*[sqSz sqSz]],...
                            'FaceColor',obj.categories.colors(cix(1,i),:),'Curvature',1,...
                            'EdgeColor',obj.categories.colors(cix(1,i),:))
                  rectangle('Position',[boxPos+[csep (i-0.5)*sqSz] 0.4*[sqSz sqSz]],...
                            'FaceColor',obj.categories.colors(cix(1,i),:),'Curvature',1,...
                            'EdgeColor',obj.categories.colors(cix(1,i),:))
                  rectangle('Position',[boxPos+[(i-0.9)*sqSz csep] 0.4*[sqSz sqSz]],...
                            'FaceColor',obj.categories.colors(cix(2,i),:),'Curvature',1,...
                            'EdgeColor',obj.categories.colors(cix(2,i),:))
                  rectangle('Position',[boxPos+[csep (i-0.9)*sqSz] 0.4*[sqSz sqSz]],...
                            'FaceColor',obj.categories.colors(cix(2,i),:),'Curvature',1,...
                            'EdgeColor',obj.categories.colors(cix(2,i),:))
               end
            end
         end

         %% DRAW ACTUAL BOXES
         tt = linspace(-pi/2,pi/2,50)+pi;
         if     strcmpi(side,'left'),   lside =  1;
         elseif strcmpi(side,'right'),  lside = -1;
         end
         if     strcmpi(boxpart,'lowert'), ix = flip(ix);
         elseif strcmpi(boxpart,'diagonal') % plot as circles!
            ix = repmat(1:nchoosek(n,2),[2 1]);
            if FIRSTORDER, colors = mat2cell(obj.categories.colors,ones(n,1)); end
         end
         for i = 1:N
            if ~all(colors{i}==1)
               switch lower(side)
                  case {'left';'right'}
                     if circleNotSquare(i)
                        patchx = -0.5 + ix(1,i) + lside*cos(tt)*sizes(i)*sqScl/2;
                        patchy = -0.5 + ix(2,i) + lside*sin(tt)*sizes(i)*sqScl/2;
                     else
                        patchx = -0.5 + ix(1,i) + [0 -lside -lside 0]*sizes(i)*sqScl/2;
                        patchy = -0.5 + ix(2,i) + [1 1 -1 -1]*sizes(i)*sqScl/2;
                     end
                     patch(boxPos(1)+sqSz*patchx,boxPos(2)+sqSz*patchy,colors{i},...
                           'EdgeColor',edgecolors{i},'LineWidth',edgewidth);
                  otherwise
                     rectangle('Position',[boxPos+sqSz*(ix(:,i)'-(1+sqScl*sizes(i))/2) ...
                                           sizes(i)*sqScl*[sqSz sqSz]],'Curvature',circleNotSquare(i),...
                               'FaceColor',colors{i},'EdgeColor',edgecolors{i},'LineWidth',edgewidth)
               end
            end
         end
         axis equal ij off
      end

      function showSig(obj,varargin)
         % Parse inputs
         for v = 1:nargin-1
            if isa(varargin{v},'SetOfHyps'), sig = varargin{v}.sig;
            elseif ishandle(varargin{v}),    ax  = varargin{v};
            elseif ischar(  varargin{v}),    LGND = varargin{v};
            else                             sig = varargin{v};
            end
         end
         if ~exist('ax'  ,'var') || isempty(ax ) , ax   = gca; axis ij off; end
         if ~exist('sig' ,'var') || isempty(sig) , sig  = obj.sig;          end
         if exist('LGND','var') && ~isempty(LGND) && strcmpi(LGND,'legend'), LGND = true;
         else                                                                LGND = false;
         end

         nSigLevs = obj.nSigLevs; % 3 means [0.95 0.99 0.999] levels

         % Convert to sigmas significance, maxed to nSig(=3) and floored
         %sig = structfun(@(x) erfinv(x)*sqrt(2),sig,'UniformOutput',false);
         sigThresh = structfun(@(x) x>0.95,sig,'UniformOutput',false);
         for i = 2:nSigLevs
            for f = fieldnames(sig)'; f=f{1};
               sigThresh.(f) = sigThresh.(f) + double(sig.(f)>(1-10^-i));
            end
         end
         sigThresh = structfun(@(x) max(0,x/nSigLevs-1e-5),sigThresh,'UniformOutput',false);

         n = numel(obj.radii);
         % Time to get a-plottin'
         obj.shapesInABox(statsmat(sigThresh.ov-sigThresh.ma,sigThresh.di,sigThresh.ra),ax)
         if ~isempty(sigThresh.ra), title('Significant differences')
         else                       title('Significant values')
         end
         if LGND, obj.showSigLegend(ax); end
      end

      function showValues(obj,varargin)
         % Parse inputs
         for v = 1:nargin-1
            if isa(varargin{v},'SetOfHyps'), hi = varargin{v};
            elseif ishandle(varargin{v}),    ax = varargin{v};
            end
         end
         if ~exist('ax','var') || isempty(ax), ax = gca; axis ij off; end
         if ~exist('hi','var') || isempty(hi), hi = obj; obj = [];    end

         himat = statsmat( -hi.margins, hi.dists, hi.radii);

         if ~isempty(obj)
         lomat = statsmat(-obj.margins,obj.dists,obj.radii);
         % Normalize both to same scale
         maxval = max(abs([himat(:);lomat(:)]));

          hi.shapesInABox(himat/maxval,'left', ax)
         obj.shapesInABox(lomat/maxval,'right',ax)
         title('Values comparison')
         else
          hi.shapesInABox(himat/max(abs(himat(:))),'size',ax)
         title('Values')
         end
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


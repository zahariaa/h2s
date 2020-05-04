classdef Hypersphere < handle
   properties (SetObservable, GetObservable, AbortSet)
      centers    % [n x d] matrix: n centers of d-dimensional hyperspheres
      radii      % [1 x n] vector: radii of each of n hyperspheres
      categories(1,1) % a Categories object (see help Categories)
   end

   methods
      function obj = Hypersphere(centers,radii,varargin)
      % Constructor for Hypersphere object for hypersphere2sphere, SetOfHyps
      % e.g.:
      % hyp = Hypersphere(centers,radii,<categories>)                % (1)
      % hyp = Hypersphere('estimate',points,categories,<extraargs>)  % (2)
      % hyp = Hypersphere(hypset)                                    % (3)
      % hyp = Hypersphere(hypstruct)                                 % (4)
      % hyp = Hypersphere({hyps},radii,<extraargs>)                  % (5)
      % 
      % Constructor input options:
      %    (1) centers is an [n x d] numeric matrix, and radii is a [1 x n]
      %       vector, where n is the number of (hyper-)spheres, d is the
      %       number of dimensions. If a Categories object is provided, it
      %       is copied into the Hypersphere object. Otherwise, a dummy
      %       Categories object is automatically generated.
      %       This is the native format of Hypersphere.
      %    (2) If the first argument is 'estimate', then the following inputs
      %       (including the optional <extraargs>) are passed to
      %       estimateHypersphere. That function, as called here, outputs
      %       Hypersphere object(s). points and categories are required
      %       inputs; other optional valid inputs are the number of bootstraps
      %       and the 'stratified' or 'permute' inputs to determine whether
      %       to sample with or without replacement during bootstrapping.
      %    (3) A SetOfHyps object input is stripped down to a Hypersphere
      %       object output.
      %    (4) Mainly for compatibility with old code, this takes a struct with
      %       the same fields as a Hypersphere object and converts it to a
      %       Hypersphere object, creating a dummy Categories object if none
      %       exists.
      %    (5) centers is a cell of centers matrices, and radii is a vector
      %       of radii. This recursively calls Hypersphere on the individual
      %       centers cell elements and the corresponding radii vector
      %       elements. This makes sense for bootstrapped estimates of
      %       centers/radii.
      % 
      % N.B.: unlike SetOfHyps, volume, dists, margins, and overlap are all
      %    computed on-demand, every time (and only when) those methods are
      %    called.
      % 
      % Properties:
      %    centers    - [n x d] matrix: n centers of d-dimensional hyperspheres
      %    radii      - [1 x n] vector: radii of each of n hyperspheres
      %    categories - a Categories object (see help Categories)
      % 
      % Methods:
      %   Hypersphere.select
      %   Hypersphere.concat
      %   Hypersphere.unconcat
      %   Hypersphere.meanAndMerge
      %   Hypersphere.show
      %   Hypersphere.camera
      %   Hypersphere.movie
      %   Hypersphere.volume
      %   Hypersphere.dists
      %   Hypersphere.margins
      %   Hypersphere.overlap
      % 
      % 2018-06-07 AZ Created
      % 
      % See also SETOFHYPS, CATEGORIES, ESTIMATEHYPERSPHERE

         % create dummy object (for initializing multiple objects, e.g. option (5))
         if nargin==0; return; end

         if isa(centers,'SetOfHyps')                           % input option (3)
            obj = Hypersphere(centers.centers,centers.radii,...
                              centers.categories);
            return;
         elseif isa(centers,'Hypersphere'), obj = centers; return;
         elseif ischar(centers) && strcmpi(centers,'estimate') % input option (2)
            obj = estimateHypersphere(radii,varargin{:});
            return
         elseif isstruct(centers)                              % input option (4)
            % Helper for older struct-based code
            if isfield(centers,'categories')
               obj = Hypersphere(centers.centers,centers.radii,...
                                 centers.categories);
            elseif isa(radii,'Categories')
               obj = Hypersphere(centers.centers,centers.radii,radii);
            else
               obj = Hypersphere(centers.centers,centers.radii);
            end
            return
         elseif iscell(centers)                                % input option (5)
            % Recurse if multiple centers provided (e.g., bootstrapping)
            obj = repmat(Hypersphere(),[numel(centers) 1]);
            for i = 1:numel(centers)
               obj(i) = Hypersphere(centers{i},radii{i},varargin{:});
            end
            return
         end

         % input option (1)
         obj.centers = centers;
         obj.radii = radii;

         for v = 1:numel(varargin)
            if isa(varargin{v},'Categories')
               obj.categories = varargin{v};
            end
         end
         % Auto-generate dummy Categories if none exists by now
         if ~isa(obj.categories,'Categories')
            obj.categories = Categories(numel(obj.radii));
         end
      end

%% TODO: add in input options to choose which estimator to use!!!!!!!!!!!

      function newobj = select(self,i)
      % Hypersphere.select: outputs a Hypersphere object that has been
      %    subsampled to have one or more hyperspheres, indexed by input i
      % e.g.:
      % fewerHyps = allHyps.select(i)
      %    where i can be a logical vector or list of indices
         if isa(self,'SetOfHyps')
            newobj =   SetOfHyps(self.centers(i,:),...
                                 self.radii(i),    ...
                                 self.categories.select(i));
         else
            newobj = Hypersphere(self.centers(i,:),...
                                 self.radii(i),    ...
                                 self.categories.select(i));
         end
      end

      function self = concat(self,cats)
      % Hypersphere.concat: collapses an array of Hypersphere objects into
      %    a single Hypersphere object, with centers and radii concatenated,
      %    i.e., an [n x 1] array of Hypersphere objects with one hypersphere
      %.   each turns into one Hypersphere object that has n hyperspheres.
      % 
      % Optional input argument:
      %    cats: a categories object to embed in the concatenated object
      % 
      % e.g.:
      % oneBigHyp = aFewHypsArray.concat
      % oneBigHyp = aFewHypsArray.concat(categoriesForUnifiedHyp)
         if ~exist('cats','var'),   cats = [];   end
         if isa(self,'SetOfHyps')
            self =   SetOfHyps(cat(1,self.centers),[self.radii],cats);
         else
            self = Hypersphere(cat(1,self.centers),[self.radii],cats);
         end
      end

      function hyps = unconcat(self)
      % Hypersphere.unconcat: takes a single Hypersphere object with multiple
      %    centers/radii and breaks it out into an array of Hypersphere objects
      %    with a single hypersphere in each. Useful for JSON export
      % e.g.:
      % arrayOfIndividualHyps = hypsallinone.unconcat;
      % 
      % SEE ALSO HYPERSPHERE.MERGE, SETOFHYPS.TOJSON
         n = numel(self.radii);
         for i = 1:n
            hyps(i) = self.select(i);
         end
      end

      function self = meanAndMerge(self,FORCE)
      % Hypersphere.meanAndMerge: takes an array of Hypersphere (or SetOfHyps)
      %    objects as input. Outputs a single SetOfHyps object, where the
      %    centers and radii are the means of those in the object array. 95%
      %    confidence intervals are computed and their range values are stored
      %    in self.ci.centers and self.ci.radii. The object array originally
      %    input to this function is stored for posterity in self.ci.bootstraps,
      %    and h2s significance tests (from SetOfHyps.significance) are run,
      %    with the results embedded in self.sig and self.sigdiff.
      % e.g.:
      % hypset = arrayOfHyps.meanAndMerge
      % 
      % SEE ALSO SETOFHYPS.SIGNIFICANCE
         ci = [];
         % uses 1st object's categories
         n = numel(self);
         if n > 1 % assumes Hypersphere was bootstrapped
            if exist('FORCE','var') && FORCE
               cperms = true(n,1);
            else
               cperms = arrayfun(@(x) x.categories.ispermuted,self)';
            end
            ci.bootstraps = self(cperms);
            ci.centers = prctile(cat(3,self(cperms).centers),[2.5 97.5],3);
            ci.radii   = prctile(vertcat(self.radii),[2.5 97.5])';
            if all(cperms) % Do we always want this?
               self = Hypersphere(mean(cat(3,self(cperms).centers),3),...
                                  mean(cat(1,self(cperms).radii)),...
                                  self(1).categories);
            else % preserve best (unpermuted) estimate if exists (should be self(1))
               self = self(find(~cperms,1,'first'));
            end
         end
         self = SetOfHyps(self,ci);

         % might as well populate significance fields if we have the bootstraps
         if ~isempty(ci) && isempty(self.sig)
            self = self.significance(ci);
         end
      end

      function varargout = show(self,varargin)
      % Hypersphere.show: visualizes a Hypersphere object as a set of circles (for
      %    a 2D input) or spheres (for 3D input). If the dimensionality of the
      %    input is larger than 3, h2s is automatically run (with no arguments)
      %    and the result visualized (though the h2s-reduced SetOfHyps object is
      %    not saved unless requested).
      % 
      %    A line indicating the error of the largest single statistic of
      %    interest is plotted to the side of the visualization. Additional lines
      %    are plotted in the visualization that indicate when an overlap or
      %    margin in the original (high-dimensional) space has been visualized as
      %    the opposite (i.e., a margin or overlap, respectively). If requested,
      %    the handle(s) for these maximum error and overlap error lines are
      %    outputted.
      % 
      % IMPORTANT(!) N.B.: if a 3D (sphere) visualization figure is resized or
      %    zoomed in matlab, the error line will not resize accordingly and thus
      %    will be inaccurate.
      % 
      % e.g.:
      % hyps.show
      % hyps.show(axhandle, SETCAMERA, titleString, patchDetail)
      % hyps.show([],'')  % uses current axes (gca), with no title
      % reducedSetOfHyps = hyps.show % where hyps are dimensionality 4+; h2s ran
      % 
      % Optional inputs (argument order doesn't matter):
      %    axhandle (DEFAULT = gca): specifies the axis handle where the
      %       visualization should be rendered
      %    SETCAMERA (3D sphere visualization only, DEFAULT = true):
      %       automatically calls self.camera to standardize 3D viewpoint.
      %    titleString (DEFAULT = 'error <= [maxerror]'): specifies what to print
      %       as the plot title. The default setting prints the actual maximum
      %       error across all the individual statistics of interest.
      %    patchDetail (DEFAULT = 100): determines how much detail the circle or
      %       sphere has.
      % 
      % SEE ALSO SHOWMODEL, HYPERSPHERE.CAMERA, SETOFHYPS.PLOTOVERLAPERRORS,
      % SETOFHYPS.H2S
         maxerror = NaN;
         if isa(self,'SetOfHyps'), maxerror = max(self.error); end

         [~,XY] = showModel(self,varargin{:});

         switch size(self.centers,2)
            case 2;
               %% Draw max error bar
               if ~isnan(maxerror)
                  maxXY = max(XY);
                  maxline = plot(maxXY([1 1]) - [0 maxerror],...
                                 maxXY([2 2]),'k-','LineWidth',4);
               end
            case 3;
               %% Draw max error bar
               ax = gca;
               ax.Units = 'normalized';

               if ~isnan(maxerror)
               maxline = annotation('line',[1-maxerror/diff(ax.XLim) 1]*ax.Position(3)+ax.Position(1),...
                                           ax.Position([2 2]),'LineWidth',4,'Tag','maxline');
               end
            otherwise
               varargout = SetOfHyps(self).h2s.show(varargin);
               return
         end
         axis ij tight equal vis3d off % puts error bar at bottom right

         %% DRAW SIGN FLIP ERROR LINES
         if isa(self,'SetOfHyps'), ovlines = self.plotOverlapErrors; end

         if nargout > 0, varargout = {maxline,ovlines}; end
      end

      function camera(self,camSettings)
      % Hypersphere.camera: Sets camera, lighting, and surface properties (for a
      %    3D h2s visualization). Settings can be passed in a struct, though if
      %    visualized can come in different forms, but this function assumes
      %    no settings are explicitly passed to this, Hypersphere.camera applies
      %    standard settings.
      % 
      % e.g.:
      % hyps.camera
      % hyps.camera(camSettings)
      % 
      % Optional input:
      %    camSettings: a struct with fields corresponding to standard matlab
      %       axis camera properties, each containing the corresponding value
      %       the setting should be set to.
      %       e.g.:
      %          camSettings.CameraPosition = [1 1 1];
      %          camSettings.CameraTarget = [0 0 0];
      % 
      % SEE ALSO SHOWMODEL, HYPERSPHERE.SHOW, HYPERSPHERE.CAMERACALC
         if ~exist('camSettings','var') || isempty(camSettings)
            camSettings = self.cameraCalc;
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

      function movie(self,varargin)
      % Hypersphere.movie: Makes a dynamic h2s visualization, and (optionally)
      %    saves it to an .avi movie. Hypersphere objects that can be dynamically
      %    visualized can come in different forms, but this function assumes
      %    you've run h2s on all time points simultaneously, so that when it
      %    generates the frames, they are all in a common space.
      %    estimateHypersphere handles this automatically; so does h2s with the
      %    'joint' option.
      % 
      %    N.B.: if objects in self are high dimensional, h2s is automatically
      %       run with no arguments first, then Hypersphere.movie continues on
      %       that. Output of h2s is not saved.
      % 
      % e.g.:
      % bunchohyps = Hypersphere('estimate',X(:,:,1:end),cats); % use 3D tensor
      % bunchohyps.movie
      % bunchohyps.movie(times)
      % bunchohyps.movie(-100:10:700,'ms','dyn.avi') % saves movie to dyn.avi
      % bunchohyps.movie(times,'ms',staticFrameAxes) % plots static frames
      %                                              % in axes given
      % 
      % Optional inputs:
      %   times (DEFAULT = 1:nframes): Numeric labels for each frame, e.g., the
      %      time of that frame.
      %   timeUnits (DEFAULT = 'ms' if times given, [] otherwise): A string
      %      representing the units of the time labels above. This one string is
      %      appended to the time labels above in every frame, e.g., 1 ms, 2 ms...
      %   saveFilename (DEFAULT = []): If given, must end with '.avi'. A 20fps,
      %      uncompressed .avi movie will be saved to the path given.
      %   staticFrameAxes (DEFAULT: gca): Axes handle(s) can specify where the
      %      movie is rendered (if one handle is given), or where static frames
      %      are rendered (if multiple frames are given). In the latter case, the
      %      number of frames and axes handles should match.
      % 
      % SEE ALSO ESTIMATEHYPERSPHERE, HYPERSPHERE.CAMERA, CATEGORIES.LEGEND
         dimLow = size(self(1).centers,2);
         if dimLow > 3
            self   = self.h2s;
            dimLow = size(self(1).centers,2);
         end

         nframes   = numel(self);
         STATIC    = false;
         SAVE      = false;
         timeUnits = [];
         for v = 2:nargin
            if all(ishandle(varargin{v-1})), SAVE  = false; STATIC = true;
                                             ax    = varargin{v-1};
            elseif isnumeric(varargin{v-1}), times = varargin{v-1};
            else
               if ~isempty(strfind(varargin{v-1},'.avi'))
                  SAVE = varargin{v-1};
               else
                  timeUnits = varargin{v-1};
               end
            end
         end
         
         if ~exist('ax'  ,'var') || numel(ax)==1
            ax = repmat(gca,[1 nframes]);
         elseif numel(ax) ~= nframes
            warning('Mismatch between number of frames and axes handles provided.')
         end
         % Default frame times and label
         if exist('times','var') && ~isempty(times)
            if isempty(timeUnits), timeUnits = 'ms'; end
         else
            times = 1:nframes;
         end

         % Calculate axis bounds
         axbounds = NaN(1,dimLow*2);
         for i = 1:nframes
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
            vidObj           = VideoWriter(SAVE);
            vidObj.FrameRate = fps;
            vidObj.Quality   = 100;
            open(vidObj);
         end

         %% Do stuff for every frame
         ann = [];
         for i = 1:nframes
            axtivate(ax(i)); cla;
            % DRAW PLOT
            self(i).show(false,[],'');
            if dimLow == 3,   self.camera(camSettings);   end
            axis(axbounds)
         
            %% Generate title string
            extratxt = sprintf('\\color{black}%g %s',times(i),timeUnits);
            if STATIC
               title(extratxt,'Units','normalized','Position',[0.5 0.2 0]);
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

      %% FUNCTIONS TO COMPUTE VOLUME, DISTS, MARGINS, OVERLAPS
      % (they call the hidden static functions below)
      function V = volume(self)
      % Hypersphere.volume: Compute volume(s) of hypersphere(s)
      % e.g.:
      % V = hyp.volume
         if numel(self) > 1
            V = cell2mat_concat(arrayfun(@volume,self,'UniformOutput',false));
            return
         end
         V = self.calcVolume(self.centers,self.radii);
      end

      function D = dists(self)
      % Hypersphere.dists: Compute distance(s) between pair(s) of hyperspheres
      %    along their center-connection lines.
      % e.g.:
      % D = hyp.dists
         if numel(self) > 1
            D = cell2mat_concat(arrayfun(@dists,self,'UniformOutput',false));
            return
         end
         D = self.calcDists(self.centers,self.radii);
      end

      function M = margins(self)
      % Hypersphere.margins: Compute margin(s) between pair(s) of hyperspheres.
      %    Margins are defined as the distance along the center-connection line
      %    between two hyperspheres, minus the radii. Margins are exactly the 
      %    same as negative overlaps.
      % e.g.:
      % M = hyp.margins
         if numel(self) > 1
            M = cell2mat_concat(arrayfun(@margins,self,'UniformOutput',false));
            return
         end
         M = self.calcMargins(self.radii,self.dists);
      end

      function ov = overlap(self)
      % Hypersphere.overlap: Compute overlap(s) between pair(s) of hyperspheres.
      %    Overlaps are defined as the distance along the center-connection line
      %    between two hyperspheres, minus the radii. Overlaps are exactly the 
      %    same as negative margins.
      % e.g.:
      % OV = hyp.overlap
         ov = -self.margins;
      end

      %% SET FUNCTION TO VALIDATE CENTERS (they must be numeric, real, and finite)
      function self = set.centers(self,centers)
         arguments
            self
            centers {mustBeNumeric, mustBeFinite, mustBeReal}
         end
         % Ensure consistent formatting
         if isempty(self.radii)
            self.centers = centers;
         else
            if     size(centers,1) == numel(self.radii)
               self.centers = centers;
            elseif size(centers,2) == numel(self.radii)
               self.centers = centers';
            else
               error('size of radii and centers don''t match')
            end
         end
      end
      %% SET FUNCTION TO VALIDATE RADII (they must be numeric, nonnegative, and
      %%    finite)
      function self = set.radii(self,radii)
         arguments
            self
            radii {mustBeNumeric, mustBeFinite, mustBeNonnegative}
         end
         if ~isempty(radii) && ~isempty(self.centers) ...
            && numel(radii) ~= size(self.centers,1)
            error('size of radii and centers don''t match')
         end
         self.radii = radii(:)';
      end
      %% SET FUNCTION TO VALIDATE CATEGORIES (they must be a Categories class
      %%    object)
      function self = set.categories(self,categories)
         arguments
            self
            categories {mustBeCategoriesClass}
         end
         self.categories = categories;
      end
   end

   methods (Hidden = true)
      function camSettings = cameraCalc(self)
      % Hypersphere.cameraCalc: Calculates a "consistent" camera position (for a
      %    3D h2s visualization) by rotating everything so first 2 PC's are in
      %    view and centering the camera on the centroid.
      % 
      % e.g.:
      % hyps.cameraCalc
      % 
      % SEE ALSO HYPERSPHERE.CAMERA, HYPERSPHERE.SHOW, SHOWMODEL
         self = self.concat;
         
         % set view angle orthogonal to the PC12 plane of the locations
         eigenvecs = pca(self.centers);
         if size(eigenvecs,2)<2, normal = [0 0 0]';
         else                    normal = null(eigenvecs(:,1:2)');
         end
         
         if corr(self.centers*normal,self.radii')<0
            normal = -normal;
         end
         
         objectsCentroid = mean(self.centers,1);
         camDist = 50*sqrt(sum(std(self.centers).^2));
         camSettings.CameraPosition = objectsCentroid'-normal*camDist;
         camSettings.CameraTarget   = objectsCentroid;
      end
   end

   methods(Static, Hidden = true)
      % Made these functions static so SetOfHyps can use them
      function V = calcVolume(centers,radii)
      % Compute volumes of each hypersphere
         d = size(centers,2);
         V = (radii.^d).*(pi.^(d/2))./gamma((d/2)+1);
      end

      function D = calcDists(centers,radii)
      % Compute distances between pairs of hyperspheres (along center-connection
      %    lines)
         n = numel(radii);
         D = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               D(i) = sqrt(sum(diff(centers([a b],:)).^2));
            end
         end
      end

      function M = calcMargins(radii,dists)
      % Compute margins between pairs of hyperspheres
         n = numel(radii);
         M = zeros(1,(n^2-n)/2);
         i=0;
         for a = 1:n-1
            for b = a+1:n, i=i+1;
               M(i) = dists(i) - radii(a) - radii(b);
            end
         end
      end
   end
end

% Helper (validation) function
function mustBeCategoriesClass(c)
   if ~isa(c,'Categories')
      error('categories must be Categories class');
   end
end

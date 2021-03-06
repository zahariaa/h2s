classdef Categories
   properties
      labels                       % [1 x c] cell of strings; names of categories
      colors  = [250   138   117   % [c x 3] matrix of colors, c = # of categories
                 246   136   159
                 215   148   196
                 162   166   215
                 100   180   210
                  61   187   180
                  88   189   139
                 135   184    98
                 182   173    74
                 223   157    79]/255; 
      vectors  % [n x c] logical matrix: (n = # points, c = # categories)
   end
   properties (SetAccess = 'protected')
      ispermuted = false;
   end
   methods
      function obj = Categories(vectors,labels,colors)
      % Contructor for Categories object for Hypersphere, SetOfHyps, for
      %    use with hypersphere2sphere
      % e.g.:
      % cats = Categories(vectors)
      % cats = Categories(vectors,labels)
      % cats = Categories(vectors,labels,colors)
      % cats = Categories(vectors,[],colors)
      % 
      % The following input options are also Categories object properties.
      % vectors [required] input can be:
      %    (1) a [n x c] logical matrix, where n is the number of points and
      %       c is the number of categories. each row must have only one true
      %       value, signifiying that that point belongs to only 1 category.
      %       This is the native format for the Categories object.
      %       e.g.
      %       vectors = logical(blkdiag(true(4,1),true(5,1),true(2,1),true(9,1)));
      % 
      %    (2) a [n x 1] numeric vector, which will be converted to the format
      %       above. The number of each entry indicates the column the true
      %       value will be placed in.
      %       e.g.
      %       vectors = randi(4,20,1); % resampled (w/replacement) version of (1)
      % 
      %    (3) a cell containing the number of consecutive entries for each
      %       class
      %       e.g.
      %       vectors = {4,5,2,9};     % same as (1)
      % 
      %    (4) a structure containing vectors, labels, and colors (basically,
      %       a Categories object-like structure)
      % 
      % labels [optional]: a [1 x c] cell of strings of the category names.
      %    default = {'category 1', 'category 2', ... }
      % 
      % colors [optional]: [c x 3] numeric matrix of values in [0,1],
      %    defining the color for each category.
      %    Default: blue, red, green, orange if 4 or fewer categories given,
      %    and a circular sampling of up to 10 colors that are meant to be
      %    perceptually equidistant and at the same luminance (in L*a*b color
      %    space)
      % 
      % Methods:
      %    Categories.select
      %    Categories.internalrepmat
      %    Categories.permute
      %    Categories.slice
      %    Categories.vectorsForDistanceMatrix
      %    Categories.legend
      %    Categories.plotSamples
      %    Categories.legendText
      % 
      % 2018-06-07 AZ Created
      % 
      % See also HYPERSPHERE, SETOFHYPS

         if isstruct(vectors) && numel(vectors)==1  % input option (4)
            obj.labels  = vectors.labels;
            obj.vectors = vectors.vectors;
            if isfield(vectors,'colors')
               obj.colors  = vectors.colors;
               return
            end
         elseif islogical(vectors)                  % input option (1)
            obj.vectors = vectors;
         elseif iscell(vectors)                     % input option (3)
            % automatically generate block diagonal with each block the length
            % of each element in vectors
            cmd = sprintf('true(%u,1),',vectors{:});
            cmd = ['obj.vectors = blkdiag(' cmd(1:end-1) ');'];
            eval(cmd)
            obj.vectors = logical(obj.vectors);
         elseif isnumeric(vectors)                  % input option (2)
            if numel(vectors)==1
               obj.vectors = [];
               n = vectors;
            else
               obj.vectors = false(numel(vectors),max(vectors));
               for i = unique(vectors)'
                  obj.vectors(vectors==i,i) = true;
               end
            end
         end

         if ~isempty(obj.vectors)
            n = size(obj.vectors,2);
         end
         if ~exist('labels','var') || isempty(labels)
            obj.labels  = mat2cell([repmat('category ',[n 1]) num2str((1:n)')],ones(n,1))';
         else
            obj.labels  = labels;
         end

         if ~exist('colors','var')
            % duplicate colors if more requested than exist in c
            special = [2 5 7 10];
            if n<5;    obj.colors = obj.colors(special(1:n),:);
            else       obj.colors = obj.colors(mod((1:n)-1,10)+1,:);
            end
         else obj.colors = colors;
         end
      end

      function self = select(self,i)
      % Categories.select: outputs a Categories object that has been
      %    subsampled to have one or more categories, indexed by input i
      % e.g.:
      % fewercats = allcats.select(i)
      %    where i can be a logical vector or list of indices
         if islogical(i)
            self.labels  = self.labels(find(i));
         else
            self.labels  = self.labels(i);
         end
         self.colors  = self.colors(i,:);
         if ~isempty(self.vectors)
            self.vectors = self.vectors(:,i);
         end
      end

      function obj = internalrepmat(self,N)
      % Categories.internalrepmat: outputs a Categories object that has all
      %    its fields appended with N-1 copies of their contents.
      %    categories.vector is appropriately extended. Useful for generating
      %    multiple movie frames using the same categories object.
      % e.g.:
      % catWithStuffRepeatedNtimes = cat.internalrepmat(N)
         obj.labels  = repmat(self.labels,[1 N]);
         obj.colors  = repmat(self.colors,[N 1]);
         % Replicate vectors on block diagonal
         X = repmat({sparse(self.vectors)},[N 1]);
         obj.vectors = blkdiag(X{:});
         obj.vectors = cast(obj.vectors,class(self.vectors));
         % Put it all together
         obj = Categories(obj);
      end

      function objs = permute(self,N,STRAT_BOOTSTRAP)
      % Categories.permute: outputs a Categories object that has been
      %    subsampled to have one or more categories, indexed by input i
      % e.g.:
      % permutedcats = cats.permute
      % permutedcats = cats.permute(N)
      % permutedcats = cats.permute(N,STRAT_BOOTSTRAP)
      % 
      % N = 100 by default, is the number of bootstraps to do
      % STRAT_BOOTSTRAP = false by default, indicates whether to do a 
      %    stratified bootstrap (sampled with replacement) or a random
      %    permutation (sampled without replacement, default).
         if ~exist('N','var') || isempty(N), N=100;
         elseif N < 2,                       objs=self; return; end
         if ~exist('STRAT_BOOTSTRAP','var') || isempty(STRAT_BOOTSTRAP)
            STRAT_BOOTSTRAP = false;
         elseif numel(STRAT_BOOTSTRAP)==2
            objs = [self.permute(N,STRAT_BOOTSTRAP(1)); self.permute(N,STRAT_BOOTSTRAP(2))];
            return
         end

         [p,n] = size(self.vectors); %p=# points, n=# hyps
         vecs  = self.vectors*(1:n)';
         p     = nnz(vecs);
         if     p < 2^8 , dtype = 'uint8';
         elseif p < 2^16, dtype = 'uint16';
         elseif p < 2^32, dtype = 'uint32';
         else             dtype = 'uint64';
         end

         % Build categories objs with permuted vector identities
         self.ispermuted = true;
         if STRAT_BOOTSTRAP % resample with replacement
            self.vectors = zeros(p,n,dtype);
            objs = repmat(self,[N 1]);
            for i = 1:n
               includedvecs = find(vecs==i);
               ni = numel(includedvecs);
               for j = 1:N
                  objs(j).vectors(includedvecs,i) = includedvecs(sort(randi(ni,ni,1)));
               end
            end
         else
            self.vectors = false(p,n);
            objs = repmat(self,[N 1]);
            includedvecs = find(vecs);
            for i = 1:N
               ivec = vecs(includedvecs(randperm(p)));
               for j = 1:n
                  objs(i).vectors(ivec==j,j) = true;
               end
            end
         end
      end

      function [slicedpoints,newindices] = slice(self,points,UNIQUE)
      % Categories.slice: slices points based on self.vectors
         if islogical(self.vectors)
            ix = @(i) self.vectors(:,i);
         else
            ix = @(i) self.vectors(~~self.vectors(:,i),i);

            if exist('UNIQUE','var') && ~isempty(UNIQUE) && strcmpi(UNIQUE,'unique')
               ix = @(i) unique(ix(i));
            end
         end

         [p,n] = size(self.vectors);
         for i = 1:n
            newindices{i} = ix(i);
            slicedpoints{i} = points(newindices{i},:,:);
         end
         if nargout>1
            if islogical(newindices{1})
               newindices = cellfun(@(x) x(~~x),newindices,'UniformOutput',false);
            else
               newindices = cellfun(@(x) x-min(x)+mod(min(x),numel(x)),newindices,'UniformOutput',false);
            end
         end
      end

      function objs = leaveoneout(self)
      % Categories.leaveoneout: replicates a Categories object into p copies,
      %    but with each vectors missing one point. Useful for jackknife/leave-
      %    one-out significance testing, as is used for the margins.
      % SEE ALSO MARGINSAMPLING, SETOFHYPS.SIGNIFICANCE, HYPERSPHERE.MEANANDMERGE
         [p,n] = size(self.vectors);

         self.ispermuted = true;
         objs = repmat(self,[p 1]);
         for c = 1:p
            objs(c).vectors(c,:) = false;
         end
      end

      function ix = vectorsForDistanceMatrix(self)
      % Categories.vectorsForDistanceMatrix: converts vectors to indices to
      %    make selections from a distance matrix.
      % e.g.:
      % dists = squareform(distanceMatrix);            % extract upper triangle
      % dists(cats.select(2).vectorsForDistanceMatrix) % distances within category 2
         ix = [];
         [p,n] = size(self.vectors); %p=# points, n=# hyps

         if       isempty(self.vectors), return
         elseif ~any(self.vectors(:)>1), sels = self.vectors.*repmat((1:p)',[1 n]);
         else                            sels = self.vectors;
         end

         % find nonzero entries in all columns (DO WE CARE ABOUT INDIVIDUAL COLS?)
         uniqueEntries = unique(sels(~~sels))';
         excluded = setdiff(1:p,uniqueEntries);

         % Populate indices, take away excluded ones
         ix = true(1,nchoosek(p,2));
         ixkey = nchoosek_ix(p);
         for i = excluded
            ix(any(ixkey==i)) = false;
         end
      end

      function varargout = legend(self,pos,extratxt)
      % Categories.legend: creates a text legend, in which the category
      %    labels are rendered in their respective colors.
      % e.g.:
      % cats.legend
      % ann = cats.legend(pos)
      % [ann,anntxt] = cats.legend(pos,extratxt)
      % 
      % pos = [0.01 0.9 1 0.1] by default, placing the text at the top of
      %       the figure axis.
      % extratxt is empty by default. it's extra text to append at the end
      %       of the legend.
      % ann [optional output] is the annotation object in the figure axis
      %       containing the text
      % anntxt [optional output] is the text itself that's in the legend
         if ~exist('extratxt','var'), extratxt = []; end
         anntxt = self(1).legendText(extratxt);
         
         if ~exist('pos','var') || isempty(pos)
            pos = [0.01 0.9 1 0.1];
         end
         ann = annotation('TextBox',pos,'String',...
               anntxt,'EdgeColor','none','Interpreter','tex',...
               'HorizontalAlignment','left','Units','normalized');

         switch nargout
            case 0; varargout = {};
            case 1; varargout = {ann};
            case 2; varargout = {ann; anntxt};
         end
      end

      function varargout = plotSamples(self,points,ax)
      % Categories.plotSamples: plots inputted points according to the colors
      %    and indices in self.categories. Can provide multiple bootstraps of
      %    points: if one (or no) axis handle is provided, this plots all
      %    bootstraps in one axis. If multiple axis handles are provided, this
      %    plots the first nax individual bootstraps of the points in each of
      %    the nax axes. By default, if points are more than 3D, only the first
      %    3 dimensions are plotted.
      % e.g.:
      % ax = cats.plotSamples(points)
      % cats.plotSamples(points,ax)
      % 
      % Required input:
      %    points ([n x d x nboots] tensor): the points to be plotted. If d>3,
      %       only the 1st 3 dimensions are plotted. If nboots>1 and nax<=1, all
      %       points are plotted in the same plot. If nax>1, then the first nax
      %       matrices in the points tensor are plotted in the nax axes handles
      %       provided.
      % Optional input:
      %    ax (DEFAULT = gca): axis handle(s) for the plots. If multiple axis
      %       handles are provided, the first nax matrices in the points
      %       tensor (points(:,:,iax)) are plotted separately in each of the nax
      %       axes handles provided.
      % 
      % SEE ALSO HYPERSPHERE.PLOTSAMPLES

         if ~exist('ax','var') || isempty(ax), ax = gca; end
         nax    = numel(ax);
         nboots = size(points,3);

         for i = 1:numel(self.labels)
            cat = self.select(i);
            if cat.ispermuted
               cat.vectors = unique(find(cat.vectors));
            end

            for b = 1:nboots
               if b==1 || (b>1 && nax>1), axtivate(ax(b)); end

               if size(points,2) > 2
                  plot3(points(cat.vectors,1,b),points(cat.vectors,2,b),...
                        points(cat.vectors,3,b),'wo','MarkerFaceColor',cat.colors);
               else
                  plot( points(cat.vectors,1,b),points(cat.vectors,2,b),...
                                                'wo','MarkerFaceColor',cat.colors);
               end

               if b==nax && nax>1, break; end
            end
         end
         for a = 1:nax
            axtivate(ax(a))
            axis equal off
         end

         if nargout, varargout = {ax}; end
      end

   % end
   % methods(Access = 'private')

      function txt = legendText(self,extratxt)
      % Categories.legendText: creates a text string for use in
      %    Categories.legend. The string contains the category labels
      %    rendered in their respective colors.
      % e.g.:
      % txt = cats.legendText
      % txt = cats.legendText(extratxt)
      % 
      % extratxt = nothing by default. Adds additional text at end of legend
         if ~exist('extratxt','var') || isempty(extratxt)
            extratxt = '';
         end
         txt = '\fontsize{16}{';
         for i = 1:numel(self(1).labels)
            txt = [txt sprintf('\\color[rgb]{%1.2f %1.2f %1.2f}%s ',...
                             self(1).colors(i,:),...
                             self(1).labels{i})           ];
         end
         txt = [txt extratxt '}'];
      end
   end
end


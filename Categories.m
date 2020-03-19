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
      %    Categories.permute
      %    Categories.legend
      %    Categories.legendText
      % 
      % 2018-06-07 AZ Created
      % 
      % See also HYPERSPHERE, SETOFHYPS

         if isstruct(vectors) && numel(vectors)==1  % input option (4)
            obj.labels  = vectors.labels;
            obj.colors  = vectors.colors;
            obj.vectors = vectors.vectors;
            return
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
         end

         [p,n] = size(self.vectors); %p=# points, n=# hyps
         vecs = self.vectors*(1:n)';
         p    = nnz(vecs);

         % Build categories objs with permuted vector identities
         if STRAT_BOOTSTRAP % resample with replacement
            self.vectors = zeros(p,n);
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
            for i = 1:N
               includedvecs = find(vecs);
               ivec = vecs(includedvecs(randperm(p)));
               for j = 1:p
                  objs(i).vectors(includedvecs(j),ivec(j)) = true;
               end
            end
         end
      end

      function itis = ispermuted(self)
      % Categories.ispermuted: returns a logical, or vector of logicals,
      %    for each Categories object this is run on, indicating whether
      %    it has been permuted or not
      % e.g.:
      % if cats.ispermuted, disp('cats is permuted'); end
         if numel(self) > 1 % recurse
            itis = arrayfun(@ispermuted, self);
            return
         elseif isempty(self.vectors)
            itis = false;
            return
         end

         % actual function
         if isnumeric(self.vectors)
            itis = all(sum(self.vectors,2));
         else
            itis = numel(unique(self.vectors))==size(self.vectors,1);
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


classdef Categories
   properties
      labels
      colors
      vectors = true;
   end
   methods
      function obj = Categories(vectors,labels,colors)
         if numel(vectors)==1
            if isstruct(vectors)
               obj.labels  = vectors.labels;
               obj.colors  = vectors.colors;
               obj.vectors = vectors.vectors;
               return
            else obj.vectors = logical(eye(vectors));
            end
         elseif islogical(vectors)
            obj.vectors = vectors;
         elseif isnumeric(vectors)
            % automatically generate block diagonal with each block the length
            % of each element in vectors
            cmd = sprintf('true(%u,1),',vectors);
            cmd = ['obj.vectors = blkdiag(' cmd(1:end-1) ');'];
            eval(cmd)
         end
         n = size(obj.vectors,2);
         if ~exist('labels','var') || isempty(labels)
            obj.labels  = arrayfun(@(n) sprintf('category %u',n),1:n,'UniformOutput',false);
         else
            obj.labels  = labels;
         end
         if ~exist('colors','var')
            if n <= 4
               c = [0 0 0; 0.8 0 0; 0 0.8 0; 0 0 0.8];
            elseif verLessThan('matlab','9.4')
               c = [250   138   117
                    246   136   159
                    215   148   196
                    162   166   215
                    100   180   210
                     61   187   180
                     88   189   139
                    135   184    98
                    182   173    74
                    223   157    79]/255;
            else c = lines;
            end
            obj.colors = c(1:n,:);
         else obj.colors = colors;
         end
      end
      function obj = select(obj,i)
         obj.labels  = obj.labels(i);
         obj.colors  = obj.colors(i,:);
         obj.vectors = obj.vectors(:,i);
      end
      function objs = permute(obj,N,STRAT_BOOTSTRAP)
         if ~exist('N','var') || isempty(N), N=100;
         elseif N < 2,                       objs=obj; return; end
         if ~exist('STRAT_BOOTSTRAP','var') || isempty(STRAT_BOOTSTRAP)
            STRAT_BOOTSTRAP = false;
         end

         [p,n] = size(obj.vectors); %p=# points, n=# hyps
         vecs = obj.vectors*(1:n)';
         p    = nnz(vecs);
         % Build categories objs with permuted vector identities
         if STRAT_BOOTSTRAP % resample with replacement
            obj.vectors = zeros(p,n);
            objs = repmat(obj,[N 1]);
            for i = 1:n
               includedvecs = find(vecs==i);
               ni = numel(includedvecs);
               for j = 1:N
                  objs(j).vectors(includedvecs,i) = includedvecs(sort(randi(ni,ni,1)));
               end
            end
         else
            obj.vectors = false(p,n);
            objs = repmat(obj,[N 1]);
            for i = 1:N
               includedvecs = find(vecs);
               ivec = vecs(includedvecs(randperm(p)));
               for j = 1:p
                  objs(i).vectors(includedvecs(j),ivec(j)) = true;
               end
            end
         end
      end
%       function obj = internalrepmat(self,N)
%          obj.labels  = repmat(self.labels,[1 N]);
%          obj.colors  = repmat(self.colors,[N 1]);
%          % Replicate vectors on block diagonal
%          cmd = ['obj.vectors = blkdiag(' repmat('self.vectors,',1,N)];
%          eval([cmd(1:end-1) ');']);
%          % Put it all together
%          obj = Categories(obj);
%       end
      function txt = legendText(self,closingText)
         if ~exist('closingText','var') || isempty(closingText)
            closingText = '}';
         end
         txt = '\fontsize{16}{';
         for i = 1:numel(self(1).labels)
            txt = [txt sprintf('\\color[rgb]{%1.2f %1.2f %1.2f}%s ',...
                             self(1).colors(i,:),...
                             self(1).labels{i})           ];
         end
         txt = [txt closingText];
      end
      function [ann,anntxt] = legend(self,pos,extratxt)
         if ~exist('extratxt','var'), extratxt = []; end
         anntxt = self(1).legendText(extratxt);
         
         if ~exist('pos','var') || isempty(pos)
            pos = [0.01 0.9 1 0.1];
         end
         ann = annotation('TextBox',pos,'String',...
               anntxt,'EdgeColor','none','Interpreter','tex',...
               'HorizontalAlignment','left','Units','normalized');
      end
   end
end


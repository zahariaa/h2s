classdef Categories
   properties
      labels
      colors
      vectors = true;
   end
   methods
      function obj = Categories(vectors,labels,colors)
         if numel(vectors)==1
            obj.vectors = logical(eye(vectors));
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
            else c = colorcube(n);
            end
            obj.colors = c(1:n,:);
         end
      end
      function obj = select(obj,i)
         obj.labels  = obj.labels(i);
         obj.colors  = obj.colors(i,:);
         obj.vectors = obj.vectors(:,i);
      end
   end
end


classdef cvindex
   properties
      nCV     % number of folds
      n       % number of samples
      rp      % cell of random permutation of indices organized by fold
   end
   methods
      function obj = cvindex(n,nCV) % Constructor
      % Note: can pass bootstrappedIndices in first argument instead of n
         if numel(n) > 1 % bootstrappedIndices passed as first input
            if islogical(n)
               bootstrappedIndices = find(n);
            else
               bootstrappedIndices = n(~~n); % strip zeros
            end

            uniquebi = unique(bootstrappedIndices);
            density  = hist(bootstrappedIndices,single(uniquebi));
            n        = numel(uniquebi);
         end

         % Construct random permutations (rp) of indices and separate into nCV folds
         rp = reshape_sloppy(randperm(n),nCV);
         % Strip NaN's and convert to cell
         rp = cellfun(@(x) x(~isnan(x)), num2cell(rp,1),'UniformOutput',false);

         % Density-based expansion of indices (re-applying bootstrap)
         if exist('density','var')
            % Initialize density-based expanded indices
            obj.rp = cellfun(@(x) NaN(sum(density(x)),1),rp,'UniformOutput',false);
            % Fill in indices
            for iCV = 1:nCV
               ix = density(rp{iCV});
               j = 1;
               for i = 1:numel(ix)
                  obj.rp{iCV}(j:j+ix(i)-1) = uniquebi(rp{iCV}(i));
                  j = j+ix(i);
               end
            end
            % unit test
            if any(intersect(obj.rp{:}))
               error('same bootstrap index appears in multiple folds')
            end
         else
            obj.rp = rp;
         end

         % Populate object
         obj.nCV = nCV;
         obj.n   = n;
      end
      function ix = train(obj,i)        % Output training set for i'th fold
         a = 1:obj.nCV;
         ix = vectify(obj.rp{a(~ismembc(a,i))}); % faster than rp{setdiff(a,i)}
      end
      function ix = test(obj,i)         % Output test set for i'th fold
         ix = obj.rp{i};
      end
      function out = crossvalidate(obj,fcn,varargin)
      % Calculate fcn on training data in all folds
         X = varargin{1};
         for iCV = 1:obj.nCV
            out{iCV} = fcn(X(obj.train(iCV),:));
         end
      end
   end
end

%% DEBUG/DEMO
% n = 1000; d = 8; nCV = 5;
% X = randnball(n,d);
% cv = cvindex(n,nCV);
% radii = @(X) fminunc(@(r) max(sqrt(sum((X-repmat(r(:)',[size(X,1) 1])).^2,2))),ones(1,d));
% out = cv.crossvalidate(radii,X);


classdef cvindex
   properties
      nCV     % number of folds
      n       % number of samples
      rp      % matrix of random permutation of indices organized by fold
   end
   methods
      function obj = cvindex(n,nCV)     % Constructor
         obj.nCV = nCV;
         obj.n   = n;
         obj.rp  = reshape_sloppy(randperm(n),nCV);
      end
      function ix = train(obj,i)        % Output training set for i'th fold
         a = 1:obj.nCV;
         ix = vectify(obj.rp(:,a(~ismembc(a,i)))); % faster than setdiff(a,i))
         ix = ix(~isnan(ix));
      end
      function ix = test(obj,i)         % Output test set for i'th fold
         ix = obj.rp(:,i);
         ix = ix(~isnan(ix));
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


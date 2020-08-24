function bootstrap_test
% bootstrap_test: simplified code for testing bootstrapping and significance
%    testing in 1D case

nsims  = 10000;
nboots = 1000;
n      = 5;    % 2^n = number of samples

gt = Hypersphere([0;0],[1 1]);
[points,gt] = gt.sample(2.^[n n],nsims);

centers = NaN(2,nsims);
dists   = NaN(nsims,1);
loc_cv = repmat({[NaN;NaN]},[nsims 2]);

for isim = 1:nsims
   %hyps(isim) = Hypersphere.estimate(points,gt.categories);
   centers(1,isim) = mean(points(gt.categories.vectors(:,1),:,isim));
   centers(2,isim) = mean(points(gt.categories.vectors(:,2),:,isim));
   %dists(isim)     = Hypersphere.calcDists(centers(:,isim));

   % 2-fold cross-validated distances
   for i = 1:2
      loc_cv{isim,i} = cvCenters(points(gt.categories.vectors(:,i),:,isim),2);
   end

   %for iboot = 1:nboots
   %end
end
dists   = diff(centers);
dists2  = sign(dists).*(dists.^2);
distsCV = cvCenters2cvSqDists(loc_cv);


newfigure('dists',[3 1]);
axtivate(1); histogram(dists);        xlabel('center diffs');
axtivate(2); histogram(dists2);       xlabel('signed center distances')
axtivate(3); histogram([distsCV{:}]); xlabel('CV''d distances')
matchy('x')

keyboard
end


function loc_cv = cvCenters(points,nCVs)
   [n,d,f] = size(points);
   for icv = 1:nCVs
      cv(icv) = cvindex(n,2);
      loc_cv{icv} = [mean(points(cv(icv).train(1),:));
                     mean(points(cv(icv).train(2),:)) ];
   end
   loc_cv = mean(cell2mat(loc_cv),2);
end
   


function cvdists = cvCenters2cvSqDists(loc_cv)
% Compute SQUARED cross-validated distances from a cell array of 2-fold cross-
%    validated centers (a standard output of estimateHypersphere).
   [f,n] = size(loc_cv);
   if f > 1 % recurse
      for i = 1:f
         cvdists{i} = cvCenters2cvSqDists(loc_cv(i,:));
      end
      return
   end

   cvdists = NaN(nchoosek(n,2),1);
   i=0;
   for a = 1:n-1
      for b = a+1:n, i=i+1;
         cvdists(i) =  (loc_cv{a}(1,:)-loc_cv{b}(1,:))...
                      *(loc_cv{a}(2,:)-loc_cv{b}(2,:))';
      end
   end
   % cvdists = sign(cvdists).*sqrt(abs(cvdists));
   return
end


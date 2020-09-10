function bootstrap_test
% bootstrap_test: simplified code for testing bootstrapping and significance
%    testing in 1D case

nsims  = 1000;
nboots = 1000;
d      = 1;
n      = 5;    % 2^n = number of samples
nCVs   = 1;    % # of CV folds to average over for CV distance estimation

%% helper functions
smallertail     = @(x) min(cat(3,1-x,x),[],3);
ciprctile       = @(x,N) sum([-Inf(1,size(x,2));x;Inf(1,size(x,2))]>N)/(size(x,1)+2);
ciprctileSmTail = @(x) smallertail(ciprctile(x));


%% Set up simulations
gt = Hypersphere(zeros(2,2.^d),[1 1]);
[points,gt] = gt.sample(2.^[n n],nsims);

%% Initialize
centers   = NaN(2,2^d,nsims);
loc_cv    = repmat({NaN(2,2^d,nCVs)},[nsims 2]);
bsloc_cv  = repmat({NaN(2,2^d,nCVs)},[nsims nboots 2]);

%% Run simulations
for isim = 1:nsims
   [hyps(isim),loc_cv(isim,:)] = Hypersphere.estimate(points(:,:,isim),gt.categories);

   %% SIGNIFICANCE via bootstrapping
   bs = gt.categories.permute(nboots,false);
   %if isim==1
   parfor iboot = 1:nboots
      [bootstraps(isim,iboot),bsloc_cv(isim,iboot,:)] = ...
         Hypersphere.estimate(points(:,:,isim),bs(iboot));
   end
   %end
   stationarycounter(isim,nsims)
end
% full data measures
centers = cat(3,hyps.centers);
dists     = diff(centers);
if d==0
   dists2 = squeeze(sign(dists).*(dists.^2));
else
   dists2 = squeeze(sum(dists.^2));
end
distsCV   = cell2mat_concat(cvCenters2cvSqDists(loc_cv));
%distsCV   = sign(distsCV).*sqrt(abs(distsCV));

% bootstrapped measures
bscenters = permute(reshape(cat(3,bootstraps.centers),[2 2^d nsims nboots]),[1 2 4 3]);
bsdists   = squeeze(diff(bscenters));
if d==0
   bsdists2  = squeeze(sign(bsdists).*(bsdists.^2));
else
   bsdists2  = squeeze(sum(bsdists.^2));
end
bsdistsCV = cellify(cvCenters2cvSqDists(reshape(bsloc_cv,[nsims*nboots 2])));
bsdistsCV = reshape(bsdistsCV,[nsims nboots]);
%bsdistsCV = sign(bsdistsCV).*sqrt(abs(bsdistsCV));


%% PLOT
i=1;
newfigure('dists',[3 1]);
axtivate([1 1]); histcompare(bsdists(:,i,:),dists);
xlabel('center diffs');
axtivate([2 1]); histcompare(bsdists2(i,:),dists2);
xlabel('signed center distances')
axtivate([3 1]); histcompare(bsdistsCV(i,:),distsCV);
xlabel('cross-validated distances')

matchy('x')

%% SIGNIFICANCE TESTING
h=NaN(1,nsims);p=NaN(1,nsims);
for i=1:nsims
   %p(i) = ciprctile(bsdists2(i,:)',dists2(i));
   p(i) = ciprctile(bsdistsCV(i,:)',distsCV(i));
end
h = p<=0.05;
mean(h)

keyboard
end


%% MORE HELPER FUNCTIONS
function centers = calcCenters(points,cat)
   %centers = cell2mat_concat(cellfun(@mean,cat.slice(points),'UniformOutput',false));
   hyps = cellfun(@Hypersphere.estimate,cat.slice(points));
   centers = vertcat(hyps.centers);
end


function loc_cv = cvCenters(points,nCVs,cat)
   %pointsorig=points;
   %points = cat.slice(points,'unique');

   d       = size(points,2);
   nCats   = size(cat.vectors,2);
   loc_cv  = repmat({NaN(2,d,nCVs)},[1 nCats]);
   
   for i = 1:nCats
   for icv = 1:nCVs
      cv = cvindex(cat.select(i).vectors,2);
      %loc_cv{i}(:,:,icv) = [mean(points(cv.train(1),:));
      %                      mean(points(cv.train(2),:)) ];
      loc_cv{i}(:,:,icv) = cell2mat_concat(cv.crossvalidate(@mean,points));
      %keyboard
   end
   end
end


function cvdists = cvCenters2cvSqDists(loc_cv)
% Compute SQUARED cross-validated distances from a cell array of 2-fold cross-
%    validated centers (a standard output of estimateHypersphere).
   [f,n] = size(loc_cv);
   nCVs  = size(loc_cv{1},3);
   if f > 1 % recurse
      for i = 1:f
         cvdists{i} = cvCenters2cvSqDists(loc_cv(i,:));
      end
      return
   end

   cvdists = NaN(nchoosek(n,2),nCVs);
   i=0;
   for a = 1:n-1
      for b = a+1:n, i=i+1;
         for icv = 1:nCVs
            cvdists(i,icv) =  (loc_cv{a}(1,:,icv)-loc_cv{b}(1,:,icv))...
                             *(loc_cv{a}(2,:,icv)-loc_cv{b}(2,:,icv))';
         end
      end
   end
   cvdists = mean(cvdists,2);
   % cvdists = sign(cvdists).*sqrt(abs(cvdists));
end


function histhandles = histcompare(black,red)
   histhandles(1) = histogram(black,'Normalization','probability','FaceColor',...
                             [0 0 0],'EdgeColor','none');
   histhandles(2) = histogram( red ,'Normalization','probability','FaceColor',...
                             [1 0 0],'EdgeColor','none','BinEdges',...
                             histhandles(1).BinEdges);
end


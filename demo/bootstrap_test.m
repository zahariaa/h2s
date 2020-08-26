function bootstrap_test
% bootstrap_test: simplified code for testing bootstrapping and significance
%    testing in 1D case

nsims  = 1000;
nboots = 1000;
n      = 5;    % 2^n = number of samples
nCVs   = 1;    % # of CV folds to average over for CV distance estimation

%% helper functions
smallertail     = @(x) min(cat(3,1-x,x),[],3);
ciprctile       = @(x) sum([-Inf(1,size(x,2));x;Inf(1,size(x,2))]>0)/(size(x,1)+2);
ciprctileSmTail = @(x) smallertail(ciprctile(x));


%% Set up simulations
gt = Hypersphere([0;0],[1 1]);
[points,gt] = gt.sample(2.^[n n],nsims);

%% Initialize
centers   = NaN(2,nsims);
loc_cv    = repmat({NaN(2,1,nCVs)},[nsims 2]);
bscenters = NaN(2,nsims,nboots);
bsloc_cv  = repmat({NaN(2,1,nCVs)},[nsims nboots 2]);

%% Run simulations
for isim = 1:nsims
   %hyps(isim) = Hypersphere.estimate(points,gt.categories);
   centers(:,isim) = calcCenters(points(:,:,isim),gt.categories);
   %dists(isim)     = Hypersphere.calcDists(centers(:,isim));

   % 2-fold cross-validated distances
   for i = 1:2
      loc_cv{isim,i} = cvCenters(points(gt.categories.vectors(:,i),:,isim),nCVs);
   end

   %% SIGNIFICANCE via bootstrapping
   bs = gt.categories.permute(nboots,true);
   %if isim==1
   parfor iboot = 1:nboots
      bspoints = bs(iboot).slice(points(:,:,isim));
      bscenters(:,isim,iboot) = calcCenters(points(:,:,isim),bs(iboot));

      % 2-fold cross-validated distances
      bsloc_cv(isim,iboot,:) = cellfun(@(x) cvCenters(x,nCVs),bspoints,'UniformOutput',false);
   end
   %end
   stationarycounter(isim,nsims)
end
% full data measures
dists     = diff(centers);
dists2    = sign(dists).*(dists.^2);
distsCV   = cellify(cvCenters2cvSqDists(loc_cv));
%distsCV   = sign(distsCV).*sqrt(abs(distsCV));

% bootstrapped measures
bsdists   = diff(bscenters);
bsdists2  = sign(bsdists).*(bsdists.^2);
bsdistsCV = cellify(cvCenters2cvSqDists(reshape(bsloc_cv,[nsims*nboots 2])));
bsdistsCV = reshape(bsdistsCV,[nsims nboots]);
%bsdistsCV = sign(bsdistsCV).*sqrt(abs(bsdistsCV));


%% PLOT
i=1;
newfigure('dists',[3 1]);
axtivate([1 1]);
hdb = histogram(bsdists(:,i,:),'Normalization','probability','FaceColor',[0 0 0],'EdgeColor','none');
hda = histogram(  dists       ,'Normalization','probability','FaceColor',[1 0 0],'EdgeColor','none',...
                'BinEdges',hdb.BinEdges);
xlabel('center diffs');
axtivate([2 1]);
hd2b = histogram(bsdists2(:,i,:),'Normalization','probability','FaceColor',[0 0 0],'EdgeColor','none');
hd2a = histogram(  dists2     ,'Normalization','probability','FaceColor',[1 0 0],'EdgeColor','none',...
                 'BinEdges',hd2b.BinEdges);
xlabel('signed center distances')
axtivate([3 1]);
hcvb = histogram(bsdistsCV(i,:),'Normalization','probability','FaceColor',[0 0 0],'EdgeColor','none');
hcva = histogram(  distsCV     ,'Normalization','probability','FaceColor',[1 0 0],'EdgeColor','none',...
                 'BinEdges',hcvb.BinEdges);
xlabel('cross-validated distances')

matchy('x')

%% SIGNIFICANCE TESTING
h=NaN(1,nsims);p=NaN(1,nsims);
for i=1:nsims
   p(i) = ciprctile(squeeze(bsdists(:,i,:)));
end
h = p<=0.05;
mean(h)

keyboard
end


%% MORE HELPER FUNCTIONS
function centers = calcCenters(points,cat)
   centers = cellfun(@mean,cat.slice(points))';
end


function loc_cv = cvCenters(points,nCVs)
   [n,d,f] = size(points);
   loc_cv  = NaN(2,1,nCVs);
   for icv = 1:nCVs
      cv(icv) = cvindex(n,2);
      loc_cv(:,:,icv) = [mean(points(cv(icv).train(1),:));
                         mean(points(cv(icv).train(2),:)) ];
      %loc_cv(:,:,icv) = cellify(cvindex(n,2).crossvalidate(@mean,points));
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

   cvdists = NaN(nchoosek(n,2),1,nCVs);
   i=0;
   for a = 1:n-1
      for b = a+1:n, i=i+1;
         for icv = 1:nCVs
            cvdists(i,:,icv) =  (loc_cv{a}(1,:,icv)-loc_cv{b}(1,:,icv))...
                               *(loc_cv{a}(2,:,icv)-loc_cv{b}(2,:,icv))';
         end
      end
   end
   cvdists = mean(cvdists,3);
   % cvdists = sign(cvdists).*sqrt(abs(cvdists));
   return
end


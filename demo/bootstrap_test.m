function bootstrap_test
% bootstrap_test: simplified code for testing bootstrapping and significance
%    testing in 1D case

nsims  = 1000;
nboots = 1000;
d      = 1;
n      = 5;    % 2^n = number of samples
nCVs   = 1;    % # of CV folds to average over for CV distance estimation

%% helper functions
ciprctile       = @(x,N) sum([-Inf(1,size(x,2));x;Inf(1,size(x,2))]>N)/(size(x,1)+2);


%% Set up simulations
gt = Hypersphere(zeros(2,2.^d),[1 1]);
[points,gt] = gt.sample(2.^[n n],nsims);

%% Run simulations
for isim = 1:nsims
   hyps(isim) = Hypersphere.estimate(points(:,:,isim),gt.categories);%,nboots,'permute');

   %% SIGNIFICANCE via bootstrapping
   bs = gt.categories.permute(nboots,false);
   %if isim==1
   parfor iboot = 1:nboots
      bootstraps(isim,iboot) = Hypersphere.estimate(points(:,:,isim),bs(iboot));
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
distsCV   = [hyps.distsCV];

% bootstrapped measures
bscenters = permute(reshape(cat(3,bootstraps.centers),[2 2^d nsims nboots]),[1 2 4 3]);
bsdists   = squeeze(diff(bscenters));
if d==0
   bsdists2  = squeeze(sign(bsdists).*(bsdists.^2));
else
   bsdists2  = squeeze(sum(bsdists.^2));
end
bsdistsCV = reshape([bootstraps.distsCV],[nsims nboots]);


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
function histhandles = histcompare(black,red)
   histhandles(1) = histogram(black,'Normalization','probability','FaceColor',...
                             [0 0 0],'EdgeColor','none');
   histhandles(2) = histogram( red ,'Normalization','probability','FaceColor',...
                             [1 0 0],'EdgeColor','none','BinEdges',...
                             histhandles(1).BinEdges);
end


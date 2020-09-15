function pout = bootstrap_test(d)
% bootstrap_test: simplified code for testing bootstrapping and significance
%    testing in 1D case

if ~exist('d','var') || isempty(d)
   d = 1;
end

nsims  = 1000;
nboots = 1000;
n      = 5;    % 2^n = number of samples
nCVs   = 1;    % # of CV folds to average over for CV distance estimation

%% helper functions
ciprctile       = @(x,t) sum([-Inf(1,size(x,2));x;Inf(1,size(x,2))]>t)/(size(x,1)+2);


%% Set up simulations
gt = Hypersphere(zeros(2,2.^d),[1 1]);
[points,gt] = gt.sample(2.^[n n],nsims);
bhyps = repmat(Hypersphere(),[nsims nboots]);
pout = NaN(nsims,nchoosek(2,2));

%% Run simulations
for isim = 1:nsims
   hyptmp = Hypersphere.estimate(points(:,:,isim),gt.categories);%,nboots,'calcStats');
   hyps(isim) = hyptmp(1);
   
   aa = SetOfHyps(hyptmp(1)).significance(points(:,:,isim),nboots);
   bhyps(isim,:) = aa.ci.permutations;
   pout(isim,:) = aa.sig.di;
   stationarycounter(isim,nsims)
end

% full data measures
centers   = cat(3,hyps.centers);
dists     = diff(centers);
if d==0
   dists2 = squeeze(sign(dists).*(dists.^2));
else
   dists2 = squeeze(sum(dists.^2));
end
distsCV   = [hyps.distsCV];

% bootstrapped measures
bscenters = permute(reshape(cat(3,bhyps.centers),[2 2^d nsims nboots]),[1 2 4 3]);
bsdists   = squeeze(diff(bscenters));
if d==0
   bsdists2  = squeeze(sign(bsdists).*(bsdists.^2));
else
   bsdists2  = squeeze(sum(bsdists.^2));
end
bsdistsCV = reshape([bhyps.distsCV],[nsims nboots]);


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
mean(pout<=0.05)

%keyboard
end


%% MORE HELPER FUNCTIONS
function histhandles = histcompare(black,red)
   histhandles(1) = histogram(black,'Normalization','probability','FaceColor',...
                             [0 0 0],'EdgeColor','none');
   histhandles(2) = histogram( red ,'Normalization','probability','FaceColor',...
                             [1 0 0],'EdgeColor','none','BinEdges',...
                             histhandles(1).BinEdges);
end


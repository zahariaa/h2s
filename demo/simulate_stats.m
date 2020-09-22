function [sigtest,estimates,hyps,groundtruth] = simulate_stats(drn,s,estimator,nsims,nboots)
% simulate_stats: simulate null distributions from different h2s statistical
%    tests with a few special cases

% # of simulations
if ~exist('nsims' ,'var') || isempty(nsims ), nsims  = 1000; end
% # of bootstraps during statistical testing
if ~exist('nboots','var') || isempty(nboots), nboots = 1000; end
if ~exist('s'     ,'var') || isempty(   s  ),      s = 1;    end
if ~exist('estimator','var') || isempty(estimator),estimator=[]; end
% s is the scenario chosen from the following:
% (1) 2 balls with 0 overlap.
% (2) 2 balls with 0 distance.
% (3) 2 balls with same (non-zero) radius.
% (4) 3 balls, 2 with same (non-zero) overlaps
% (5) 3 balls, 2 with same (non-zero) distances
% Can combine scenarios:
% (a) is (1) and (2) in 3 balls. Balls 1 and 2 have 0 distance, ball 3 placed
%    to have 0 overlap with ball 1.
% ... though this gets annoying if playing with r- and n-ratios.

%% TODO: ALL TESTS WITH DIFFERENT ESTIMATORS?
%% Testing as a function of d, n, n-ratio, and r-ratio (2^these):
sigthresh = 0.05;
ds =  1:7; % # of dimensions
rs = -1:3; % radius ratios
ns =  6:8; % # of samples to test 
nd = numel(ds);
nr = numel(rs);
nn = numel(ns);

STANDALONE = exist('drn','var') && ~isempty(drn);
if     ~STANDALONE % do nothing
elseif isnumeric(drn) && numel(drn)==1
   [n,r,d] = ind2sub([nn nr nd],drn);
   ds = ds(d);
   ns = ns(n);
   rs = rs(r);
else
   [ds,rs,ns] = dealvec(drn);
end
nd = numel(ds);
nr = numel(rs);
nn = numel(ns);

sigfield = {'sig',     'sigdiff'};
measures = {'overlap', 'distance', 'radius difference', ...
            'overlap difference', 'distance difference'};
mnames   = {'overlap', 'distsCV', 'radii', ...
            'overlap', 'dists'};
testtype = {'bootstrap' ,'permute'     ,'bootstrap' ,'bootstrap' ,'bootstrap' };
testrslt = {'bootstraps','permutations','bootstraps','bootstraps','bootstraps'};
% % debug
% for s = 1:5
%    testScenario(s,2.^ds,2.^rs,measures,mnames);
% end
% keyboard

nm  = numel(measures);
nc2 = 1+2*double(s>2);

basedir = '';%'/moto/nklab/users/az2522/';
simfolder = '20200922/';
simfolder = [basedir 'data/statsim/' simfolder];
if ~exist(simfolder,'dir'), mkdir(simfolder); end
simfile = @(s,d,r,n) sprintf('%s%ss%g_d%g_r%g_n%g.mat',simfolder,estimator,s,d,r,n);
% Initialize data for saving (could probably change from cell to tensor)
if ~STANDALONE % Collect all simulations, make figures
   sigtest   = repmat({false(nn,nsims       )},[nm nd nr]);
   bootsamps = repmat({  NaN(nn,nsims,nboots)},[nm nd nr]);
   hyps      = repmat({repmat(Hypersphere(),[nn nsims])},[nm nd nr]);
end
% else do only one simulation, don't make figures

DISPLAYED = false;
FINISHED  = false;
[di,ri,ni]= deal(0);
for d = 1:nd
   for r = 1:nr
      % Generate points for the different scenarios
      gt = generateScenario(s,2^ds(d),2^rs(r));

      for n = 1:nn
         if exist(simfile(s,ds(d),rs(r),ns(n)),'file')
            [hyp,gt,sigtmp,bootmp] = hyp_load(simfile,s,ds(d),rs(r),ns(n),DISPLAYED);
            %% TODO: check nboots/nsims
         else
            if ~DISPLAYED
               fprintf('Simulating points and estimating hyperspheres...      ')
            end
            [points,gt] = gt.sample(2.^[ns(n) ns(n)*ones(1,1+double(s>2))],nsims);
            % Simulate points
            hyp = Hypersphere.estimate(points,gt.categories,'independent',estimator);%.meanAndMerge(true);
            % save in-progress fits
            savtmp = struct('hyp',hyp,'gt',gt,'n',ns(n),'d',ds(d),'r',rs(r));
            save(simfile(s,ds(d),rs(r),ns(n)),'-struct','savtmp')
            stationarycounter([d r n],[nd nr nn])
         end

         if ~exist('sigtmp','var') || isempty(sigtmp)
            if ~DISPLAYED
               fprintf('\nSignificance testing...      ')
               DISPLAYED = true;
            end
            
            if ~exist('points','var') % regenerate points
               gt.resetRandStream
               [points,gt] = gt.sample(2.^[ns(n) ns(n)*ones(1,1+double(s>2))],nsims);
            end
            % Assess significance on samples
            sigtmp = NaN(nsims,nc2);
            bootmp = NaN(nsims,nc2,nboots);
            for b = 1:nsims
               hyptmp = SetOfHyps(hyp(b)).significance(points(:,:,b),nboots,testtype{s},estimator);
               sigtmp(b,:) = hyptmp.(sigfield{1+double(s>2)}).(mnames{s}(1:2));
               bootmp(b,:,:) = cat(1,hyptmp.ci.(testrslt{s}).(mnames{s}))';
               stationarycounter(b,nsims)
            end
            sigtmp = sigtmp(:,1+2*double(s==3));      % pick 3rd comparison only when s==3
            switch s
               case {1,2}; bootmp = bootmp(:,1,:);
               case 3;     bootmp = diff(bootmp(:,2:3,:),[],2);   % difference of  last two measures
               case {4,5}; bootmp = diff(bootmp(:,1:2,:),[],2);   % difference of first two measures
            end
            % save in-progress fits
            savtmp = struct('hyp',hyp,'gt',gt,'sigtmp',sigtmp,'bootmp',bootmp,'n',ns(n),'d',ds(d),'r',rs(r));
            save(simfile(s,ds(d),rs(r),ns(n)),'-struct','savtmp')
            stationarycounter([d r n],[nd nr nn])
            FINISHED = true;
         end

         if STANDALONE && FINISHED
            fprintf('Simulations complete.\n');
            return
         elseif ~STANDALONE % collect data
            hyps{s,d,r}(n,:)        = hyp;
            groundtruth{s,d,r}      = gt;
            sigtest{s,d,r}(n,:)     = sigtmp;
            bootsamps{s,d,r}(n,:,:) = bootmp;
         end
         clear sigtmp
%       keyboard
      end
   end
end
if STANDALONE, fprintf('Simulations complete.\n'); return; end

%% Extract data: collect estimates (e.g., overlaps, distances)
estimates = repmat({NaN(nn,nsims  )},[nm nd nr]);
bootprc   = repmat({NaN(nn,nsims,2)},[nm nd nr]);
for d = 1:nd
   for r = 1:nr
      for n = 1:nn
         switch s
            case {1,2}; estimates{s,d,r}(n,:) =             cat(1,hyps{s,d,r}(n,:).(mnames{s}));
            case   3  ; estimates{s,d,r}(n,:) = diff(indexm(cat(1,hyps{s,d,r}(n,:).(mnames{s})),[],2:3),[],2);
            case {4,5}; estimates{s,d,r}(n,:) = diff(indexm(cat(1,hyps{s,d,r}(n,:).(mnames{s})),[],1:2),[],2);
         end
         if s==1
         bootprc{s,d,r}(n,:,:) = prctile([-Inf(nsims,1) squeeze(bootsamps{s,d,r}(n,:,:))]',...
                                         [sigthresh*100 100])';
         else
         bootprc{s,d,r}(n,:,:) = prctile([-Inf(nsims,1) squeeze(bootsamps{s,d,r}(n,:,:)) Inf(nsims,1)]',...
                                         [0 100] + [1 -1]*sigthresh*100/2)';
         end
      end
   end
end
switch s
   case {1,4}; normfactor = @(r) min(2^rs(r).^[1 -1]);
   case {2,5}; normfactor = @(r) min(2^rs(r).^[0 -1]);
   case   3  ; normfactor = @(r)     2^-rs(r);
end

%% Are the signifiance tests significantly performing correctly?
ptests = NaN(nd,nn,nr,nm);
for d = 1:nd
   for r = 1:nr
      for n = 1:nn
         sigcount = [nsims 0] + [-1 1]*sum(sigtest{s,d,r}(n,:));
         ptests(d,n,r,s) = testCountDistribution_chi2(sigcount,[1-sigthresh sigthresh]);
      end
   end
end


%% Plot analyses
[d,r,n] = deal(1);
nShown = [1 nn  1 nn];
dShown = [1  1 nd nd];
rShown = [1 nr  1 nr];

fh = newfigure([3 4],sprintf('%u%s',s,measures{s}));
ax = fh.a.h([9 10]);

for i = 1:3
   sampsz = (2^(i+4))*ones(1,2+double(s>2));
   generateScenario(s,2,2^rs(rShown(i))).plotSamples(sampsz,fh.a.h(i));
   if s<3, title(sprintf('Samples: n_1= %u, n_2= %u',         sampsz))
   else    title(sprintf('Samples: n_1= %u, n_2= %u, n_3= %u',sampsz))
   end
end

axtivate(4)
plot(ns,100*squeeze(ptests(:,:,r,s)),'-o')
plot(ns([1 end]),100*sigthresh*[1 1],'--k')
if numel(ns)>1
   xlim(ns([1 end])); logAxis(2)
end
ylim([0 max(ylim)])
xlabel('samples')
ylabel('Significance of false positives compared to 5%')
title('Test performance')

axtivate(5)
plot([0 nsims],[0 0],'k--')
for i = 1:nsims
   if sigtest{s,d,r}(n,i), linecol = [0   0   0  ];
   else                    linecol = [0.5 0.5 0.5];
   end
   plot([i i],squeeze(bootprc{s,d,r}(n,i,:)),'-','Color',linecol,'LineWidth',4)
   plot(i,estimates{s,d,r}(n,i),'wo','MarkerFaceColor','r','LineWidth',1.5)
end
ylabel(sprintf('%s\nestimates (red)\nconfidence intervals (black = significant)',measures{s}))
xlabel('Simulation #')
title({'Estimate (red)' 'boostrapped CI (black/grey)'})

axtivate(6)
hb = histogram(bootsamps{s,d,r}(n,:,:),100,'Normalization','pdf',...
               'Orientation','horizontal','FaceColor',[0 0 0],'EdgeColor','none');
he = histogram(estimates{s,d,r}(n,:),'BinEdges',hb.BinEdges,'Normalization','pdf',...
               'Orientation','horizontal','FaceColor',[1 0 0],'EdgeColor','none');
ylabel(measures{s})
xlabel('pdf')
title({'Boostrapped estimates from' 'all simulations (black) and' 'estimates (red)'})
matchy(fh.a.h(5:6),'y')

%% TODO: different colors for different conditions

for i = 3:4
   axtivate(ax(i-2))
   % histogram of all estimates (newly calculated for new dimensionality)
   hb = histogram(bootsamps{s,dShown(i),r}(nShown(i),1,:),...
               'Normalization','probability','FaceColor',[0 0 0],'EdgeColor','none');
   he = histogram(estimates{s,dShown(i),r}(nShown(i),:),'BinEdges',hb.BinEdges,...
               'Normalization','probability','FaceColor',[1 0 0],'EdgeColor','none');
   plot([1 1]*estimates{s,dShown(i),r}(nShown(i),1),[0 0.25],'w-' ,'LineWidth',3)
   plot([1 1]*estimates{s,dShown(i),r}(nShown(i),1),[0 0.25],'r--','LineWidth',2)
   xlabel(measures{s})
   ylabel('probability')
   title(sprintf('Estimates (n = %u, d = %u)',2^ns(nShown(i)),2^ds(dShown(i))))
end
matchy(ax)

for d = nd:-1:1
   plotErrorPatch(rs,squeeze(indexm(cat(3,estimates{s,d,:}),n)) ...
                     .*(ones(nsims,1)*arrayfun(normfactor,1:nr)),...
                  fh.a.h(7),[1 d*[1 1]/(nd+1)])%,'sem')
end
xlim(rs([1 end])); xticks(rs); xticklabels(2.^rs);
xlabel('radius ratio')
ylabel(measures{s})
title('Estimate vs r-ratio')

for r = nr:-1:1
   plotErrorPatch(ds,squeeze(indexm(cat(3,estimates{s,:,r}),n))*normfactor(r),...
                  fh.a.h(8),[r*[1 1]/(nr+1) 1])%,'sem')
end
xlim(ds([1 end])); xticks(ds); xticklabels(2.^ds);
xlabel('dimensions')
ylabel(measures{s})
title('Estimate vs d')

axtivate(11)
for d = 1:nd
   a = permute(indexm(cat(3,sigtest{s,d,:}),n),[2 3 1]);
   plot(rs,100*mean(a),'-','LineWidth',2,'Color',[1 [d d]/(nd+1)])
end
plot(rs([1 end]),100*sigthresh*[1 1],'--k')
xlim(rs([1 end])); xticks(rs); xticklabels(2.^rs);
xlabel('radius ratio')
ylabel('False positive rate (%)')
title('False positive rate vs r-ratio')

axtivate(12)
for r = 1:nr
   a = permute(indexm(cat(3,sigtest{s,:,r}),n),[2 3 1]);
   plot(ds,100*mean(a),'-','LineWidth',2,'Color',[[r r]/(nr+1) 1])
end
plot(ds([1 end]),100*sigthresh*[1 1],'--k')
xlim(ds([1 end])); xticks(ds); xticklabels(2.^ds);
xlabel('dimensions')
ylabel('False positive rate (%)')
title('False positive rate vs d')
printFig;

% keyboard
return
end



%% helper functions
function hyp = generateScenario(s,d,r)
% generateScenario: Creates a hypersphere object to serve as 'ground truth' in
%    the simulations above.
% Inputs:
%    s: (scalar) the scenario number to choose
%    d: (scalar) the dimensionality of the hyperspheres
%    r: (scalar) the radius of one hypersphere (the other(s) is/are 1)
   switch s
      case 1 % OVERLAP. Null distribution: 0 overlap.
         hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1)],[r 1]);
      case 2 % INTER-CENTER DISTANCE. Null distribution: 0 distance.
         hyp = Hypersphere( zeros(2,d),[r 1]);
      case {3,4,5} % (Need 3+ balls.)
         % RADIUS DIFFERENCE. Null distribution: radii are the same.
         % OVERLAP DIFFERENCE. Null distribution: overlaps are the same.
         % INTER-CENTER DISTANCE DIFFERENCE. Null distribution: distances are the same.
         hyp = Hypersphere([zeros(1,d);r+0.75 zeros(1,d-1);0 r+0.75 zeros(1,d-2)],[r 1 1]);
         % hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1);ones(1,d)],[r 1 1]);
      % case 4 % OVERLAP DIFFERENCE. Null distribution: overlaps are the same. (Need 3+ balls.)
         % hyp = Hypersphere([zeros(1,d);3*r/2 zeros(1,d-1);0 3*r/2 zeros(1,d-2)],[r 1 1]);
   end
end

% debug function
function fh = testScenario(s,ds,rs,measures,mnames)
   mnames([2 5]) = {'dists'}; % because we are running on ground truth hyps
   nd = numel(ds);
   nr = numel(rs);

   fh = newfigure([nd nr],sprintf('%u%s',s,measures{s}));
   fh = newfigure([nd nr],sprintf('%u%s_values',s,measures{s}),fh);
   for d = 1:nd
      for r = 1:nr
         stationarycounter([d r],[nd nr])
         hyp = generateScenario(s,ds(d),rs(r));
         axtivate([d r],fh.f(1)); hyp.snip.show
         axtivate([d r],fh.f(2)); SetOfHyps(hyp).showValues
         % debug
         valcheck = hyp.(mnames{s});
         if     numel(valcheck)==1 && valcheck==0                                 % do nothing
         elseif numel(valcheck) >1 && any(histcounts(valcheck,numel(valcheck))>1) % do nothing
         else   keyboard
         end
      end
      drawnow
   end
   return
end

function [hyp,gt,sigtmp,bootmp] = hyp_load(simfile,s,d,r,n,DISPLAYED)
   checkFile(simfile(s,d,r,n),d,r,n)
   sigtmp = []; bootmp = [];

   load(simfile(s,d,r,n),'hyp','gt','sigtmp','bootmp');
   if ~DISPLAYED
      fprintf('Loaded simulations.\n')
   end
end

function checkFile(simfile,d,r,n)
   if exist(simfile,'file')
      A = load(simfile,'d','r','n');
      for fn = fieldnames(A)';fn=fn{1};
         if eval([fn ' ~= A.' fn])
            error('parameter mismatch') % TODO: reconcile instead
         end
      end
   end
end


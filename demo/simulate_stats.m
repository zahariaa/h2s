function [sigtest,estimates,hyps,groundtruth] = simulate_stats(s,nsims,nboots,d,r,n)
% simulate_stats: simulate null distributions from different h2s statistical
%    tests with a few special cases

% # of simulations
if ~exist('nsims' ,'var') || isempty(nsims ), nsims  = 100; end
% # of bootstraps during statistical testing
if ~exist('nboots','var') || isempty(nboots), nboots = 100; end
if ~exist('s'     ,'var') || isempty(   s  ),      s = 1;   end
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
%% Testing as a function of d, n, n-ratio, and r-ratio:
ns = 2.^( 6:8); % # samples to test
ds = 2.^( 1:5); % # dimensions
rs = 2.^(-1:3); % # radii

sigthresh = 0.05;

sigfield = {'sig',     'sigdiff'};
measures = {'overlap', 'distance', 'radius difference', ...
            'overlap difference', 'distance difference'};
mnames   = {'overlap', 'distances2CrossValidated', 'radii', ...
            'overlap', 'distances2CrossValidated'};
% % debug
% for s = 1:5
%    testScenario(s,ds,rs,measures,mnames);
% end
% keyboard

nn = numel(ns);
nd = numel(ds);
nr = numel(rs);
nm = numel(measures);

% Initialize data for saving (could probably change from cell to tensor)
nc2       = 1+2*double(s>2);
sigtest   = repmat({NaN(nn,nsims       )},[nm nd nr]);
bootsamps = repmat({NaN(nn,nsims,nboots)},[nm nd nr]);
hyps      = repmat({repmat(Hypersphere(),[nn nsims])},[nm nd nr]);

if ~exist('d','var')
   % maybe input should be ds, etc?
   % Collect all simulations, make figures
   [d,r,n]    = deal(1);
   STANDALONE = false;
else
   %% TODO: do only one simulation, don't make figures
   % ds = ds(d); rs = rs(r); ns = ns(n);
   % [nd,nr,nn] = deal(1);
   STANDALONE = true;
end

simfile = @(s,d,r,n) sprintf('statsim_s%g_d%g_r%g_n%g.mat',...
                             s,log2(ds(d)),log2(rs(r)),log2(ns(n)));
% check if existing ns,ds,rs match saved ones
checkFile(simfile(s,d,r,n),ns,ds,rs)
DISPLAYED = false;
FINISHED  = false;
for d = 1:nd
   for r = 1:nr
      % Generate points for the different scenarios
      gt = generateScenario(s,ds(d),rs(r));

      for n = 1:nn
         if exist(simfile(s,d,r,n),'file')
            [hyp,gt,sigtmp,bootmp] = hyp_load(simfile(s,d,r,n),ns,ds,rs,DISPLAYED);
         else
            if ~DISPLAYED
               fprintf('Simulating points and estimating hyperspheres...      ')
            end
            [points,gt] = gt.sample([ns(n) ns(n)*ones(1,1+double(s>2))],nsims);
            % Simulate points
            hyp = Hypersphere.estimate(points,gt.categories,'independent');%.meanAndMerge(true);
            % save in-progress fits
            save(simfile(s,d,r,n),'hyp','gt','ns','ds','rs')
            stationarycounter([d r n],[nd nr nn])
         end

         if ~exist('sigtmp','var') || isempty(sigtmp)
            if ~DISPLAYED
               fprintf('\nSignificance testing...      ')
               DISPLAYED = true;
            end
            
            if ~exist('points','var') % regenerate points
               gt.resetRandStream
               [points,gt] = gt.sample([ns(n) ns(n)*ones(1,1+double(s>2))],nsims);
            end
            % Assess significance on samples
            sigtmp = NaN(nsims,nc2);
            bootmp = NaN(nboots,nc2,nsims);
            parfor b = 1:nsims
               hyptmp = SetOfHyps(hyp(b)).significance(points(:,:,b),nboots);
               sigtmp(b,:) = hyptmp.(sigfield{1+double(s>2)}).(mnames{s}(1:2));
               bootmp(:,:,b) = cat(1,hyptmp.ci.bootstraps.(mnames{s}));
            end
            sigtmp = sigtmp(:,1+2*double(s==3));      % pick 3rd comparison only when s==3
            switch s
               case {1,2}; bootmp = bootmp(:,1,:);
               case 3;     bootmp = diff(bootmp(:,2:3,:),[],2);   % difference of  last two measures
               case {4,5}; bootmp = diff(bootmp(:,1:2,:),[],2);   % difference of first two measures
            end
            % save in-progress fits
            save(simfile(s,d,r,n),'hyp','gt','sigtmp','bootmp','ns','ds','rs')
            stationarycounter([d r n],[nd nr nn])
            FINISHED = true;
         end

         if STANDALONE && FINISHED
            return
         elseif ~STANDALONE % collect data
            hyps{s,d,r}(n,:)        = hyp;
            groundtruth{s,d,r}      = gt;
            sigtest{s,d,r}(n,:)     = sigtmp;
            bootsamps{s,d,r}(n,:,:) = bootmp;
         end
         clear sigtmp
      end
   end
end

if STANDALONE, return; end

%% Extract data: collect estimates (e.g., overlaps, distances)
estimates = repmat({NaN(nn,nsims  )},[nm nd nr]);
bootprc   = repmat({NaN(nn,nsims,2)},[nm nd nr]);
for d = 1:nd
   for r = 1:nr
      for n = 1:nn
         switch s
            case {1,2}; estimates{s,d,r}(n,:) =                   hyps{s,d,r}(n,:).(mnames{s});
            case 3;     estimates{s,d,r}(n,:) = diff(indexm(cat(1,hyps{s,d,r}(n,:).(mnames{s})),[],2:3),[],2);
            case {4,5}; estimates{s,d,r}(n,:) = diff(indexm(cat(1,hyps{s,d,r}(n,:).(mnames{s})),[],1:2),[],2);
         end
         bootprc{s,d,r}(n,:,:) = prctile([-Inf(1,nsims);squeeze([bootsamps{s,d,r}(n,:,:)]);Inf(1,nsims)],...
                                         [0 100] + [1 -1]*sigthresh*100/2)';
      end
   end
end

%% Are the signifiance tests significantly performing correctly?
ptests = NaN(nd,nn,nr,nm);
for d = 1:nd
   for r = 1:nr
      for n = 1:nn
         sigcount = [nsims 0] + [-1 1]*sum(1-sigtest{s,d,r}(n,:)<=sigthresh/2);
         ptests(d,n,r,s) = testCountDistribution_chi2(sigcount,[1-sigthresh/2 sigthresh/2]);
      end
   end
end


%% Plot analyses
nShown = [1 nn  1 nn];
dShown = [1  1 nd nd];
rShown = 1;

fh = newfigure([3 4],sprintf('%u%s',s,measures{s}));
ax = fh.a.h([9 10]);

for i = 1:3
   sampsz = (2^(i+4))*ones(1,2+double(s>2));
   generateScenario(s,2,rs(rShown)).plotSamples(sampsz,fh.a.h(i));
   if s<3, title(sprintf('Samples: n_1= %u, n_2= %u',         sampsz))
   else    title(sprintf('Samples: n_1= %u, n_2= %u, n_3= %u',sampsz))
   end
end

axtivate(4)
plot(log2(ns),100*squeeze(ptests(:,:,r,s)))
plot(log2(ns([1 end])),100*sigthresh*[1 1],'--k')
logAxis(2)
xlabel('samples')
ylabel('Significance of false positives compared to 5%')
title('Test performance')

axtivate(5)
plot([0 nsims],[0 0],'k--')
d=1;n=1;r=1;
for i = 1:nsims
   if bootprc{s,d,r}(n,i,1)<0 && bootprc{s,d,r}(n,i,2)>0, linecol = [0.5 0.5 0.5];
   else                                                   linecol = [0   0   0  ];
   end
   plot([i i],squeeze(bootprc{s,d,r}(n,i,:)),'-','Color',linecol,'LineWidth',4)
   plot(i,estimates{s,d,r}(n,i),'wo','MarkerFaceColor','r','LineWidth',1.5)
end
ylabel(sprintf('%s\nestimates (red)\nconfidence intervals (black = significant)',measures{s}))
xlabel('Simulation #')
title({'Estimate (red)' 'boostrapped CI (black/grey)'})

axtivate(6)
[counts,bins] = hist(estimates{s,d,r}(n,:));
barh(bins,counts/nsims,'FaceColor',[1 0 0],'EdgeColor','none')%[1 1 1])
[counts,bins] = hist(bootsamps{s,d,r}(n,:),100);
barh(bins,counts/(nsims*nboots),'FaceColor',[0 0 0],'EdgeColor','none')%[1 1 1])
ylabel(measures{s})
xlabel('frequency')
title({'Boostrapped estimates from' 'all simulations (black) and' 'true estimates (red)'})
matchy(fh.a.h(5:6),'y')

%% TODO: different colors for different conditions

for i = 3:4
   axtivate(ax(i-2))
   [counts,bins] = hist(estimates{s,dShown(i),r}(nShown(i),:));
   bar(bins,counts/nsims,'FaceColor',[0 0 0],'EdgeColor',[1 1 1])
   plot([0 0],[0 0.25],'r-','LineWidth',2)
   xlabel(measures{s})
   ylabel('frequency')
   title(sprintf('Estimates (n = %u, d = %u)',ns(nShown(i)),ds(dShown(i))))
end
matchy(ax)

for i = 2:3
   plotErrorPatch(log2(ns),estimates{s,dShown(i)}',fh.a.h(7),[0 0 0])%,'sem')
end
logAxis(2)
xlabel('samples')
ylabel(measures{s})
title('Estimate vs n')

for i = 2:3
   plotErrorPatch(log2(ds),squeeze(indexm(cat(3,estimates{s,:}),nShown(i))),fh.a.h(8),[0 0 0])%,'sem')
end
logAxis(2)
xlabel('dimensions')
ylabel(measures{s})
title('Estimate vs d')

axtivate(11)
for i = 2:3
   plot(log2(ns),100*mean(sigtest{s,dShown(i)}<=sigthresh/2,2),'-k','LineWidth',2)
end
plot(log2(ns([1 end])),100*sigthresh*[1 1],'--k')
logAxis(2)
xlabel('samples')
ylabel('False positive rate (%)')
title('False positive rate vs n')

axtivate(12)
for i = 2:3
   plot(log2(ds),100*squeeze(mean(indexm(cat(3,sigtest{s,:}),nShown(i))<=sigthresh/2,2)),'-k','LineWidth',2)
end
plot(log2(ds([1 end])),100*sigthresh*[1 1],'--k')
logAxis(2)
xlabel('dimensions')
ylabel('False positive rate (%)')
title('False positive rate vs d')
printFig;

keyboard
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

function [hyp,gt,sigtmp,bootmp] = hyp_load(simfile,ns,ds,rs,DISPLAYED)
   checkFile(simfile,ns,ds,rs)
   sigtmp = []; bootmp = [];

   load(simfile,'hyp','gt','sigtmp','bootmp');
   if ~DISPLAYED
      fprintf('Loaded simulations.\n')
   end
end

function checkFile(simfile,ns,ds,rs)
   if exist(simfile,'file')
      A = load(simfile,'ns','ds','rs');
      for fn = fieldnames(A)';fn=fn{1};
         if eval([fn '(1) ~= A.' fn '(1)'])
            error('parameter mismatch') % TODO: reconcile instead
         end
      end
   end
end
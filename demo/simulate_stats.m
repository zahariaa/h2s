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
if s==3, rs = 0; end % DISABLE RADIUS RATIOS FOR RADIUS DIFFERENCE TEST
ns =  4:10; % # of samples to test
nd = numel(ds);
nr = numel(rs);
nn = numel(ns);

if strcmpi(estimator,'mcmc')
   testtype = repmat({'bootstrap' },[1 7]);
   testrslt = repmat({'bootstraps'},[1 7]);
else
   testtype = [{'bootstrap'  'permute'     } repmat({'bootstrap' },[1 4]) {'bootstrap' }];
   testrslt = [{'bootstraps' 'permutations'} repmat({'bootstraps'},[1 4]) {'jackknives'}];
end
measures = {'overlap', 'distance', 'radius difference', ...
            'overlap difference', 'distance difference','overlap difference',...
            'margin'};
mnames   = {'overlap', 'dists', 'radii', ...
            'overlap', 'dists','overlap','margins'};
if isempty(strfind(measures{s},'difference')), sigOrDiff = 'sig';
else                                           sigOrDiff = 'sigdiff';
end
if strcmp(mnames{s},'dists'), cvdists = 'cvdists';
else                          cvdists = 'nocvdists';
end

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
nm = numel(measures);

% % debug
% for s = 1:6
%    testScenario(s,2.^ds,2.^rs,measures,mnames);
% end
% keyboard

basedir = '';%'/moto/nklab/users/az2522/';
simfolder = '20210203/';%'20210126/';%'20200929/';%'20201215/';%'20200929/';
simfolder = [basedir 'data/statsim/' simfolder];
if ~exist(simfolder,'dir'), mkdir(simfolder); end
simfile = @(s,d,r,n) sprintf('%s%ss%g_d%g_r%g_n%g.mat',simfolder,estimator,s,d,r,n);
% Initialize data for saving (could probably change from cell to tensor)
if ~STANDALONE % Collect all simulations, make figures
   sigtest   = repmat({false(nn,nsims       )},[nm nd nr]);
   bootsamps = repmat({  NaN(nsims,nboots)},[nm nd nr nn]);
   hyps      = repmat({repmat(Hypersphere(),[nn nsims])},[nm nd nr]);
end
% else do only one simulation, don't make figures

DISPLAYED = false;
FINISHED  = false;
for d = 1:nd
   for r = 1:nr
      % Generate points for the different scenarios
      gt     = generateScenario(s,2^ds(d),2^rs(r));
      nballs = numel(gt.radii);
      if s==3, nc2 = nballs;
      else     nc2 = nchoosek(nballs,2);
      end
      nc2c2  = nchoosek(max(2,nc2),2);

      for n = 1:nn
         stationarycounter([d r n],[nd nr nn])

         if strcmpi(testrslt{s},'jackknives')
            nboots = nballs*2^ns(n);
         end

         if exist(simfile(s,ds(d),rs(r),ns(n)),'file')
            try [hyp,gt,sigtmp,bootmp] = hyp_load(simfile,s,ds(d),rs(r),ns(n),DISPLAYED);
            catch fprintf('FILE CORRUPT\n'); continue
            end
         else
            if ~STANDALONE
               fprintf('FILE MISSING\n')
               continue
            end
            if ~DISPLAYED
               fprintf('Simulating points and estimating hyperspheres...      ')
            end
            [points,gt] = gt.sample(2.^[ns(n) ns(n)*ones(1,nballs-1)],nsims);
            % Simulate points
            hyp = Hypersphere.estimate(points,gt.categories,'independent',estimator,cvdists);
            % save in-progress fits
            savtmp = struct('hyp',hyp,'gt',gt,'n',ns(n),'d',ds(d),'r',rs(r));
            save(simfile(s,ds(d),rs(r),ns(n)),'-struct','savtmp')
         end

         if exist('sigtmp','var') && ~isempty(sigtmp)
            startSim = find(isnan(sigtmp(:,1)),1);
         end
         if ~exist('sigtmp','var') || isempty(sigtmp) || size(sigtmp,2)>1 || (~isempty(startSim) && startSim < nsims)
            if ~STANDALONE
               fprintf('FILE INCOMPLETE\n')
               continue
            end
            stationarycounter([d r n],[nd nr nn]);fprintf('       \n');
            if ~DISPLAYED
               fprintf('\nSignificance testing...      ')
               DISPLAYED = true;
            end
            
            if ~exist('points','var') % regenerate points
               gt.resetRandStream
               [points,gt] = gt.sample(2.^[ns(n) ns(n)*ones(1,nballs-1)],nsims);
            end
            % Assess significance on samples
            if ~exist('sigtmp','var') || isempty(sigtmp)
               sigtmp    = NaN(nsims,nc2c2);
               %sigptmp   = NaN(nsims,nc2c2);
               bootmp    = NaN(nsims,nc2,nboots);
            end
            startSim = find(isnan(sigtmp(:,1)),1);
            for b = startSim:nsims
               hyptmp = SetOfHyps(hyp(b)).significance(points(:,:,b),nboots,testtype{s},estimator,cvdists);
               sigtmp(b,:)  = hyptmp.(sigOrDiff).(mnames{s}(1:2));
               %sigptmp(b,:) = hyptmp.([sigOrDiff 'p']).(mnames{s}(1:2));
               bootmp(b,:,:) = cat(1,hyptmp.ci.(testrslt{s}).(mnames{s}))';
               stationarycounter(b,nsims)
               if mod(b,50)==0
                  %keyboard
                  savtmp = struct('hyp',hyp,'gt',gt,'sigtmp',sigtmp,'bootmp',bootmp,'n',ns(n),'d',ds(d),'r',rs(r));
                  save(simfile(s,ds(d),rs(r),ns(n)),'-struct','savtmp')
                  fprintf('\nSaved checkpoint          ')
               end
            end
            sigtmp = sigtmp(:,1+4*(s==6)); % pick 5th for s=6
            switch s
               case {1,2,7}; bootmp = bootmp(:,1,:);
               case {3,4,5}; bootmp = diff(bootmp(:, 1:2 ,:),[],2); % difference of first two measures
               case    6   ; bootmp = diff(bootmp(:,[1 6],:),[],2); % difference of first & last measures
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
            hyps{s,d,r}(n,:)    = hyp;
            groundtruth{s,d,r}  = gt;
            sigtest{s,d,r}(n,:) = sigtmp;
            bootsamps{s,d,r,n}  = bootmp;
         end
         clear sigtmp
      end
   end
end
if STANDALONE, fprintf('Simulations complete.\n'); return; end

%% Extract data: collect estimates (e.g., overlaps, distances)
estimates = repmat({NaN(nn,nsims  )},[nm nd nr]);
bootprc   = repmat({NaN(nn,nsims,2)},[nm nd nr]);
varestms  = NaN(nm,nd,nr,nn);
varboots  = repmat({NaN(nsims,1)},[nm nd nr nn]);
for d = 1:nd
   for r = 1:nr
      for n = 1:nn
         try
         switch s
            case {1,2,7}; estimates{s,d,r}(n,:) =             cat(1,hyps{s,d,r}(n,:).(mnames{s}));
            case {3,4,5}; estimates{s,d,r}(n,:) = diff(indexm(cat(1,hyps{s,d,r}(n,:).(mnames{s})),[],1:2),[],2);
            case    6   ; estimates{s,d,r}(n,:) = diff(indexm(cat(1,hyps{s,d,r}(n,:).(mnames{s})),[],[1 6]),[],2);
         end
         varestms(s,d,r,n) = var(estimates{s,d,r}(n,:));
         varboots{s,d,r,n} = var(bootsamps{s,d,r,n});
         if strcmp(sigOrDiff,'sig')
         bootprc{s,d,r}(n,:,:) = prctile([-Inf(nsims,1) squeeze(bootsamps{s,d,r,n})]',...
                                         [sigthresh*100 100])';
         else
         bootprc{s,d,r}(n,:,:) = prctile([-Inf(nsims,1) squeeze(bootsamps{s,d,r,n}) Inf(nsims,1)]',...
                                         [0 100] + [1 -1]*sigthresh*100/2)';
         end
         catch
         end
      end
   end
end
switch mnames{s} % normalization factor, depending on statistical test
   case {'overlap','margins'}, normfactor = @(r) min(2^rs(r).^[1 -1]);
   case  'dists',              normfactor = @(r) min(2^rs(r).^[0 -1]);
   case  'radii',              normfactor = @(r)     2^-rs(r);
end

%% Plot analyses
d = 1;
r = find(rs==0);
n = min(5,numel(ns));
nShown = [1 nn  1 nn];
dShown = [1  1 nd nd];
rShown = [1 nr  1 nr];

fp.DESTINATION = '1col'; fp.pap.sz = [7.5 3.5];
figsetup;
fp.dots{2} = 3; fp.nticks = [-1 4];

fh = newfigure([2 6],sprintf('%u%s',s,measures{s}),fp);
sampax = fh.a.h([1 7]);

showCIs(estimates{s,d,r}(n,:),squeeze(bootprc{s,d,r}(n,:,:)),sigtest{s,d,r}(n,:),...
        measures{s},2^ds(d),2^ns(n),fh.a.h(2))

axtivate(3)
if s==2, bootcentered = squeeze(bootsamps{s,d,r,n});
else     bootcentered = squeeze(bootsamps{s,d,r,n})-repmat(mean(bootsamps{s,d,r,n},3),[1 nboots]);
end
hb = histogram(bootcentered,100,'Normalization','pdf',...
               'Orientation','horizontal','FaceColor','none','EdgeColor',[0 0 0],'DisplayStyle','stairs');
he = histogram(estimates{s,d,r}(n,:),'BinEdges',hb.BinEdges,'Normalization','pdf',...
               'Orientation','horizontal','FaceColor','none','EdgeColor',[0.5 0.5 0.5],'DisplayStyle','stairs');
plot([0 max([hb.Values he.Values])],[0 0],'k-','Tag','phyaxwidth')
xlabel('pdf')
matchy(fh.a.h(2:3),'y')

showVar(varboots{s,d,r,n},varestms(s,d,r,n),2^ds(d),2^ns(n),'variance',fh.a.h(8))

axtivate(9)
for d = nd:-1:1
   for r = nr:-1:1
      for n = nn:-1:1
         plot(log10(varestms(s,d,r,n)),log10(mean(varboots{s,d,r,n})),'o','Color',[1 d*[1 1]/(nd+1)])
      end
   end
end
plot([min(axis) max(axis)],[min(axis) max(axis)],'k--')
xlabel('variance (full data)')
ylabel(sprintf('variance (%s)',testrslt{s}))
logAxis('xy')

if s~=3
for d = nd:-1:1
   plotErrorPatch(rs,squeeze(indexm(cat(3,estimates{s,d,:}),n)) ...
                     .*(ones(nsims,1)*arrayfun(normfactor,1:nr)),...
                  fh.a.h(4),[1 d*[1 1]/(nd+1)])%,'sem')
end
if s~=2, plot(rs([1 end]),[0 0],'k-','Tag','phyaxwidth'); end
xlim(rs([1 end])); xticks(rs); xticklabels(2.^rs);
ylabel(measures{s})
end

if s~=3
for r = nr:-1:1
   plotErrorPatch(ds,squeeze(indexm(cat(3,estimates{s,:,r}),n))*normfactor(r),...
                  fh.a.h(5),[r*[1 1]/(nr+1) 1])%,'sem')
end
r = find(rs==0);
else
for n = nn:-1:1
   plotErrorPatch(ds,squeeze(indexm(cat(3,estimates{s,:,r}),n))*normfactor(r),...
                  fh.a.h(5),[n*[1 1]/(nn+1) 1])%,'sem')
end
n = min(5,numel(ns));
end
if s~=2, plot(ds([1 end]),[0 0],'k-','Tag','phyaxwidth'); end
xlim(ds([1 end])); xticks(ds); xticklabels(2.^ds);
title(sprintf('%s, %u simulations each',measures{s},nsims))

for d = nd:-1:1
   plotErrorPatch(ns,cat(3,estimates{s,d,r})',...
                  fh.a.h(6),[1 d*[1 1]/(nd+1)])%,'sem')
end
if s~=2, plot(ns([1 end]),[0 0],'k-','Tag','phyaxwidth'); end
xlim(ns([1 end])); xticks(ns); xticklabels(2.^ns);
ylabel(measures{s})

if s~=3
axtivate(10)
for d = nd:-1:1
   a = permute(indexm(cat(3,sigtest{s,d,:}),n),[2 3 1]);
   plot(rs,100*mean(a),'-','LineWidth',2,'Color',[1 [d d]/(nd+1)])
end
plot(rs([1 end]),100*sigthresh*[1 1],'--k')
axis([rs([1 end]) 0 max(ylim)]); xticks(rs); xticklabels(2.^rs);
xlabel('radius ratio')
ylabel('false positive rate (%)')
end

axtivate(11)
if s~=3
for r = nr:-1:1
   a = permute(indexm(cat(3,sigtest{s,:,r}),n),[2 3 1]);
   plot(ds,100*mean(a),'-','LineWidth',2,'Color',[[r r]/(nr+1) 1])
end
r = find(rs==0);
else
for n = nn:-1:1
   a = permute(indexm(cat(3,sigtest{s,:,r}),n),[2 3 1]);
   plot(ds,100*mean(a),'-','LineWidth',2,'Color',[[n n]/(nn+1) 1])
end
n = min(5,numel(ns));
end
plot(ds([1 end]),100*sigthresh*[1 1],'--k')
axis([ds([1 end]) 0 max(ylim)]); xticks(ds); xticklabels(2.^ds);
xlabel('dimensions')
title(sprintf('false positive rate, %u simulations each',nsims))

axtivate(12)
for d = nd:-1:1
   a = cat(3,sigtest{s,d,r})';
   plot(ns,100*mean(a),'-','LineWidth',2,'Color',[1 [d d]/(nd+1)])
end
plot(ns([1 end]),100*sigthresh*[1 1],'--k')
axis([ns([1 end]) 0 max(ylim)]); xticks(ns); xticklabels(2.^ns);
xlabel('samples')

matchy(fh.a.h(10:12),'y')

%keyboard
batchPlotRefine(fh,fp);
%set(findall(fh.a.h(2).Children,'Color',linecol),'LineWidth',4)
hideAxes(sampax,'x','y');
hideAxes(sampax,'x'); %% TODO: need to fix above line
hideAxes(fh.a.h([3 5 11 12]),'y');
hideAxes(fh.a.h(4:6),'x');
printFig;


fh = newfigure([1 1]*ceil(sqrt(numel(ns))),fh,sprintf('%u%sCIs',s,measures{s}));
for n = 1:numel(ns)
   showCIs(estimates{s,d,r}(n,:),squeeze(bootprc{s,d,r}(n,:,:)),sigtest{s,d,r}(n,:),...
           measures{s},2^ds(d),2^ns(n),fh.a(end).h(n))
end

fh = newfigure([1 1]*ceil(sqrt(numel(ns))),fh,sprintf('%u%svars',s,measures{s}));
for n = 1:numel(ns)
   showVar(varboots{s,d,r,n},varestms(s,d,r,n),2^ds(d),2^ns(n),'variance',fh.a(end).h(n))
end

fh = newfigure([1 1]*ceil(sqrt(numel(ns))),fh,sprintf('%u%svarsIndividual',s,measures{s}));
% find 0:9:100th percentile estimates
n = numel(ns);
y = prctile(estimates{s,d,r}(n,:),linspace(0,100,numel(fh.a(end).h)));
for i = 1:numel(fh.a(end).h)
   yix = minix(abs(estimates{s,d,r}(n,:)-y(i)));
   showVar(squeeze(bootsamps{s,d,r,n}(yix,:,:)),estimates{s,d,r}(n,yix),...
           2^ds(d),2^ns(n),measures{s},fh.a(end).h(i));
end
printFig(fh.f(2:end))

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
      case {1,3,7} % OVERLAP. Null distribution: 0 overlap.
         hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1)],[r 1]);
      case 2 % INTER-CENTER DISTANCE. Null distribution: 0 distance.
         hyp = Hypersphere( zeros(2,d),[r 1]);
      case {4,5} % (Need 3+ balls.)
         % RADIUS DIFFERENCE. Null distribution: radii are the same.
         % OVERLAP DIFFERENCE. Null distribution: overlaps are the same.
         % INTER-CENTER DISTANCE DIFFERENCE. Null distribution: distances are the same.
         hyp = Hypersphere([zeros(1,d);r+0.75 zeros(1,d-1);0 r+0.75 zeros(1,d-2)],[r 1 1]);
         % hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1);ones(1,d)],[r 1 1]);
      % case 4 % OVERLAP DIFFERENCE. Null distribution: overlaps are the same. (Need 3+ balls.)
         % hyp = Hypersphere([zeros(1,d);3*r/2 zeros(1,d-1);0 3*r/2 zeros(1,d-2)],[r 1 1]);
      case 6
         hyp = Hypersphere([[-1 -r;-1 0.75; 1 r; 1 -0.75] zeros(4,d-2)],[r 1 r 1]);
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

function showCIs(estimates,bootprc,sigtest,measure,d,n,ax)
   if ~exist('ax','var') || isempty(ax), ax = gca; end
   axtivate(ax);

   nsims = numel(estimates);
   for i = 1:nsims
      if sigtest(i), linecol = [0   0   0  ];
      else           linecol = [0.5 0.5 0.5];
      end
      plot([i i],bootprc(i,:),'-','Color',linecol,'LineWidth',2)
      plot(i,estimates(i),'wo','MarkerFaceColor',linecol,'LineWidth',1.5)
   end
   plot([0 nsims],[0 0],'k-')
   ylabel(sprintf('%s estimates',measure))
   xlabel('simulation #')
   title(sprintf('all simulations (%ud, %u samples/ea)',d,n))
end

function showVar(varboots,varestms,d,n,xlab,ax)
   if ~exist('ax','var') || isempty(ax), ax = gca; end
   axtivate(ax);

   hb = histogram(varboots,'Normalization','probability','FaceColor',[0 0 0],'EdgeColor','none');
   plot([1 1]*varestms,[0 0.25],'w-' ,'LineWidth',3)
   plot([1 1]*varestms,[0 0.25],'k--','LineWidth',2)
   ylabel('probability')
   xlabel(xlab)
   title(sprintf('(%ud, %u samples/ea)',d,n))
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


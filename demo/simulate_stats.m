function [sigtest,estimates,hyps,groundtruth] = simulate_stats(nsamps)
% simulate_stats: simulate null distributions from different h2s statistical
%    tests with a few special cases

if ~exist('nsamps','var') || isempty(nsamps), nsamps = 100; end

%% TODO: ALL TESTS WITH DIFFERENT ESTIMATORS?
nboots = 100;   % # of bootstraps during statistical testing
%% Testing as a function of d, n, n-ratio, and r-ratio:
ns = 2.^( 6:8); % # samples to test
ds = 2.^( 1:5); % # dimensions
rs = 2.^(-1:3); % # radii

r  = 2;

sigfield = {'sig',     'sigdiff'};
measures = {'overlap', 'distance', 'radius difference', ...
            'overlap difference', 'distance difference'};
mnames   = {'overlap', 'distances2CrossValidated', 'radii', ...
            'overlap', 'distances2CrossValidated'};
% % debug
% for s = 1:5
% 	testScenario(s,ds,rs,measures,mnames);
% end
% keyboard

nn = numel(ns);
nd = numel(ds);
nm = numel(measures);

s  = 3;         % scenario chosen from the following:
% Different scenarios needing testing:
% (1) 2 balls with 0 overlap.
% (2) 2 balls with 0 distance.
% (3) 2 balls with same (non-zero) radius.
% (4) 3 balls, 2 with same (non-zero) overlaps
% (5) 3 balls, 2 with same (non-zero) distances
% Can combine scenarios:
% (a) is (1) and (2) in 3 balls. Balls 1 and 2 have 0 distance, ball 3 placed
%    to have 0 overlap with ball 1.
% ... though this gets annoying if playing with r- and n-ratios.

% Initialize data for saving (could probably change from cell to tensor)
nc2     = 1+2*double(s>2);
sigtest = repmat({NaN(nn,nsamps)},[nm nd]);
simfile = sprintf('statsim_%u%s.mat',s,measures{s}(1:2));

if exist(simfile,'file'), load(simfile); fprintf('Loaded simulations.\n'); end

DISPLAYED = false;
for d = 1:nd
	% Generate points for the different scenarios
	groundtruth{s,d,r} = generateScenario(s,ds(d),rs(r));

   for n = 1:nn
   	if any(isnan(vectify(sigtest{s,d}(n,:,:))))
   		if ~DISPLAYED
   			fprintf('Simulating points, estimating hyperspheres, and significance testing...      ')
   			DISPLAYED = true;
   		end
		   % Simulate points
	      [points,groundtruth{s,d,r}] = groundtruth{s,d,r}.sample([ns(n) ns(n)*ones(1,1+double(s>2))],nsamps);
	      hyps{s}(d,n,:) = Hypersphere('estimate',points,groundtruth{s,d,r}.categories,'independent');%.meanAndMerge(true);

	      % Assess significance on samples
	      tmp = NaN(nsamps,nc2);
	      parfor b = 1:nsamps
	         tmp(b,:) = SetOfHyps(hyps{s}(d,n,b)).significance(...
	                            points(:,:,b),nboots).(sigfield{1+double(s>2)}).(mnames{s}(1:2));
	      end
	      sigtest{s,d}(n,:) = tmp(:,1+2*double(s==3)); % pick 3rd comparison only when s==3
	      stationarycounter([d n],[nd nn])
	      % save in-progress fits
	      save(simfile,'hyps','sigtest','groundtruth')
	   end
   end
end

%% Extract data: collect estimates (e.g., overlaps, distances)
estimates = repmat({NaN(nn,nsamps)},[nm nd]);
for d = 1:nd
	for n = 1:nn
      switch s
      	case {1,2}; estimates{s,d}(n,:) =                     hyps{s}(d,n,:).(mnames{s});
      	case 3;     estimates{s,d}(n,:) = diff(indexm(vertcat(hyps{s}(d,n,:).(mnames{s})),[],2:3),[],2);
   		case {4,5}; estimates{s,d}(n,:) = diff(indexm(vertcat(hyps{3}(d,n,:).(mnames{s})),[],1:2),[],2);
   	end
   end
end

%% Plot analyses
nShown = [1 nn  1 nn];
dShown = [1  1 nd nd];
rShown = 2;

fh = newfigure([3 3],measures{s});
ax = fh.a.h([4 5 7 8]);

for i = 1:3
   sampsz = (2^(i+4))*ones(1,2+double(s>2));
   generateScenario(s,2,rs(rShown)).plotSamples(sampsz,fh.a.h(i));
   if s<3, title(sprintf('Samples: n_1= %u, n_2= %u',         sampsz))
   else    title(sprintf('Samples: n_1= %u, n_2= %u, n_3= %u',sampsz))
   end
end

%% TODO: different colors for different conditions

for i = 1:4
   axtivate(ax(i));
   hist(estimates{s,dShown(i)}(nShown(i),:))
   set(ax(i).Children(end),'FaceColor',[0 0 0],'EdgeColor',[1 1 1])
   plot([0 0],[0 nsamps/4],'r-','LineWidth',2)
   xlabel(measures{s})
   ylabel('count')
   title(sprintf('Estimates (n = %u, d = %u)',ns(nShown(i)),ds(dShown(i))))
end
matchy(ax)

for i = 2:3
   plotErrorPatch(log2(ns),estimates{s,dShown(i)}',fh.a.h(6),[0 0 0])%,'sem')
end
logAxis(2)
xlabel('samples')
ylabel(measures{s})
title('Estimate vs n')

axtivate(9);
for i = 2:3
   plot(log2(ns),100*mean(1-sigtest{s,dShown(i)}<=0.05,2),'-k','LineWidth',2)
end
plot(log2(ns([1 end])),[5 5],'--k')
logAxis(2)
xlabel('samples')
ylabel('False positive rate (%)')
title('False positive rate vs n')

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
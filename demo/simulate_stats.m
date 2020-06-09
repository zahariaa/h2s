function [sigtest,estimates,hyps,groundtruth] = simulate_stats(nboots)
% simulate_stats: simulate null distributions from different h2s statistical
%    tests with a few special cases

if ~exist('nboots','var') || isempty(nboots), nboots = 100; end

%% TODO: ALL TESTS WITH DIFFERENT ESTIMATORS?
%% Testing as a function of d, n, n-ratio, and r-ratio:
ns = 2.^( 6:8); % # samples to test
ds = 2.^( 1:5); % # dimensions
rs = 2.^(-1:3); % # radii
s  = 3;         % scenario chosen from the following:
% Different scenarios needing testing:
% (1) 2 balls with 0 overlap.
% (2) 2 balls with 0 distance.
% (3) 2 balls with same (non-zero) radius.
% (4) 3 balls, 2 with same (non-zero) overlaps
% (5) 3 balls, 2 with same (non-zero) distances

r  = 2; 

sigfield = {'sig',     'sigdiff'};
measures = {'overlap', 'distance', 'radius difference', ...
            'overlap difference', 'distance difference'};
mnames   = {'overlap', 'distances2CrossValidated', 'radii', ...
            'overlap', 'distances2CrossValidated'};
simfile  = sprintf('statsim_%u%s.mat',s,measures{s}(1:2));

nn = numel(ns);
nd = numel(ds);
nm = numel(measures);

% Can combine scenarios:
% (a) is (1) and (2) in 3 balls. Balls 1 and 2 have 0 distance, ball 3 placed
%    to have 0 overlap with ball 1.
% ... though this gets annoying if playing with r- and n-ratios.

% Initialize data for saving (could probably change from cell to tensor)
nc2   	 = 1+2*double(s>2);
sigtest   = repmat({NaN(nn,nboots,nc2)},[nm nd]);

if exist(simfile,'file'), load(simfile); fprintf('Loaded simulations.\n'); end

DISPLAYED = false;
for d = 1:nd
	% Generate points for the different scenarios
	groundtruth = generateScenario(s,ds(d),rs(r));

   for n = 1:nn
   	if any(isnan(vectify(sigtest{s,d}(n,:,:))))
   		if ~DISPLAYED
   			fprintf('Simulating points, estimating hyperspheres, and significance testing...      ')
   			DISPLAYED = true;
   		end
		   % Simulate points
	      [points,groundtruth] = groundtruth.sample([ns(n) ns(n)*ones(1,1+double(s>2))],nboots);
	      hyps(d,n,:) = Hypersphere('estimate',points,groundtruth.categories,'independent');%.meanAndMerge(true);

	      % Assess significance on samples
	      tmp = NaN(size(sigtest{s,d},[2 3]));
	      parfor b = 1:nboots
	         tmp(b,:) = SetOfHyps(hyps(d,n,b)).significance(...
	                            points(:,:,b),100).(sigfield{1+double(s>2)}).(mnames{s}(1:2));
	      end
		   sigtest{s,d}(n,:,:) = tmp;
	      stationarycounter([d n],[nd nn])
	      % save in-progress fits
	      save(simfile,'hyps','sigtest')
	   end
   end
end

%% Extract data: collect estimates (e.g., overlaps, distances)
estimates = repmat({NaN(nn,nboots,nc2)},[nm nd]);
for d = 1:nd
	for n = 1:nn
      switch s
      	case {1,2};   estimates{s,d}(n,:) =              hyps(d,n,:).(mnames{s});
      	case {3,4,5}; estimates{s,d}(n,:) = diff(vertcat(hyps(d,n,:).(mnames{s})),[],2);
   	end
   end
end

%% Plot analyses
nShown = [1 nn  1 nn];
dShown = [1  1 nd nd];

fh = newfigure([3 3],measures{s});
ax = fh.a.h([4 5 7 8]);

for i = 1:3
   sampsz = 2.^[i+4 i+4];
   groundtruth.plotSamples(sampsz,fh.a.h(i));
   title(sprintf('Samples: n_1= %u, n_2= %u',sampsz))
end

%% TODO: different colors for different conditions

for i = 1:4
   axtivate(ax(i));
   hist(estimates{s,dShown(i)}(nShown(i),:))
   set(ax(i).Children(end),'FaceColor',[0 0 0],'EdgeColor',[1 1 1])
   plot([0 0],[0 nboots/4],'r-','LineWidth',2)
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
      	hyp = Hypersphere( zeros(2,d),                  [r 1]);
   	case 3 % RADIUS DIFFERENCE. Null distribution: radii are the same. (Need 3+ balls.)
   		hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1);ones(1,d)],[r 1 1]);
   	% case 4 % OVERLAP DIFFERENCE. Null distribution: overlaps are the same. (Need 3+ balls.)
   	% 	hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1);ones(1,d)],[r 1 1]);
   	% case 5 % INTER-CENTER DISTANCE DIFFERENCE. Null distribution: distances are the same.
   	% 	    %    (Need 3+ balls.)
   	% 	hyp = Hypersphere([zeros(1,d);r+1 zeros(1,d-1);ones(1,d)],[r 1 1]);
   end

end
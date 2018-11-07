% demo_dynamic: generates an H2S movie based on temporal RDM
% 
% 2018-10-16 AZ Created

dimLow = 3;

%% LOAD DATA
datafile = 'RDMs_srd_weighted_S2';
X = load(['data' filesep datafile '.mat']);
X = X.(cellify(fields(X)));

%% Define categories
labels = {'faces','fruits','places','bodyparts','objects'};
cats = Categories([20 20 20 20 20],labels);
times = -100:700;

%% Compute!
h2ssavefile = ['data' filesep datafile '-h2s.mat'];
if exist(h2ssavefile,'file'), load(h2ssavefile);
else
   hi = SetOfHyps('estimate',X(:,:,1:end),cats,1);
   lo = hi.h2s(dimLow);
%% SAVE
   save(h2ssavefile,'hi','lo');
end

%% Run movie
figure; lo.movie(times,'ms',datafile)

%% Generate plots
figure;
subplot(3,1,2); hi.plotComparisons('margins',times);
ylabel('separations')
subplot(3,1,3); hi.plotComparisons('dists',times);
ylabel('distances')
xlabel('time (ms)');
subplot(3,1,1); plot(times,reshape([hi.radii],5,[])','LineWidth',1);
ylabel('spreads')
box off; xlim([min(times) max(times)]);
title(datafile,'Interpreter','none')

cats.legend([0.01 0.9 1 0.1])
%% Print
set(gcf,'Name','comparisons','Renderer','painters');
printFig([],[],'eps')


% demo_dynamic: generates an H2S movie based on temporal RDM
% 
% 2018-10-16 AZ Created

dimLow = 3;

%% LOAD DATA
datafile = 'RDMs_srd_weighted_S2';
X = load(['data' filesep datafile '.mat']);
X = X.(cellify(fields(X))); % de-struct

%% Define categories
labels = {'faces','fruits','places','bodyparts','objects'};
cats = Categories({20 20 20 20 20},labels,lines(numel(labels)));
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
figure; lo.movie(times,'ms',[datafile '.avi'])

%% Generate plots
fh = newfigure([3 1],'comparisons');
hi.plotDynamics('all',times,fh.a.h,true);
xlabel('time (ms)');
axtivate(fh.a.h(1));
title(datafile,'Interpreter','none')

%% Print
set(gcf,'Renderer','painters');
printFig(fh,[],'eps')


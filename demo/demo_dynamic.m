% demo_dynamic: generates an H2S movie based on temporal RDM
% 
% 2018-10-16 AZ Created

dimLow = 3;

%% LOAD DATA
datafile = 'RDMs_srd_weighted_S2';
X = load(['data/' datafile '.mat']);
X = X.(cellify(fields(X)));

%% Define categories
labels  = {'faces','fruits','places','bodyparts','objects'};
vectors = false(size(X,1),numel(labels));
for i = 1:numel(labels)
   vectors((i-1)*20+1:20*i,i) = true;
end
cats = Categories(vectors,labels);


%% Compute!
h = SetOfHyps('estimate',X(:,:,1:end),cats,1);
h = h.h2s(dimLow);
times = -100:700;

%% Run movie
figure; h.movie(times)


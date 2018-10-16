% demo_dynamic: generates an H2S movie based on temporal RDM
% 
% 2018-10-16 AZ Created


%% LOAD DATA
X = load('data/RDMs_srd_weighted_S1.mat');
X = X.(cellify(fields(X)));

%% Define categories
labels  = {'faces','fruits','places','bodyparts','objects'};
vectors = false(size(X,1),numel(labels));
for i = 1:numel(labels)
   vectors((i-1)*20+1:20*i,i) = true;
end
cats = Categories(vectors,labels);


%% Do stuff
figure;
for i = 1:10:size(X,3)
   h = SetOfHyps('estimate',X(:,:,i),cats).h2s;
   cla;h.show;
   title(num2str(i))
   drawnow
end



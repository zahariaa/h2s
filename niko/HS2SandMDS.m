function HS2SandMDS(pointsOrDistMat,categories,ax,titleStr,dimLow)
% create hypersphere-to-sphere visualization and show it alongside the MDS solution for comparison 
% 
% pass 4 axis handles for: h2s, PCA, MDS, t-SNE
%      3              for: h2s,      MDS, t-SNE
%      2              for: h2s,      MDS
%
% NK & AZ 2017-8

%%  preparations
if ~exist('dimLow','var'), dimLow = 3; end
[nPoints,nCats] = size(categories.vectors);
nax = numel(ax);

%% hypersphere to sphere visualization
model = SetOfHyps('estimate',pointsOrDistMat,categories).h2s(dimLow);
model.show(ax(1),titleStr);

%% other visualizations
dists = pdist(pointsOrDistMat,'Euclidean');
points2D{floor(nax/2)} = mdscale(dists,2,'criterion','metricstress');
% only run requested visualizations
if nax>=3,   points2D{nax-1} = tsne(pointsOrDistMat,[],[],min(size(pointsOrDistMat)));   end
if nax==4,   points2D{1} = simplepca(pointsOrDistMat,dimLow);   end

%% Draw points from circle/sphere to align other visualizations to
modelPointsAlign = NaN(0,dimLow);
nPointsPerCat = sum(categories.vectors);
% if all vectors are less than nPoints, add more proportionally
if sum(nPointsPerCat)<nPoints
   nPointsPerCat = nPointsPerCat + floor((nPoints-sum(nPointsPerCat))*nPointsPerCat/sum(nPointsPerCat));
   catI = 1;
   while sum(nPointsPerCat)<nPoints
      nPointsPerCat(catI) = nPointsPerCat(catI) + 1;
      catI=catI+1;
   end
end
for catI = 1:nCats
   modelPointsAlign = [modelPointsAlign;
                       randsphere(nPointsPerCat(catI),dimLow,model.radii(catI)) ...
                       + repmat(model.centers(catI,:),[nPointsPerCat(catI) 1])];
end

%% Align PCA, MDS, and t-SNE relative to H2S
if size(pointsOrDistMat,2) > 1
   for i = 1:numel(points2D)
      [~,points2D{i}] = procrustes(modelPointsAlign,points2D{i});
   end
end

%% PLOT
for i = 1:nax-1
   axtivate(ax(i+1));
   ms=5;
   for catI = 1:nCats
       plot(points2D{i}(logical(categories.vectors(:,catI)),1), ...
            points2D{i}(logical(categories.vectors(:,catI)),2),...
           'o', 'MarkerFaceColor',categories.colors(catI,:), ...
           'MarkerEdgeColor','none', 'MarkerSize',ms);
   end
   axis equal off;
end

% showRDMs(dists,300,0);

% if nCats==2
%     dists_sq = squareform(dists);
%     withinCat1dists = dists_sq(1:end/2,1:end/2);
%     withinCat1dists = withinCat1dists(~logical(eye(size(dists_sq,1)/2)));
%     withinCat2dists = dists_sq(end/2+1:end,end/2+1:end);
%     withinCat2dists = withinCat2dists(~logical(eye(size(dists_sq,1)/2)));
% 
%     betweenCat1dists = dists_sq(1:end/2,end/2+1:end);
% 
%     figure(301); set(301,'Color','w'); clf;
%     hist([withinCat1dists withinCat2dists betweenCat1dists(1:numel(withinCat1dists))']);
%     legend({'within cat 1','within cat 2','between cat'});
% end



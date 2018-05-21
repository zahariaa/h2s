function HS2SandMDS(pointsOrDistMat,categories,ax,titleStr,dimLow)
% create hypersphere-to-sphere visualization and show it alongside the MDS solution for comparison 


%%  preparations
if ~exist('dimLow','var'), dimLow = 3; end
[nPoints,nCats] = size(categories.vectors);

%% hypersphere to sphere visualization
model = hypersphere2sphere(pointsOrDistMat,categories,[],dimLow);
if     dimLow==2;   showCirclesModel(model, ax(1), titleStr);
elseif dimLow==3;   showSpheresModel(model, ax(1), titleStr);
end

%% multidimensional scaling visualization
dists = pdist(pointsOrDistMat,'Euclidean');

points2D{1} = simplepca(pointsOrDistMat,dimLow);
points2D{2} = mdscale(dists,2,'criterion','metricstress');
points2D{3} = tsne(pointsOrDistMat,[],[],min(size(pointsOrDistMat)));

% Draw points from circle/sphere to align other visualizations to
modelPointsAlign = NaN(0,dimLow);
for catI = 1:nCats
   modelPointsAlign = [modelPointsAlign;
                       randsphere(nPoints/nCats,dimLow,model.radii(catI)) ...
                       + repmat(model.centers(catI,:),[nPoints/nCats 1])];
end
% Align PCA, MDS, and t-SNE relative to H2S
if size(pointsOrDistMat,2) > 1
   for i = 1:numel(points2D)
      [~,points2D{i}] = procrustes(modelPointsAlign,points2D{i});
   end
end

for i = 1:numel(points2D)
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



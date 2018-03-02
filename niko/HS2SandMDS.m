function HS2SandMDS(pointsOrDistMat,categories,figPanelSpec,titleStr,dimLow)

% create hypersphere-to-sphere visualization and show it alongside the MDS solution for comparison 


%%  preparations
if ~exist('dimLow','var'), dimLow = 3; end
[~,nCats] = size(categories.vectors);

%% hypersphere to sphere visualization
model = hypersphere2sphere(pointsOrDistMat,categories,[],dimLow);
if     dimLow==2;   [~,XY] = showCirclesModel(model, figPanelSpec, titleStr,ceil(nPoints/nCats));
elseif dimLow==3;   showSpheresModel(model, figPanelSpec, titleStr);
end

%% multidimensional scaling visualization
dists = pdist(pointsOrDistMat,'Euclidean');

points2D{1} = simplepca(pointsOrDistMat,dimLow);
points2D{2} = mdscale(dists,2,'criterion','metricstress');
points2D{3} = tsne(pointsOrDistMat,[],[],min(size(pointsOrDistMat)));

% Align MDS and t-SNE relative to PCA
if size(pointsOrDistMat,2) > 1
   [~,points2D{1}] = procrustes(XY,points2D{1});
   [~,points2D{2}] = procrustes(XY,points2D{2});
   [~,points2D{3}] = procrustes(XY,points2D{3});
end

for i = 1:numel(points2D)
   subplot(figPanelSpec(2), figPanelSpec(3), figPanelSpec(4)+i); hold on;
   ms=5;
   for catI = 1:nCats
       plot(points2D{i}(logical(categories.vectors(:,catI)),1), points2D{i}(logical(categories.vectors(:,catI)),2),...
           'o', 'MarkerFaceColor',categories.colors(catI,:), 'MarkerEdgeColor','none', 'MarkerSize',ms);
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



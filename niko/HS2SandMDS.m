function HS2SandMDS(pointsOrDistMat,categories,figPanelSpec,titleStr,psOutput,pdfOutput)

% create hypersphere-to-sphere visualization and show it alongside the MDS solution for comparison 


%%  preparations

[~,nCats] = size(categories.vectors);

%% hypersphere to sphere visualization
model = hypersphere2sphere(pointsOrDistMat,categories);
showSpheresModel(model, figPanelSpec, titleStr)


%% multidimensional scaling visualization
dists = pdist(pointsOrDistMat,'Euclidean');

mdsPoints_2d = mdscale(dists,2,'criterion','metricstress');

subplot(figPanelSpec(2), figPanelSpec(3), figPanelSpec(4)+1); hold on;
ms=5;
for catI = 1:nCats
    plot(mdsPoints_2d(logical(categories.vectors(:,catI)),1), mdsPoints_2d(logical(categories.vectors(:,catI)),2),...
        'o', 'MarkerFaceColor',categories.colors(catI,:), 'MarkerEdgeColor','none', 'MarkerSize',ms);
end
axis equal off;


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



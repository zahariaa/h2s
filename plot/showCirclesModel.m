function [objHs,XY] = showCirclesModel(model, figPanelSpec, titleStr,nCircPts)

%% control variables
if ~exist('nCircPts','var') || isempty(nCircPts),   nCircPts = 100;   end
opacity = 0.2;

%% preparations
if numel(figPanelSpec)==4,   figurePanel(figPanelSpec);
elseif numel(figPanelSpec)==1 && ishandle(figPanelSpec), axtivate(figPanelSpec)
end
if ~exist('titleStr','var')
    titleStr = {'\bfhypersphere2sphere',any2str('error <= ',ceil(model.error*100),'%')};
end
title(titleStr);

%% draw circles
[objHs,XY] = drawCircle(model.centers(:,1:2),model.radii,model.categories.colors,[],opacity,nCircPts);
% xlabel({any2str('# samples/category: ',model.nSamplesPerCat)});
axis equal tight off;

return


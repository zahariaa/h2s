function showSpheresModel(model, figPanelSpec, titleStr)

%% control variables
patchDetail = 100;
opacity = 0.5;


%% preparations
if numel(figPanelSpec)==4,   figurePanel(figPanelSpec);
elseif numel(figPanelSpec)==1 && ishandle(figPanelSpec), axtivate(figPanelSpec)
end

if ~exist('titleStr','var')
    titleStr = {'\bfhypersphere2sphere',any2str('error <= ',ceil(model.error*100),'%')}
end
title(titleStr);


%% draw spheres
[~,nCats] = size(model.categories.vectors);
covs = arrayfun(@(r) eye(3)*r^2,model.radii,'UniformOutput',false);
objHs = draw3dEllipsoid(model.centers,covs,model.categories.colors,patchDetail,opacity);

axis equal tight off;
set(gcf,'Renderer','OpenGL');


%% set view and lighting
% set view angle orthogonal to the PC12 plane of the locations
eigenvecs = princomp(model.centers);
normal = null(eigenvecs(:,1:2)');

if corr(model.centers*normal,model.radii')<0
    normal = -normal;
end

objectsCentroid = mean(model.centers,1);
camDist = 50*sqrt(sum(std(model.centers).^2));
camPos = objectsCentroid'-normal*camDist;

set(gca,'CameraPosition',camPos,'CameraTarget',objectsCentroid);

%delete(findall(gca,'type','light'));

set(objHs,'AmbientStrength',  0,...
          'SpecularStrength', 1,...
          'DiffuseStrength',  1);


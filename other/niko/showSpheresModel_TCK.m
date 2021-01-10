function showSpheresModel(model, figPanelSpec, titleStr)

%% control variables
patchDetail = 100;
opacity = 0.7;


%% preparations
figurePanel(figPanelSpec);
if ~exist('titleStr','var')
    titleStr = {'\bfhypersphere2sphere',any2str('error <= ',ceil(model.error*100),'%')}
end
title(titleStr);


%% draw spheres
[nPoints nCats] = size(model.categories.vectors);
objHs = nan(nCats,1);
for catI = 1:nCats
    objHs(catI) = draw3dEllipsoid(model.centers(catI,:),eye(3)*model.radii(catI)^2,model.categories.colors(catI,:),patchDetail,opacity);
    % plot3(model.samples(:,1,catI),model.samples(:,2,catI),model.samples(:,3,catI),'.k');
end
% xlabel({any2str('# samples/category: ',model.nSamplesPerCat)});

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

%hack to have fixed viewpoint
camPos= [-17.9148, -13.0353, 143.1015];
objectsCentroid = 1.0e-15.*[-0.2220,0.0139,0];
set(gca,'CameraPosition',camPos,'CameraTarget',objectsCentroid);

%camlight('headlight');
camlight('left');
%camlight('right');
lighting phong;
%delete(findall(gca,'type','light'));

set(objHs,'AmbientStrength',  0,...
          'SpecularStrength', 1,...
          'DiffuseStrength',  1);


      
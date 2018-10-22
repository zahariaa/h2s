function showSpheresModel(varargin)

%% control variables
patchDetail = 100;
opacity = 0.5;

%% preparations
SETCAMERA = true;
model = varargin{1};
vv = 1;
for v = 2:nargin
   if     islogical(varargin{v}),                   SETCAMERA = varargin{v};
   elseif ishandle( varargin{v}) && numel(varargin{v})==1, ax = varargin{v};
   elseif isnumeric(varargin{v}) && numel(varargin{v})==4
      ax = figurePanel(varargin{v});
   else
      if vv==1 && ~exist('ax','var'), ax       = varargin{v}; vv=2;
      else                            titleStr = varargin{v};
      end
   end
end

if ~exist('ax','var') || isempty(ax), ax = gca; end
axtivate(ax);

if ~exist('titleStr','var')
    titleStr = any2str('error <= ',ceil(model.error*100),'%');
end
title(titleStr);


%% draw spheres
[~,nCats] = size(model.categories.vectors);
covs = arrayfun(@(r) eye(3)*r^2,model.radii,'UniformOutput',false);
objHs = draw3dEllipsoid(model.centers,covs,model.categories.colors,patchDetail,opacity);

axis equal tight off;
set(get(ax,'Parent'),'Renderer','OpenGL');


%% set view and lighting
% set view angle orthogonal to the PC12 plane of the locations
eigenvecs = pca(model.centers);
if size(eigenvecs,2)<2, normal = [0 0 0]';
else                    normal = null(eigenvecs(:,1:2)');
end

if corr(model.centers*normal,model.radii')<0
    normal = -normal;
end

objectsCentroid = mean(model.centers,1);
camDist = 50*sqrt(sum(std(model.centers).^2));
camPos = objectsCentroid'-normal*camDist;

set(gca,'CameraPosition',camPos,'CameraTarget',objectsCentroid);
% Force lighting to camlight headlight
delete(findall(gca,'type','light')); camlight
set(objHs,'AmbientStrength',  0,...
          'SpecularStrength', 1,...
          'DiffuseStrength',  1);


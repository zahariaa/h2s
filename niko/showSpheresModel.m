function showSpheresModel(varargin)

%% control variables
patchDetail = 100;
opacity = 0.4;

%% preparations
SETCAMERA = true;
model = varargin{1};
vv = 1;
for v = 2:nargin
   if     islogical(varargin{v}),    SETCAMERA = varargin{v};
   elseif  ~isempty(varargin{v}) && ishandle(varargin{v}) && numel(varargin{v})==1
      ax = varargin{v};
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
covs = arrayfun(@(r) eye(3)*r^2,model.radii,'UniformOutput',false);
objHs = draw3dEllipsoid(model.centers,covs,model.categories.colors,patchDetail,opacity);

axis equal tight off;
set(get(ax,'Parent'),'Renderer','OpenGL');


%% set view and lighting
if SETCAMERA,  SetOfHyps(model.centers,model.radii).camera;   end

return
end


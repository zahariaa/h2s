function [objHs,XY] = showModel(varargin)

model  = varargin{1};
dimLow = size(model.centers,2);

%% control variables
patchDetail = 100;
if dimLow == 2
  opacity = 0.2;
elseif dimLow == 3
  opacity = 0.5;
  SETCAMERA = true;
end

%% preparations
vv = 1;
for v = 2:nargin
   if     islogical(varargin{v}),    SETCAMERA = varargin{v};
   elseif  ~isempty(varargin{v}) && ishandle(varargin{v}) && numel(varargin{v})==1
      ax = varargin{v};
   elseif isnumeric(varargin{v})
      if     numel(varargin{v})==4,   ax = figurePanel(varargin{v});
      elseif numel(varargin{v})==1,   patchDetail = varargin{v};
      end
   else
      if vv==1 && ~exist('ax','var'), ax       = varargin{v}; vv=2;
      else                            titleStr = varargin{v};
      end
   end
end

if ~exist('ax','var') || isempty(ax), ax = gca; end
axtivate(ax);

if ~exist('titleStr','var')
   maxerror = max(model.error)/max([model.radii model.dists]);
   titleStr = any2str('error <= ',ceil(maxerror*100),'%');
end
title(titleStr);

if dimLow == 2
  %% draw circles
  [objHs,XY] = drawCircle(model.centers(:,1:2),model.radii,model.categories.colors,[],opacity,patchDetail);
  % xlabel({any2str('# samples/category: ',model.nSamplesPerCat)});
elseif dimLow == 3
  %% draw spheres
  covs = arrayfun(@(r) eye(3)*r^2,model.radii,'UniformOutput',false);
  [objHs,XY] = draw3dEllipsoid(model.centers,covs,model.categories.colors,patchDetail,opacity);
  set(get(ax,'Parent'),'Renderer','OpenGL');
end

axis equal tight off;


%% set view and lighting
if dimLow == 3 && SETCAMERA,  model.camera;   end

return
end

function [objHs,XY] = showModel(varargin)
% showModel: visualizes a Hypersphere or SetOfHyps object as a set of circles
%    (for a 2D input) or spheres (for 3D input). If the dimensionality of the
%    input is larger than 3.
% 
% e.g.:
% showModel(model)
% showModel(model,'A very nice model')
% showModel(model,500)                % super hi-res!
% showModel(model,false)              % don't use standard camera angle
% [objHs,XY] = showModel(model);      % keep object handles and corrdinates
% model.show;                         % works for Hypersphere/SetOfHyps objects
% 
% Required (first) input argument:
%    model: a Hypersphere or SetOfHyps object.
% 
% Optional inputs:
%    ax (DEFAULT = gca): Axes handle where visualization is placed.
%    titleStr (DEFAULT = 'error <= [maxerror]'): String to put in plot title.
%       The default setting prints the actual maximum error across all the
%       individual statistics of interest.
%    patchDetail (DEFAULT = 100): Level of detail, as measured by the number of
%       vertices for circle/ellipsoid.
%    SETCAMERA (APPlIES TO 3d ONLY, DEFAULT = true): Logical flag which
%       determines if a standard camera angle is applied (via Hypersphere.camera)
% 
% Optional outputs:
%    objHs: object handles of circles or ellipsoids.
%    XY: coordinates of circles/ellipsoids. Used by Hypersphere.show to properly
%       scale error bar.
% 
% IMPORTANT(!) N.B.: if a 3D (sphere) visualization figure is resized or zoomed
%    in matlab, the error line will not resize accordingly and thus will be
%    inaccurate.
% 
% SEE ALSO DRAWCIRCLE, DRAWELLIPSOID, HYPERSPHERE.SHOW, HYPERSPHERE.CAMERA

model  = varargin{1};
dimLow = size(model.centers,2);

%% control variables
patchDetail = 100;
if     dimLow == 2, opacity = 0.2;
elseif dimLow == 3, opacity = 0.5; SETCAMERA = true;
end

%% preparations
for v = 2:nargin
   if     islogical(varargin{v}),     SETCAMERA      = varargin{v};
   elseif  ishandle(varargin{v}),     ax             = varargin{v};
   elseif    ischar(varargin{v}),     titleStr       = varargin{v};
   elseif isnumeric(varargin{v})
      if      numel(varargin{v})==4,  ax = figurePanel(varargin{v}); % for old niko code compatibility
      elseif  numel(varargin{v})==1,  patchDetail    = varargin{v};
      elseif  numel(varargin{v})==0,  ax             = varargin{v};
      end
   end
end

if ~exist('ax','var') || isempty(ax), ax = gca; end
axtivate(ax);

if ~exist('titleStr','var') && isprop(model,'error') && ~isnan(sum(model.error))
     maxerror = max(model.error)/max([model.radii model.dists]);
     titleStr = any2str('error <= ',ceil(maxerror*100),'%');
else titleStr = [];
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

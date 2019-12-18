function [h,IM] = scatterDensity(varargin)
% Scatterplot of translucent points, to show density
% 
% Examples:
%          scatterDensity(X     )        % X is an [N x 2] matrix
%          scatterDensity(x,y   )        % x & y are N-vectors
% [h,IM] = scatterDensity(x,y,ax)
% [~,IM] = scatterDensity(x,y,sz_pixels) % precalculate IM with size sz_pixels
%          scatterDensity(x,y,IM,ax)     % plot precalulcated IM
%          scatterDensity(x,y,true)      % force axes to be square
%          scatterDensity(X,true,'-interactive') % interactive mode with sq ax
% 
% sz_pixels makes IM such that: fliplr(size(IM)) == sz_pixels
% i.e., sz_pixels = [ylength xlength]
% 
% 2013-02-13 AZ Created
% 2013-03-21 AZ Added sumPointsImage mex capability
% 2013-04-02 AZ IM output; crappy implementation for IM input
%               ax can be sz_pixels (for precalculating IM for later plotting)
% 2013-04-08 AZ Added interactive mode
% 2017-02-21 AZ Added square axis mode

%% Parameters
circrad       =  3; % individual circle radius, in pixels
circsumthresh = 20; % number of overlapping circles to equal black

%% Initialization
sz = size(varargin{1});
if sz(1)>sz(2) && sz(2)==2
   x = varargin{1}(:,1);
   y = varargin{1}(:,2);
   nin = 2;
else
   x = varargin{1};
   y = varargin{2};
   nin = 3;
end

INTERACTIVE = false;
SQUARE      = false;
if strcmpi(varargin{end},'-interactive')
   INTERACTIVE = true;
   varargin{end} = [];
end

IM = [];
ax = gca;

for i = nin:nargin
   if numel(varargin{i}) >= 2
      IM = varargin{i};      
   elseif numel(varargin{i}) == 1
      % ax given, do nothing
      if ishandle(varargin{i})
         ax = varargin{i};
      elseif islogical(varargin{i})
         SQUARE  = varargin{i};
      else
         dotsize = varargin{i};
      end
   end
end

if numel(IM) == 2
   sz_pixels = IM;
   IM = [];
else
   sz_pixels = [];
end


%% Axis setup
nx = numel(x);

% Find most extreme points
xn =   min(x);
xx =   max(x);
yn =   min(y);
yx =   max(y);

set(gcf,'CurrentAxes',ax);
set( ax,'Units','Pixels');
if isempty(sz_pixels)
   sz_pixels    = get(ax,'Position');
   sz_pixels    = round(sz_pixels(3:4));
end

% create image matrix with 10% buffer on all sides
if SQUARE
   xn = min([xn yn]);
   yn = min([xn yn]);
   xx = max([xx yx]);
   yx = max([xx yx]);
end
xspan = xx-xn;
yspan = yx-yn;
sz_plotunits = 1.2.*[xspan yspan];

% Generate matrices that contain x and y coords of the midpoint of each pixel
allx  = linspace(xn-0.1*xspan,xx+0.1*xspan,sz_pixels(1)+1);
ally  = linspace(yn-0.1*yspan,yx+0.1*yspan,sz_pixels(2)+1);
allx  = conv(allx,[0.5 0.5],'valid');
ally  = conv(ally,[0.5 0.5],'valid');

%% Calculate image
if isempty(IM)
   % r = [max(abs(x(:))) max(abs(y(:)))]/50;
   % Set small circle patch radius to be smaller of: 1/50th axis or circrad pixels
   if ~exist('dotsize','var') || isempty(dotsize)
      r = sz_plotunits/50;
      if any(r > circrad*(sz_plotunits./sz_pixels(1:2)))
         r = circrad*(sz_plotunits./sz_pixels(1:2));
      end
   else
      r = sz_plotunits*dotsize;
   end

   % Calculate image
   try
      % Mex version (50% faster)
      IM = sumPointsImage(allx,ally,x,y,r);
   catch err
      warning('%s',err.identifier);
      % If mex file not available, do matlab code
      [XI,YI] = meshgrid(allx,ally);

      % IMlog = false([sz_pixels([2 1]) nx]);
      % parfor i = 1:nx
      %    % Elliptical bound based on axis limits, meant for plotting circles
      %    IMlog(:,:,i) = ((XI - x(i)).^2)/r(1)^2 + ((YI - y(i)).^2)/r(2)^2 < 1;
      % end
      % IM = sum(IMlog,3);

      % Low-memory version, not much slower
      IM = zeros(sz_pixels([2 1]),'single');
      for i = 1:nx
         stationarycounter(i,nx);
         % Elliptical bound based on axis limits, meant for plotting circles
         IM = IM + single(((XI - x(i)).^2)/r(1)^2 + ((YI - y(i)).^2)/r(2)^2 < 1);
      end

      % NOTE: one-line matrix version (code deleted) is WAY SLOWER
   end
end

%% Render Plot
try %#ok<TRYNC>
set(gcf,'CurrentAxes',ax,'Tag','ScatterDensity')%,'ResizeFcn',@figResize);
h = imagesc(allx([1 end]),ally([1 end]),-IM);
colormap gray;caxis([-circsumthresh 0]);
end

if INTERACTIVE
   uicontrol('Parent',get(ax,'Parent'),'Style','slider','Value',-circsumthresh,...
             'Min',-nx,'Max',-1,'SliderStep',[10 100]./nx,...
             'Units','normalized','Position', [0.92 0.1 0.05 0.8],...
             'Callback',@setthresh);
   uicontrol('Style','edit','String', num2str(circsumthresh),'Callback',@setthresh,...
             'Units','normalized','Position',[0.905 0.91 0.07 0.05],'Background',[1 1 1]);
   uicontrol('Style','text','String', num2str(nx),...
             'Units','normalized','Position',[0.91 0.03 0.07 0.05],'Background',[1 1 1]);
end

% Setting units to normalized allows for resizing
set(gcf,'Units','normalized','PaperUnits','normalized');set(ax,'Units','normalized');
axis xy
if SQUARE,   axis equal;   end

% suppress output if not requested
if nargout == 0; clear h; end

end

%% Callback functions

% Set threshold slider
function setthresh(src,~)
   % ax = get(findobj(get(src,'Parent'),'tag','ScatterDensity'));
   ax   = findobj(get(src,'Parent'),'type' ,'axes');
   htxt = findobj(get(src,'Parent'),'style','edit');
   
   switch get(src,'style')
      case 'slider'
         newval = get(src,'Value');
      case 'edit'
         newval = -str2double(get(src,'String'));
   end
   set(ax,'Clim',[newval  0]);
   set(htxt,'String',sprintf('%0.0f',-newval));
end


% % Figure resize function
% function figResize(src,evt)
% 	fpos = get(src,'Children');
%    keyboard;
% 	set(botPanel,'Position',...
% 		[1/20 1/20 fpos(3)-.1 fpos(4)*8/35])
% end

% % Bottom panel resize function
% function botPanelResize(src,evt)
% 	bpos = get(botPanel,'Position');
% 	set(plotButton,'Position',...
% 		[bpos(3)*10/120 bpos(4)*2/8 bpos(3)*24/120 2])
% end

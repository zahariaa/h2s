function [h,XY] = drawCircle(center,radius,color,border,alpha,n)
% h=drawCircle(center=[0 0],radius=1,color='k',border=2,alpha=0.2,n=100)
% -  center must be Nx2
% -  accepts multiple inputs and draws multiple circles
% 
% e.g.:
% [h,XY] = drawCircle;                % draws a black circle at [0 0] with
%                                     % radius = 1, outputs figure handle and
%                                     % XY values of circle
% drawCircle([1 2])
% drawCircle([1 2],3)
% drawCircle([1 2],3,'r')
% drawCircle([1 2],3,[1 0 0])         % same as previous
% drawCircle([1 2],3,'r',0)           % opaque circle
% drawCircle([1 2],3,'r',[],500)      % hi-res!
% drawCircle([1 2],3,'r',[],[],0)     % no border
% drawCircle([1 2;3 4;5 6],3:5,'krg') % 3 circles, in black, red, and green
% drawCircle([1 2;3 4;5 6])           % 3 circles, all black, all radius 1, but
%                                     % with specified centers
% 
% All inputs are optional:
%    center (DEFAULT = [0 0]): [n x 2] matrix of center coordinate(s) for
%       n number of circle(s). n should match n in radius/color argument(s), or
%       equal 1, in which case, that center coordinate is used for all circles.
%    radius (DEFAULT = 1): n-vector of radii for n number of circles. n should
%       match  n should match n in radius argument, or equal 1,
%       in which case, that center coordinate is used for all circles.
%    color (DEFAULT = [0 0 0]): [n x 3] matrix or n-vector of strings for the
%       color(s) of the circle(s). n should match n in center/radius argument(s),
%       or equal 1, in which case, that center coordinate is used for all circles.
%    alpha (DEFAULT = 0.2): The number, in [0,1], determining the circle(s)'
%       FaceAlpha/opacity. Same value applied to all circles.
%    n (DEFAULT = 100): Level of patch detail, as measured by the number of
%       vertices for the circle.
%    border (DEFAULT = 2): The number (in pts) determining the circle border line
%       width. Same value applied to all circles.
% 
% 2018-02-22 AZ Created
% 
% SEE ALSO DRAW3DELLIPSOID, SHOWMODEL

%% Initialize defaults
if ~exist('center','var') || isempty(center);   center = [0 0];     end
if ~exist('radius','var') || isempty(radius);   radius = 1;         end
if ~exist('color' ,'var') || isempty(color );   color  = [0 0 0];   end
if ~exist('alpha' ,'var') || isempty(alpha );   alpha  = 0.2;       end
if ~exist('n'     ,'var') || isempty(n     );   n      = 100;       end
if ~exist('border','var') || isempty(border);   border = 2;         end

%% Recurse if given multiple inputs
nr = numel(radius);
nc = numel(center)/2;
nl = numel(color );
if isnumeric(color);  nl    = nl/3;     % color is a [1x3] vector
else                  color = color(:); % color is a string
end
mn = max([nr nc nl]);
if mn > 1
   % Embed inputs in cell; duplicate variables when necessary
   if nc> 1;   center = num2cell(center,2);
   else        center = repmat({ center(:)'},[mn 1]);
   end
   if nr> 1;   radius = num2cell(radius);
   else        radius = repmat({ radius    },[mn 1]);
   end
   if nl> 1;    color = num2cell(color,2);
   else         color = repmat({ color(:)' },[mn 1]);
   end
   if numel(border) < mn;   border = repmat({border},[mn 1]);   end
   if numel(alpha ) < mn;   alpha  = repmat({alpha },[mn 1]);   end
   if numel(n     ) < mn;   n      = repmat({n     },[mn 1]);   end
   % RECURSE!
   [h,XY] = cellfun( @drawCircle, center,radius(:),color,border,alpha,n,...
                    'UniformOutput',false');
   h  = vertcat( h{:});
   XY = vertcat(XY{:});
   return
end

%% Generate plot
if border==0;   opts = {'LineStyle','none','FaceAlpha',alpha,'EdgeColor',color};
else            opts = {'LineWidth',border,'FaceAlpha',alpha,'EdgeColor',color};
end

t = linspace(0,2*pi,n);
XY = radius*[cos(t(:)) sin(t(:))] + repmat(center(:)',[n 1]);
hold on
h = fill(XY(:,1),XY(:,2),color);
set(h,opts{:});
axis equal off

return


%% Debug
figure;drawCircle([1 2;3 4;5 6],3:5,'krg',0);

function h = drawCircle(center,radius,color,border,alpha,n)
% h = drawCircle(center,radius,<color='k'>,<border=2>,<alpha=0.2>,<n=100>)
% -  center must be Nx2
% -  accepts multiple inputs and draws multiple circles
% 
% 2018-02-22 AZ Created

%% Initialize defaults
if ~exist('color' ,'var') || isempty(color );   color  = [0 0 0];   end
if ~exist('border','var') || isempty(border);   border = 2;         end
if ~exist('alpha' ,'var') || isempty(alpha );   alpha  = 0.2;       end
if ~exist('n'     ,'var') || isempty(n     );   n      = 100;       end

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
   else        center = repmat({center(:)'},[mn 1]);
   end
   if nr> 1;   radius = num2cell(radius);
   else        radius = repmat({radius     },[mn 1]);
   end
   if nl> 1;    color = num2cell(color,2);
   else         color = repmat({ color(:)'},[mn 1]);
   end
   if numel(border) < mn;   border = repmat({border},[mn 1]);   end
   if numel(alpha ) < mn;   alpha  = repmat({alpha },[mn 1]);   end
   if numel(n     ) < mn;   n      = repmat({n     },[mn 1]);   end
   % RECURSE!
   h = cellfun( @drawCircle, center,radius(:),color,border,alpha,n);
   return
end

%% Generate plot
if border==0;   opts = {'LineStyle','none','FaceAlpha',alpha};
else            opts = {'LineWidth',border,'FaceAlpha',alpha};
end

t = linspace(0,2*pi,n);
h = fill(radius*cos(t)+center(1),radius*sin(t)+center(2),color);
set(h,opts{:});
axis equal off
hold on

return

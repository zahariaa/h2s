function [h,XY] = drawCircle(center,radius,color,border,alpha,n)
% h=drawCircle(center=[0 0],radius=1,color='k',border=2,alpha=0.2,n=100)
% -  center must be Nx2
% -  accepts multiple inputs and draws multiple circles
% 
% 2018-02-22 AZ Created

%% Initialize defaults
if ~exist('center','var') || isempty(center);   center = [0 0];     end
if ~exist('radius','var') || isempty(radius);   radius = 1;         end
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
   else        radius = repmat({radius    },[mn 1]);
   end
   if nl> 1;    color = num2cell(color,2);
   else         color = repmat({ color(:)'},[mn 1]);
   end
   if numel(border) < mn;   border = repmat({border},[mn 1]);   end
   if numel(alpha ) < mn;   alpha  = repmat({alpha },[mn 1]);   end
   if numel(n     ) < mn;   n      = repmat({n     },[mn 1]);   end
   % RECURSE!
   [h,XY] = cellfun( @drawCircle, center,radius(:),color,border,alpha,n,...
                     'UniformOutput',false');
    h  = vertcat (h{:});   XY = vertcat(XY{:});
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

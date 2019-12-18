function draw3Daxes(varargin)
% draw3Daxes(ax=gca,axlims=axis,origin=axlims([1 3 5]))
% all inputs are optional and can be provided in any order
% 
% 2018-05-17 AZ Created

%% Parse inputs
for v = 1:nargin
   switch numel(varargin{v})
      case 1;   if ishandle(varargin{v}),   ax = varargin{v};   end
      case 3;   origin = varargin{v};
      case 6;   axlims = varargin{v};
   end
end

%% Assign default values
if ~exist('ax',    'var'),     ax = gca;               end
axtivate(ax);
currentlims = [xlim ylim zlim];
currentlims = currentlims + [diff(currentlims(1:2))*0.001 0 ...
              diff(currentlims(3:4))*0.001 0 diff(currentlims(5:6))*0.001 0];
if ~exist('axlims','var'), axlims = currentlims;       end
if ~exist('origin','var'), origin = axlims([1 3 5]);   end

%% Draw axes
plot3(axlims(1:2),origin([2 2]),origin([3 3]),'k-','LineWidth',1)
plot3(origin([1 1]),axlims(3:4),origin([3 3]),'k-','LineWidth',1)
plot3(origin([1 1]),origin([2 2]),axlims(5:6),'k-','LineWidth',1)


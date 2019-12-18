function logAxis(varargin)
% logAxis(xyz='x',base=10,ax=gca)
% inputs can be provided in any order
% 
% 2018-05-04 AZ Created (again?)

%% Preliminaries
for v = 1:nargin
   if     ishandle( varargin{v}),   ax   = varargin{v};
   elseif isnumeric(varargin{v}),   base = varargin{v};
   elseif ischar(   varargin{v}),   xy   = varargin{v};
   end
end

if ~exist('xy'  ,'var') || isempty(xy  ),   xy   = 'x';   end
if ~exist('base','var') || isempty(base),   base =  10;   end
if ~exist('ax'  ,'var') || isempty(ax  ),   ax   = gca;   end

%% Loop over axes
for xy = xy
   tickLabels = getTickLabels(ax,xy);
   % set ticks
   newLabels  = base.^tickLabels;
   setTickLabels(ax,xy,newLabels)
end

return


function matchy(varargin)
% matchy.m: Force y-axis limits to match for a series of axes passed as
% inputs
% 
% Ex.
% matchy
% matchy(ax)
% matchy('XLim')
% matchy('x')
% matchy('x','y')
% matchy({'x','y'})
% matchy(ax,{'x','y'})
% matchy(ax,{'XLim','YLim','ZLim'})
% matchy('c') % matches color axes
% 
% 2013-04-16 AZ Reverse-engineered to MC's code
% 2014-01-13 AZ Added type input to choose which axes to match
%               Now can pass figure handle to match all axes in figure
% 2014-07-16 AZ Allows multiple types
% 2020-05-28 AZ Made input handling way more flexible/intuitive
%               Added support for color axis matching

%% Parse inputs
type = {};
for v = 1:nargin
	if ishandle(varargin{v}), ax = varargin{v};
	else
		if ~iscell(varargin{v}), varargin{v} = {varargin{v}}; end
		for i = 1:numel(varargin{v})
			switch( upper(varargin{v}{i}) )
				case {'X','Y','Z','C'}, type = [type [varargin{v}{i} 'Lim']];
				otherwise 					type = [type  varargin{v}{i}       ];
			end
		end
	end
end

%% Set defaults
if isempty(type), type = {'XLim','YLim'}; end
if ~exist('ax','var') || isempty(ax), ax = gcf; end

if strcmpi(get(ax,'type'),'figure')
   ax = get(ax,'Children');
   ax = ax(strcmpi(get(ax,'Type'),'axes')); % choose only axes-type childrens
end


%% Initialize
na = numel(ax);
yy = zeros(na,2);

%% Main loop
for i = 1:numel(type)
	for a = 1:na
	   yy(a,:) = get(ax(a),type{i});
	end

	newy = [min(yy(:,1)) max(yy(:,2))];

	for a = 1:na
	   set(ax(a),type{i},newy);
	end
end

function matchy(ax,type)
% matchy.m: Force y-axis limits to match for a series of axes passed as
% inputs
% 
% Ex.
% matchy
% matchy(ax)
% matchy([], 'XLim')
% matchy(ax, 'XLim')
% matchy(ax,{'XLim','YLim'})
% matchy(ax,{'XLim','YLim','ZLim'})
% 
% 2013-04-16 AZ Reverse-engineered to MC's code
% 2014-01-13 AZ Added type input to choose which axes to match
%               Now can pass figure handle to match all axes in figure
% 2014-07-16 AZ Allows multiple types

if ~exist('ax','var') || isempty(ax)
   ax = gcf;
end

if strcmpi(get(ax,'type'),'figure')
   ax = get(ax,'Children');
   ax = ax(strcmpi(get(ax,'Type'),'axes')); % choose only axes-type childrens
end

if ~exist('type','var') || isempty(type)
   cellfun(@(x) matchy(ax,x),{'XLim','YLim'})
   return
elseif iscell(type)
   cellfun(@(x) matchy(ax,x),type)
   return
end

na = numel(ax);

yy = zeros(na,2);
for a = 1:na
   yy(a,:) = get(ax(a),type);
end

newy = [min(yy(:,1)) max(yy(:,2))];

for a = 1:na
   set(ax(a),type,newy);
end
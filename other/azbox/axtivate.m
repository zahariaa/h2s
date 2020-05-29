function axtivate(a,varargin)
% axtivate(a,f): activates axis in fh.a(f).h(a)
% assumes f is the last one available if not specified.
% pulls fh in workspace from calling function.
% 
% SEE ALSO: NEWFIGURE, BATCHPLOTREFINE, FIGSETUP
% 
% 2016-08-04 AZ Created

vstart = 1;

if ishandle(a(1)) && get(a(1),'Parent')~=0
   for ax = a(:)'
      fig = get(ax,'Parent');
   end
elseif evalin('caller','exist(''fh'',''var'')')
   fh = evalin('caller','fh');
   if isempty(varargin) || isempty(varargin{1}) || ~isnumeric(varargin{1})
   	  f = numel(fh.f);
   else f = varargin{1}; vstart = 2;
   end
   fig = fh.f(f);
   ax  = fh.a(f).h(a);
else
	fig = gcf;
	ax  = gca;
end

set(0,'CurrentFigure',fig);
set(fig,'CurrentAxes',ax);
hold on;

for v = vstart:2:numel(varargin)
	ax.(varargin{v}) = varargin{v+1};
end

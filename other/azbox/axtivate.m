function varargout = axtivate(a,varargin)
% axtivate(a,f): activates axis in fh.a(f).h(a)
% assumes f is the last one available if not specified.
% pulls fh in workspace from calling function.
% 
% ex.:
% axtivate(2)     % activates f=gcf;       f:subplot(2)   or fh.a(f).h(2)
% axtivate(2,3)   % activates f=figure(3); f:subplot(2)   or fh.a(3).h(2)
% axtivate([2 4]) % activates f=gcf;       f:subplot(2,4) or fh.a(f).h(a), a=subplot(2,4)
% axtivate(2,'XScale','log') % activates 2 and makes x-axis on a log scale
% arrayfun(@(a) axtivate(a,'XScale','log'),1:5) % does above for subplots 1:5
% 
% SEE ALSO: NEWFIGURE, BATCHPLOTREFINE, FIGSETUP
% 
% 2016-08-04 AZ Created

vstart = 1;

if ishandle(a(1)) && get(a(1),'Parent')~=0
   for ax = a(:)'
      fig = get(ax,'Parent');
   end
elseif isnumeric(a)
   % Determine parent figure
   if ~isempty(varargin) && ~isempty(varargin{1}) && ~ischar(varargin{1})
      if isnumeric(varargin{1}) || ishandle(varargin{1}), fig = varargin{1};
      else   error('uhoh')
      end
      vstart = 2;
   else fig = gcf;
   end

   % Address proper child axis
   kids = flipud(get(fig,'Children'));
   switch numel(a)
      case 1; a = [1 a];
      case 2  % address subplots with pair of indices
         % determine subplot arrangement
         pos = vertcat(kids.Position);
         nx  = numel(unique(pos(:,1)));
         ny  = numel(unique(pos(:,2)));
         kids= reshape(kids,[nx ny]);
   end
   ax  = kids(a(2),a(1));
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

% Apply property (name,value) pairs to activated axis
for v = vstart:2:numel(varargin)
	ax.(varargin{v}) = varargin{v+1};
end

if nargout, varargout = {ax}; end

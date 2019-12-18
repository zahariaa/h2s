function axtivate(a,f)
% axtivate(a,f): activates axis in fh.a(f).h(a)
% assumes f is the last one available if not specified.
% pulls fh in workspace from calling function.
% 
% SEE ALSO: NEWFIGURE, BATCHPLOTREFINE, FIGSETUP
% 
% 2016-08-04 AZ Created

if ishandle(a(1)) && get(a(1),'Parent')~=0
   for ax = a(:)'
      fig = get(ax,'Parent');
      set(0,'CurrentFigure',fig);
      set(fig,'CurrentAxes',ax);
   end
elseif evalin('caller','exist(''fh'',''var'')')
   fh = evalin('caller','fh');
   if ~exist('f','var') || isempty(f);   f = numel(fh.f);   end

   set(fh.f(f),'CurrentAxes',fh.a(f).h(a));
end
hold on;


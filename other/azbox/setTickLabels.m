function setTickLabels(ax,xy,tickLabels)
% setTickLabels(ax=gca,xy,tickLabels)
%
% 2018-05 AZ Created

if ~exist('ax','var') || isempty(ax),   ax = gca;   end

%% Recurse
if numel(ax) > 1
   arrayfun(@(AX) setTickLabels(AX,xy,tickLabels),ax);   return;
end

%% Detect if axesSeparate was run
sepAxEl    = separateAxisElements(ax,['phyplot_' xy 'ticklabel']);

%% write ticks
if isempty(sepAxEl),   set(ax,[xy 'TickLabel'],num2str(tickLabels));
elseif HG2             arrayfun(@num2str,tickLabels,'UniformOutput',false)
else                   arrayfun(@(a,s) set(a,'String',s),sepAxEl,tickLabels);
end

return


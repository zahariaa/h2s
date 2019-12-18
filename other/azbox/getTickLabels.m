function tickLabels = getTickLabels(ax,xy)
% tickLabels = getTickLabels(ax=gca,xy='xy')
% 
% 2018-05 AZ Created

if ~exist('ax','var') || isempty(ax),   ax = gca;   end
if ~exist('xy','var') || isempty(xy),   xy = 'xy';  end

%% Recurse
if     numel(xy) > 1 && numel(ax)==1
   tickLabels = arrayfun(@(   XY) getTickLabels(ax,XY),xy);   return;
elseif numel(ax) > 1 && numel(xy)==1
   tickLabels = arrayfun(@(AX   ) getTickLabels(AX,xy),ax);   return;
elseif numel(ax) > 1 && numel(xy) > 1
   tickLabels = arrayfun(@(AX,XY) getTickLabels(AX,XY),ax,xy);return;
end

%% Detect if axesSeparate was run
sepAxEl    = separateAxisElements(ax,['phyplot_' xy 'ticklabel']);

%% get ticks
if ~isempty(sepAxEl),   tickLabels = str2double(get(sepAxEl,'String'));
else                    tickLabels = get(ax,[xy 'TickLabel']);
   if iscell(tickLabels),  tickLabels = cellfun(@str2num,tickLabels); HG2=true;
   else                    tickLabels = str2num(tickLabels);          HG2=false;
   end
end

return



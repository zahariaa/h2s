function elements = separateAxisElements(ax,searchstr)
% separateAxisElements(ax,searchstr)
% Run this after running axesSeparate (or batchPlotRefine).
% 
% Ex. Leave only left two y-axes in plot 1 & n+1 in a 2 by n figure:
% delete(separateAxisElements(ax([2:n n+2:2*n]),'phyplot_y'))
% 
% 2017-12-20 AZ Created
% 
% SEE ALSO BATCHPLOTREFINE, AXESSEPARATE

% Recurse
if numel(ax)>1
   elements = cell2mat_concat(...
                   arrayfun(@(x) separateAxisElements(x,searchstr),ax,...
                   'UniformOutput',false) );
   return
end

kids = get(ax,'Children');
if isempty(kids);   elements = [];   return;   end

elements = kids(~cellfun(@isempty,strfind(get(kids,'Tag'),searchstr)));

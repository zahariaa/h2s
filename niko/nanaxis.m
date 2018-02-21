function nanaxis(limits)

% sets the axis limits like the built-in function axis, but allows nans,
% leaving the corresponding axis limits as they are. moreover, trailing
% nans can be omitted.

%% find old axis limits
xl=get(gca,'XLim');
yl=get(gca,'YLim');
oldLims=[xl yl];

%% define and set new axis limits
newLims=nan(1,4);
newLims(1:numel(limits))=limits;
newLims(isnan(newLims))=oldLims(isnan(newLims));

axis(newLims);




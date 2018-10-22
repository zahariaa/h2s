function varargout = figurePanel(figPanelSpec)

% create or select figure panel
figure(figPanelSpec(1)); set(figPanelSpec(1),'Color','w');
if numel(figPanelSpec)==4
    ax = subplot(figPanelSpec(2), figPanelSpec(3), figPanelSpec(4));
end
cla; hold on;

if nargout, varargout{1} = ax; end


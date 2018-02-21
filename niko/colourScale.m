function cols=colourScale(anchorCols,nCols,monitor)
% linearly interpolates between a set of given 'anchor' colours to give
% nCols and displays them if monitor is set 



%% preparations
if ~exist('monitor','var'), monitor=false; end


%% define color scale
nAnchors=size(anchorCols,1);
cols = interp1((1:nAnchors)',anchorCols,linspace(1,nAnchors,nCols));


%% visualise
if monitor
    figure(123); clf;
    imagesc(reshape(cols,[nCols 1 3]));
end



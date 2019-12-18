function subplotResize(hands,offset,margin)
% subplotResize: collects subplots (even if not made by subplot), and
% maximizes their plot areas. assumes a uniform grid of subplots.
% 
% offset adds space to bottom & left, margin adds space on all 4 sides
% 
% Ex.
% To run this manually, execute the following before running batchPlotRefine:
% subplotResize(fh.a(end).h,0.03); fh.a(end).resize=[];
% 
% 2016-04-25 AZ Created
% 
% SEE ALSO AXESSEPARATE, BATCHPLOTREFINE, FIGSETUP, NEWFIGURE, PRINTFIG

%% Initialization
if    ~exist('margin','var') ||  isempty(margin);      margin = [0 0];
elseif numel(margin)==1,                               margin = margin*[1 1];
end
if    ~exist('offset','var') ||  isempty(offset);      offset = [0 0];    
elseif numel(offset)==1,                               offset = offset*[1 1];
end
if    ~exist('hands' ,'var') ||  isempty(hands ),      hands  = gcf;      end
% get handle type
htype    = get(hands,'Type');      if ischar(htype);   htype  = {htype};  end
% collect axes from figure
if strcmpi(htype{1},'figure');
   hands = get(hands,'Children');
   hands = findobj(hands,'type','axes'); % necessary in HG2
end

%% Figure out axes positioning by ordering axis handles
% Reorder by x position
sizes = get(hands,'Position');     if iscell(sizes);   sizes = cell2mat(sizes);   end
nx    = numel(unique(round(sizes(:,1)*1000)/1000));
ny    = numel(unique(round(sizes(:,2)*1000)/1000));
hands = hands(sortix(sizes(:,1)));

hands = reshape(hands,[ny nx]);

% Reorder by y position within same x position
sizes = get(hands,'Position');     if iscell(sizes);   sizes = cell2mat(sizes);   end
for x = unique(sizes(:,1))'
   ix        = sizes(:,1)==x;
   hands(ix) = indexm(hands(ix),sortix(sizes(ix,2)));
end
[ny,nx] = size(hands);

%% Apply axis sizes
for x = 1:nx
   for y = 1:ny
      set(hands(y,x),'Position',[[(x-1)/nx (y-1)/ny]+offset [1/nx 1/ny]-(margin+offset)])
   end
end


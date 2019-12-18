function plotTuning(p,tun,subset,colors,CENTER,PATCH,ax)
% PLOTTUNING: plots tuning curves from basic characterizations
% Requires either p struct, or tun struct with the following fields:
%    bins, respmean, respstd, spontmean, spontstd
%
% Ex.,
% plotTuning(p,tun,subset,colors,CENTER,PATCH,ax)
% plotTuning(p)
% plotTuning(p,tun)
% plotTuning(p,tun,1)                   % plots the 1st subset from p only
%                                       % (good for plaid120SS)
% plotTuning(tun)
% plotTuning(tun,[],[],colors)
% plotTuning(p,[],1)
% plotTuning(p,tun,3:4,colors)          % plots 3rd & 4th subsets from p
%                                       % with specified colors
% 
% See also PLOTERRORPATCH, EXTRACTTUNING, LOADP
% 
% 2011-08-05 AZ Created
% 2012-01-30 AZ Allows for multiple curve plotting (i.e. plaid120)
% 2012-12-11 AZ plaid120 & ori16 axis (-180,180] pref dir fix

%% Parse inputs
if ~exist('p','var')
   p = [];
elseif isfield(p,'bins') && ( ~exist('tun','var') || isempty(tun) )
   tun = p;
end

if ~exist('tun'   ,'var') && isempty(tun   );   tun = extractTuning(p,[],[],true);end
if ~exist('CENTER','var') || isempty(CENTER);   CENTER = false;        end
if ~exist('PATCH' ,'var') || isempty(PATCH );   PATCH  = [false true]; end
if ~exist('subset','var') || isempty(subset)
   if numel(tun)==1;   subset = 1:size(tun.respmean,2);
   else                subset = 1:numel(tun);
   end
end

nb = numel(subset);
% Color order
if ~exist('colors','var') || isempty(colors)
   if verLessThan('matlab','8.4')
      colors = {[0 0 0],[1 0 0],[0 1 0],[0 0 1]};
   else
      colors = {lines(nb)};
   end
elseif isnumeric(colors)
   colors = {colors};
elseif isstr(colors)
   colors = num2cell(colors);
end

%% Recurse if multiple axes given
if     ~exist('ax','var') || isempty(ax), ax = gca;
elseif  exist('ax','var') && numel(ax)>1
  if numel(tun) == numel(ax)
     arrayfun(@(x,t) plotTuning(p,t,  1,colors,CENTER,PATCH,x),ax(:),tun(:));
  else
     arrayfun(@(x,s) plotTuning(p,tun,s,colors,CENTER,PATCH,x),ax(:),subset(:));
  end
  return
end

%% Convert to preferred direction, if desired.
if CENTER
orPrefGuess = [];
try
   orPrefGuess = str2double(regexp(p.specialDescription,'([\-]{0,1}[0-9]{1,3}\.[0-9]{0,1})','match'));
end
if isempty(orPrefGuess)
   orPrefGuess = tun.bins(maxix(tun.respmean(:,1)),1);
end
tun.bins = tun.bins - orPrefGuess;
tun      = tuningWrapSort(tun);
end

% Deal.
respmean  = tun.respmean;
respstd   = tun.respstd;
spontmean = tun.spontmean;
spontstd  = tun.spontstd;
bins      = tun.bins;
% Expand if necessary
if size(bins,2) < size(respmean,2)
   bins   = bins(:,[1 1 2 2]);
end

% make sure there are enough colors
if numel(colors)==1 && nb>1
   if numel(colors{1})==3
      colors = mat2cell(makeColorSaturationSeries(colors{1},nb),ones(1,nb));
   else
      colors = mat2cell(                          colors{1}    ,ones(1,nb));
   end
end

%% Manual log scaling to avoid openGL shit
bins = logify(bins);

%% PLOT
% Spontaneous rate
if numel(PATCH)==2 && PATCH(2)
   for i = vectify(subset(end:-1:1))'
      nc = size(bins(:,1));
      plotErrorPatch(bins(:,i),spontmean*ones(nc),spontstd*ones(nc),[0.95 0.95 0.95],PATCH(2),ax);
   end
else
   axtivate(ax);  hold on;
   for i = subset(end:-1:1)'
      plot(bins(:,i),repmat(spontmean-spontstd,size(bins(:,i))),...
                          '-','Color',colors{1},'Tag','phyaxwidth');
      plot(bins(:,i),repmat(spontmean+spontstd,size(bins(:,i))),...
                          '-','Color',colors{1},'Tag','phyaxwidth');
   end
end

% Response rate
c = 1;
for i = subset(:)'
   plotErrorPatch(bins(:,i),respmean(:,i),respstd(:,i),colors{c},PATCH(1),ax);
   c = c+1;
end

%% LABEL
ylabel('Response (spikes/s)');

switch p.protocol
	case {'plaid120','ori16','ori16c','plaid120_B','ori16_B','ori16c_B','plaid120SS','plaid120SS_B'}
      if     min(bins)==0
         xlabel('Direction (deg)');
         set(ax,'XLim',[0 360],...
                 'XTick',0:90:360,'XTickLabel',[' 0 ' ;' 90' ;'180';'270';'360'])
      elseif min(bins)<0
         xlabel('Preferred direction (deg)');
         set(ax,'XLim',[-180 180],...
                 'XTick',-180:90:180,'XTickLabel',['-180' ;' -90' ;'  0 ';' 90 ';' 180'])
      end
	case {'sf11','sf11c','sf11_B','sf11c_B'}
		xlabel('Spatial frequency (c/deg)');
      set(ax,'XTick',-1:1     ,'XTickLabel',['0.1' ;'1  ' ;'10 '],'XLim',[-1 1])
	case {'tf10','tf10low','tf10high','tf10_B','tf10low_B','tf10high_B'}
		xlabel('Drift rate (Hz)');
      set(ax,'XTick', 0:2     ,'XTickLabel',['1  ';'10 ';'100']  ,'XLim',[0 2])
   case {'rvc10','rvc10_B','rvcSur','rvcOrtho'}
      xlabel('log10 Contrast');
   case {'rfsize9c','rfsize9c_B'}
      xlabel('Size');
	otherwise
		% do nothing
end

% Force y-axis to be positive
ylim(max(0,ylim))

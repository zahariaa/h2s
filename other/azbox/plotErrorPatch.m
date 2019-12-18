function varargout = plotErrorPatch(varargin)
% Plots error bars (std) as a patch behind the mean
% 
% Ex.:
% plotErrorPatch([],y)
% plotErrorPatch(x,y)
% plotErrorPatch(x,ymean,ystd)
% plotErrorPatch(x,ymean,ystd,[1 0 0])
% plotErrorPatch(x,y,[1 0 0])
% plotErrorPatch([],y,ax)
% ax = plotErrorPatch([],y)
% plotErrorPatch(x,y,false) % plots error (line) bars & dots instead of patch
% 
% 2014-01-09 AZ Created

% defaults
PATCH = true;
SEM   = false;

%% Input parsing
if isstruct(varargin{1}) && isfield(varargin{1},'bins')
   xmean = varargin{1}.bins;
   ymean = varargin{1}.respmean;
   ystd  = varargin{1}.respstd;
else
   nums   = cellfun(@isnumeric    ,varargin);
   mats   = cellfun(@ismatrix     ,varargin);
   n      = cellfun(@numel        ,varargin);
   flags  = cellfun(@islogical    ,varargin);
   chars  = cellfun(@ischar       ,varargin);
   mats2d = cellfun(@(x) (size(x,1)>1 && size(x,2)>1),varargin);
   
   if any(flags)
      PATCH = varargin{flags};
   end
   for c = find(chars(:))'
   switch lower(varargin{c})
      case 'sem';   semOrStd = @sem;
   end
   end
   n(chars) = 0;
   if ~exist('semOrStd','var'),   semOrStd = @nanstd;   end
   
   % Make sure mats are 2d and nums are 1d
   mats((~mats2d &  mats) |  flags |  chars) = false;
   nums((~mats2d & ~nums) & ~flags & ~chars) = true;
   
   % Find axis handles
   ax = find(n==1).*cellfun(@ishandle,varargin(n==1));
   ax = ax(find(ax));
   if ~isempty(ax);   nums(ax) = false;   ax = varargin{ax};   end
   
   if any(n(2:end)==3)
      color = varargin{[false n(2:end)==3]};
      mats([false n(2:end)==3]) = false;
      nums([false n(2:end)==3]) = false;
   else
      color = [0 0 1];
   end
   
   nums(mats) = false;
   
   % Fill in x.
   if     mats(1),              xmean = nanmean( varargin{1});
                                xstd  = semOrStd(varargin{1});
   elseif nums(1),              xmean =          varargin{1} ;
   end
   % Fill in y.
   if     mats(2),              ymean = nanmean( varargin{2});
                                ystd  = semOrStd(varargin{2});
   elseif nums(2),              ymean =          varargin{2} ;
   elseif isempty(varargin{2}), ymean = 1:size(xmean,2);
   end
   if numel(nums)>2 && nums(3), ystd  =          varargin{3} ;  end
   if ~exist('xmean','var') || isempty(xmean), xmean = 1:size(ymean,2); end
   
   if     ~exist('ystd','var') && exist('xstd','var'), ystd = 0*ymean;
   elseif ~exist('xstd','var') && exist('ystd','var'), xstd = 0*xmean;
   end
end
% If ymean is actually all of y, and not averaged, take the average
if numel(ymean)~=numel(xmean);  ymean = nanmean(ymean);   end
% Do we need to flip the y dimension?
dim2flip = find(size(ymean)==max(size(ymean)),1);

%% Cut NaN's
xmean = xmean(~isnan(ymean));
if numel(ystd)>1
   ystd  =  ystd(~isnan(ymean));
end
xstd  =  xstd(~isnan(ystd));
ymean = ymean(~isnan(ymean));

%% Plot stuff
if isempty(ax); ax = gca; end
set(0,'CurrentFigure',get(ax,'Parent'));
set(get(ax,'Parent'),'CurrentAxes',ax);
hold on

if PATCH
patch(cat(dim2flip,xmean+xstd,flip(xmean-xstd,dim2flip)),...
      cat(dim2flip,ymean+ystd,flip(ymean-ystd,dim2flip)),...
      color,'EdgeColor','none','FaceAlpha',1/3);
plot(xmean,ymean,'Color',color,'LineWidth',2);
else
   for i = 1:numel(xmean)
      plot([xmean(i)+xstd(i) xmean(i)-xstd(i)], [ymean(i)+ystd(i) ymean(i)-ystd(i)],'-','Color',color,'Tag','phyaxwidth');
      plot( xmean(i),ymean(i),'o','MarkerFaceColor',color,'LineWidth',2,'Color','w');
   end
end

%% Outputs
if nargout > 0
   varargout{1} = ax;
end

function varargout = axesSeparate(varargin)
% SEPARATEAXES: Replaces MATLAB axes with phyplot style axes
% 
% 2014-11-12 AZ Forked from PHYPLOT by RDK

for i = 1:nargin
   if      ishandle(varargin{i});      ax = varargin{i};
   elseif  isstruct(varargin{i});      fp = varargin{i};
   end
end

if ~exist('ax','var') || isempty(ax);  ax =         gca;  end
if ~exist('fp','var') || isempty(fp);          figsetup;  end
if ~isfield(fp,'signif') || isempty(fp.signif) || fp.signif < 0
   fp.signif = 0;
end
if ~isfield(fp,'AXSCALE');   fp.AXSCALE = [];   end

%% Defaults
linewidth.x     = 1;
linewidth.y     = 1;
fontsize   = fp.txt{8};
fontname   = fp.txt{2};
fontangle  = fp.txt{4};         % normal, italic, oblique
fontweight = fp.txt{6};         % normal, bold

% Determine spacing and tick lengths
if     strcmpi(fp.pap.u,'inches');        textsize     =      fontsize/72;
elseif strcmpi(fp.pap.u,'centimeters');   textsize     = 2.54*fontsize/72;
else                                      fontsize     =      fontsize/1.25;
                                          textsize     =      fontsize/144; %fp.pap.u=='points'
end

%% Force focus
set(0,'CurrentFigure',get(ax,'Parent'));
set(get(ax,'Parent'),'CurrentAxes',ax);

%% Detect previous run and reset axes bounds if necessary
oldrun = separateAxisElements(ax,'phyplot');
if ~isempty(oldrun),   delete(oldrun); axis tight;   end

%% Extract tick and line parameters from plot
axislimits   = axis;
axisranges   = abs(diff(axislimits));
axisranges   = axisranges(1:2:end);

for xy = 'xy'
   j = find('xy'==xy);
   ticklocations.(xy) = get(ax,[xy 'Tick']);
   ticklabels.(   xy) = getTickLabels(ax,xy);
   
   if numel(ticklocations.(xy))~=size(ticklabels.(xy),1)
      ticklocations.(xy) = axislimits((j*2-1):j*2);
   end
   
   % Detect log or linear spacing based on ticklabels
   base.(xy) = basedetect(ticklabels.(xy));

   if isempty(ticklabels.(xy)) && isempty(ticklocations.(xy))
      linewidth.(    xy) = 0;
   elseif base.(xy)==1 && size(ticklabels.(xy),1)>1
      % Cut out ticks outside of axis limits
      ticklocations.(xy) = ticklocations.(xy)(ticklocations.(xy)>=axislimits(2*(j-1)+1)  );
      ticklocations.(xy) = ticklocations.(xy)(ticklocations.(xy)<=axislimits(2*(j-1)+2)  );
   end
   if strcmpi(fp.AXSCALE,'norm')
      ticklabels.(   xy) =    normcdf(ticklocations.(xy)');
   elseif base.(xy)~=1
      ticklabels.(   xy) = base.(xy).^ticklocations.(xy)' ;
   else
      ticklabels.(   xy) =            ticklocations.(xy)' ;
   end
end

%% Convert spacings to axis coordinates, set figure spacing parameters
textsize     = textsize*axisranges; % cumulative padding, in multiples of textsize
dataaxispad  = textsize/2;          % 0-0.5x      padding
ticklength   = textsize/3;          % 0.5-1x      tick
ticklabelpad = textsize*4/4;        % 1.25-2.25x  ticklabel
% axis label is 1-totalpadding      % 2.5-3.5x    axislabel
totalpadding = textsize*2.5;        % 3.5x        subplot border

%% Resize axis to axis+padding
newaranges   = axislimits(:)' + [-1 1/3 -1 1/3].*totalpadding([1 1 2 2]);
axis(newaranges);
axis off;

axh = [];
%% DRAW X & Y AXES
% stuff that differs between x & y
ai        = struct('x',[3 2 1 2],'y',[1 1 3 4]); % indices
textalign = struct('h',struct('x','center','y','right'),'v',struct('x','top','y','middle'));
for xy = 'xy'
if linewidth.(xy)
   % main line: define x & y spans in temp variable l
   l.(xy)                       = ticklocations.(xy)([1 end]);
   l.(indexm('xy',[],'xy'~=xy)) = [axislimits(ai.(xy)(1)) axislimits(ai.(xy)(1))]-dataaxispad(ai.(xy)(2));
   axh = [axh line(l.x,l.y,'color','k','linewidth',linewidth.(xy),'Tag',['phyplot_' xy 'axis'])];
   % extra lines on either end: define spans again
   l.(xy) = [axislimits(ai.(xy)(3)) ticklocations.(xy)(1) NaN ticklocations.(xy)(end) axislimits(ai.(xy)(4))];
   l.(indexm('xy',[],'xy'~=xy)) = axislimits(ai.(xy)(1)*ones(1,5))-dataaxispad(ai.(xy)(2));
   axh = [axh line(l.x,l.y,'color','k','linewidth',linewidth.(xy),'Tag',['phyplot_' xy 'axis_extra'])];
   for i=1:min(length(ticklocations.(xy)),size(ticklabels.(xy),1))
      if (fontsize>0)
         if ~isfield(ticklabels,xy)
            if ticklabels.(xy)(i) == round(ticklabels.(xy)(i))
               tempstring = sprintf('%.0f',ticklabels.(xy)(i));
            else
               roundoff = max(0,fp.signif-floor(log10(max(abs(axislimits(ai.(xy)(3:4)))))));
               tempstring = sprintf(sprintf('%%.%uf',roundoff),ticklabels.(xy)(i));
            end
         else
            roundoff = max(0,fp.signif-floor(log10(max(abs(axislimits(ai.(xy)(3:4)))))));
            tempstring = sprintf(sprintf('%%.%uf',roundoff),ticklabels.(xy)(i));
         end
         t.(xy)                       = ticklocations.(xy)(i)*(1+0.005*any(xy=='x'));
         t.(indexm('xy',[],'xy'~=xy)) = axislimits(ai.(xy)(1))-ticklabelpad(ai.(xy)(2));
         text(t.x,t.y,tempstring,...
            'horizontalalignment',textalign.h.(xy),'fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
              'verticalalignment',textalign.v.(xy),'fontweight',fontweight,'Tag',['phyplot_' xy 'ticklabel']);
         
      end
      l.(xy)                       = ticklocations.(xy)([i i]);
      l.(indexm('xy',[],'xy'~=xy)) = [axislimits(ai.(xy)(1)) axislimits(ai.(xy)(1))-ticklength(ai.(xy)(2))]-dataaxispad(ai.(xy)(2));
      axh = [axh line(l.x,l.y,'color','k','linewidth',linewidth.(xy),'Tag',['phyplot_' xy 'tick'])]; %#ok<*AGROW>
   end
   
   % If base-10, plot minor ticks
   if base.(xy)~=1 && fp.MINORTICKS
      ticklabels.(xy) = unique([base.(xy)^axislimits(ai.(xy)(3)) ticklabels.(xy)(:)' base.(xy)^axislimits(ai.(xy)(4))]);
      MINORTICKS.(xy) = [];
      for i = 1:numel(ticklabels.(xy))-1
         interval     = base.(xy)^floor(log(ticklabels.(xy)(i))/log(base.(xy)));
         MINORTICKS.(xy) = [MINORTICKS.(xy) ticklabels.(xy)(i):interval:ticklabels.(xy)(i+1)];
      end
      MINORTICKS.(xy) = log(unique(setdiff(MINORTICKS.(xy),ticklabels.(xy))))/log(base.(xy));
      for i = 1:numel(MINORTICKS.(xy))
         l.(xy) = MINORTICKS.(xy)([i i]);
         l.(indexm('xy',[],'xy'~=xy)) = [axislimits(ai.(xy)(1)) axislimits(ai.(xy)(1))-ticklength(ai.(xy)(2))/2]-dataaxispad(ai.(xy)(2));
      axh = [axh line(l.x,l.y,'color','k','linewidth',linewidth.(xy),'Tag',['phyplot_' xy 'axis_minor'])]; %#ok<*AGROW>
      end
   end
end
end % end for xy

%% TAKE CARE OF LABELS
Xpos = get(get(ax,'Xlabel'),'Position');   Xpos(1:2) = [mean(axislimits(1:2)) newaranges(3)];
Ypos = get(get(ax,'Ylabel'),'Position');   Ypos(1:2) = [newaranges(1) mean(axislimits(3:4))];
Tpos = get(get(ax,'Title' ),'Position');   Tpos(1:2) = [mean(axislimits(1:2)) newaranges(4)];

set(get(ax,'Xlabel'),'Visible','on',     'Position',Xpos,    'FontAngle', fontangle,...
                     'FontName',fontname,'FontSize',fontsize,'FontWeight',fontweight,...
                     'VerticalAlignment','bottom','HorizontalAlignment',  'center');
set(get(ax,'Ylabel'),'Visible','on',     'Position',Ypos,    'FontAngle' ,fontangle,...
                     'FontName',fontname,'FontSize',fontsize,'FontWeight',fontweight,...
                     'VerticalAlignment','top',   'HorizontalAlignment',  'center');
set(get(ax,'Title'), 'Visible','on',     'Position',Tpos,    'FontAngle', fontangle,...
                     'FontName',fontname,'FontSize',fontsize,'FontWeight',fontweight,...
                     'VerticalAlignment','top',   'HorizontalAlignment',  'center');

if nargout; varargout{1} = axh; end


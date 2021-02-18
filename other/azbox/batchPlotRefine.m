function varargout = batchPlotRefine(fh,fp,includeFigs)
% batch apply text formatting.
% fh = figure handle struct; put all figure handles in fh.f
% fp = figure paramters struct
% run figsetup.m to initialize fp
%
% can run with both, only one, or no inputs:
% batchPlotRefine(fh,fp);
% batchPlotRefine(fh   );
% batchPlotRefine([],fp);
% batchPlotRefine;
% 
% 2011-10-06 AZ Created
% 2011-10-12 AZ Updated to run without inputs
% 2011-10-18 AZ Fixed so axis tick labels also resized
% 2012-01-19 AZ Added includeFigs argument
%               Fixed case where fh.a(f) doesn't exist
% 2012-02-15 AZ Added fp.KEEPLINEWIDTH option
%               Added tick label fix
% 2012-03-14 AZ Fixed bug: force apply changes to fp.pap to fp.pap.opts
% 2013-05-15 AZ Disabled warnings
% 2013-06-23 AZ Bugfix for HG2 compatibility

% HG2 = feature('usehg2');

%% Initialization
if ~exist('fh','var') || isempty(fh)
	fh.f = get(0,'Children');
	for f = 1:numel(fh.f)
		kids = get(fh.f(f),'Children');
		for a = 1:numel(kids)
			fh.a(f).h(a) = kids(a);
		end
	end
end

if ~exist('fp','var') || isempty(fp)
	figsetup;
else
   % NOTE: this will fail for pdf printing
	% Re-apply any changes made to fp
	fp.pap.opts = {'PaperUnits',fp.pap.u,'PaperSize',fp.pap.sz,...
               'PaperPosition',fp.pap.pos,'Position',fp.scr.pos,...
               'PaperPositionMode', 'manual','Renderer',fp.pap.opts{14}};
   if isfield(fp,'DESTINATION') && ~iscell(fp.DESTINATION); fp.DESTINATION = {fp.DESTINATION}; end
end
fp.txt{end} = fp.txtsz;

% Cut out closed figures
fh.a = fh.a(ishandle(fh.f));
fh.f = fh.f(ishandle(fh.f));

if ~exist('includeFigs','var') || isempty(includeFigs)
   includeFigs = 1:numel(fh.f);
   % Cut out completely excluded figs
   if isfield(fh.a,'ex')
      hands       = cellfun(@ishandle,{fh.a.ex},'UniformOutput',false);
      hands       = cellfun(@all,hands) .* ~cellfun(@isempty,hands);
      %toExclude   = cellfun(@all,{fh.a.ex}) .* ~cellfun(@isempty,{fh.a.ex}) .* ~hands;
      %includeFigs = includeFigs(~toExclude);
   end
end

%% cycle through all axes, find & apply

% strcmpi({fh.f(1).Children.Children.Type},'line') %HG2

alltext = findall(fh.f(includeFigs),'Type','text');
allline = findall(fh.f(includeFigs),'Type','line');
allclin = findall(fh.f(includeFigs),'Type','hggroup'); % contour lines
allpaxw = findall(fh.f(includeFigs),'Tag' ,'phyaxwidth');
allline = setdiff(allline,allpaxw);

for f = includeFigs%1:numel(fh.f(includeFigs))
   % Set paper options for figure
   nplotpos = numel(fp.scr.pos/4);
   if nplotpos > 1 && isfield(fp,'DESTINATION') && strcmpi(fp.DESTINATION{1},'onlineanalysis')
       fp.pap.opts{find(strcmp(fp.pap.opts,'Position'))+1} = fp.scr.pos(mod(f-1,nplotpos)+1,:);
   end
   set(fh.f(f),fp.pap.opts{:});
   
   % If figure has simple subplots, maximize them
   if     ~isfield(fh.a(f),'resize');   subplotResize(fh.a(f).h);
   elseif ~isempty(fh.a(f).resize  );   subplotResize(fh.a(f).h(fh.a(f).resize)); end
   
	% Set axes properties
   if f <= numel(fh.a)
      for a = 1:numel(fh.a(f).h)
         % if axis handle is in exclude, skip it
         if isfield(fh.a(f),'ex') && ~isempty(fh.a(f).ex) && any(fh.a(f).ex==fh.a(f).h(a))
            continue
         end
         % skip empty axes
         if isempty(get(fh.a(f).h(a),'Children'));   continue;   end
         
         % Change font size
         set(fh.a(f).h(a),fp.txt{end-1:end});
         
         % If axis handle is an image, reset fh.a(f).h(a) to its parent
         % (because imagesc outputs handle to image, not axis)
         if strcmpi(get(fh.a(f).h(a),'Type'),'image')
            fh.a(f).h(a) = get(fh.a(f).h(a),'Parent');
         elseif ~isfield(fp,'KEEPTICKPOS') || ~fp.KEEPTICKPOS
            % Fix axis ticks: only nticks number, max & min always present
            fh.a(f).h(a) = tickscale(fp.nticks,fh.a(f).h(a),fp.RESCALE,fp.signif);
         end
         
         % Tony-style axes
         if fp.TONY
            ud = get(fh.a(f).h(a),'UserData');
            if isempty(ud) || ~isstruct(ud) || ~isfield(ud,'Type') || ~strcmpi(ud.Type,'polarplot')
               newex = axesSeparate(fp,fh.a(f).h(a));
               if isfield(fh.a(f),'ex')
                  fh.a(f).ex = [fh.a(f).ex(:)' newex(:)'];
               else
                  fh.a(f).ex = newex(:)';
               end
            end
            if isfield(fp,'phyaxwidth') && ~isempty(fp.phyaxwidth) && fp.phyaxwidth~=1
               kids = get(fh.a(f).h(a),'Children');
               % Make axes thinner (paper)
               set(kids(~cellfun(@isempty,strfind(get(kids,'Tag'),'axis')) & ...
                        ~cellfun(@isempty,strfind(get(kids,'Tag'),'phyplot_')) | ...
                        strcmpi(get(kids,'Tag'),'phyplot_xtick') | ...
                        strcmpi(get(kids,'Tag'),'phyplot_ytick')         ),...
                        'LineWidth',fp.phyaxwidth)
            end
         else
            % Get rid of box
            try             set(fh.a(f).h(a),'Box','off');
            catch         % warning('this axis is boxless');
            end
            % Make tick direction out
            try             set(fh.a(f).h(a),'TickDir','out');
            catch         % warning('no ticks here');
            end
         
         end
      end

      % discard exceptions
      if isfield(fh.a(f),'ex')
         alltext = setdiff(alltext,fh.a(f).ex);
         allline = setdiff(allline,fh.a(f).ex);
%          allclin = setdiff(allclin,fh.a(f).ex); % AZ 2012-01-19 added
      end
            
      % Hide untouched axes
      tmp = get(fh.a(f).h,'Children');
      if ~isempty(tmp) && iscell(tmp)
         set(fh.a(f).h( cellfun(@isempty,tmp)),'Visible','off')
      elseif ~isempty(tmp) && ismatrix(tmp)
         set(fh.a(f).h(arrayfun(@isempty,tmp)),'Visible','off')
      end
   end
end


%% Exclude 3D figures from TONY-fication of dots
for f = includeFigs
   if f <= numel(fh.a)
      for a = 1:numel(fh.a(f).h)
         if any(get(fh.a(f).h(a),'ZLim') ~= [-1 1])
            if fp.TONY || ~fp.KEEPDOTSIZE
               warning('Tony-style dots disabled for 3D plot. Set fp.TONY=0 & fp.KEEPDOTSIZE=1');
            end
            % Exclude all dots in 3D plot
            fh.a(f).ex = [fh.a(f).ex(:)' findall(fh.a(f).h(a),'Marker','o')' ...
                                         findall(fh.a(f).h(a),'Marker','.')'    ];
         end
      end
   end
end

%% Set stuff
if ~isempty(allline)
   alldots = allline(~strcmpi(get(allline,'Marker'),'none'));
else
   alldots = [];
end

% If used plot(x,y,'k.'), transfer dot (line) color to marker face
if fp.TONY
   % Clothe naked dots with an edge color specified in fp.dots{4}
   nakeddots = alldots(strcmpi(get(alldots,'Marker'),'.'));
   
   if ~isempty(nakeddots)
      if numel(nakeddots)==1
         set(nakeddots,'MarkerFaceColor',get(nakeddots,'Color'))
      else
         cellfun(@(x,y) set(x,'MarkerFaceColor',y),num2cell(nakeddots),get(nakeddots,'Color'))
      end
      set(nakeddots,'Marker','o',fp.dots{3:4})
   end
   fp.dots = fp.dots([1:2 5:6]);
   
   % Make transparent 'o' dots have white facecolor
   emptydots = alldots(strcmpi(get(alldots,'MarkerFaceColor'),'none'));
   if ~isempty(emptydots)
      if numel(emptydots)==1
         set(emptydots,'MarkerFaceColor','w')
      else
         cellfun(@(x) set(x,'MarkerFaceColor','w'),num2cell(emptydots))
      end
   end
   
   % make dots pile on top of each other (need to set renderer to openGL)
   start = 1;
   for i = numel(alldots):-1:1
      finish = start + numel(get(alldots(i),'XData')) - 1;
      set(alldots(i),'ZData', start:finish);
      if i > 1 && get(alldots(i),'Parent') == get(alldots(i-1),'Parent')
         start = finish + 1;
      else
         start = 1;
      end
   end
end

littledots = alldots(strcmpi(get(alldots,'Marker'),'.'));
if fp.TONY && ~isempty(littledots)
   cellfun(@(x,y) set(x,'MarkerFaceColor',y),num2cell(littledots),get(littledots,'Color'))
   set(littledots,'Marker','o')
   % make dots pile on top of each other (need to set renderer to openGL)
   start = 1;
   for i = 1:numel(littledots)
      finish = start + numel(get(littledots(i),'XData')) - 1;
      set(littledots(i),'ZData', start:finish);
      if i < numel(littledots) && get(littledots(i),'Parent') == get(littledots(i+1),'Parent')
         start = finish + 1;
      else
         start = 1;
      end
   end
end

allline = setdiff(allline,alldots);
allline = setdiff(allline,littledots);
allline = setdiff(allline,findobj('-regexp','Tag','phyplot_.*'));
            
if ~isfield(fp,'KEEPLINEWIDTH') || ~fp.KEEPLINEWIDTH
   set(allline   ,fp.line{:});
   set(allpaxw   ,'LineWidth',fp.phyaxwidth);
   if ~isfield(get(allclin),'BarLayout') % don't change contours in bar graph
      try
         set(allclin,fp.line{:});
      catch
   %      warning('no contours here')
      end
   end
end
if ~isfield(fp,'KEEPDOTSIZE') || ~fp.KEEPDOTSIZE
   set(alldots   ,fp.dots{:});
   set(littledots,fp.dots{:});
end

set(alltext,fp.txt {:});

switch nargout
   case 1
      varargout = {fh};
   case 2
      varargout = {fh;fp};
end

end
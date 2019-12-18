function printFig(fighandle,outputDir,fn_ext,resolution,renderer,EXPORTFIG)
% printFig: Print a specific figure to a pdf
% figname comes from name property of figure
% 
% See also FIGSETUP.M BATCHPLOTREFINE.M TICKSCALE.M
%
% 2011-10-06 AZ Created
% 2011-10-11 AZ Generalized to different file types
% 2012-02-16 AZ Added renderer input option
% 2013-05-15 AZ Output filename extension more explicit
% 2013-05-24 AZ Text output more compact

% Initialization
%if ~exist('figname','var') || isempty(figname)
%	figname = get(fighandle(f),'Name');
%end
if ~exist('fighandle','var') || isempty(fighandle)
	fighandle = gcf;
elseif isstruct(fighandle) && isfield(fighandle,'f')
   fighandle = fighandle.f;
end

if ~exist('outputDir','var') || isempty(outputDir)
   outputDir = [pwd filesep];
end

if ~exist('resolution','var') || isempty(resolution)
	resolution = 300;
end

if ~exist('renderer','var') || isempty(renderer)
	renderer   = get(fighandle(1),'Renderer');
end
renderer   = lower(renderer);

resolution = sprintf('-r%u',round(abs(resolution)));
renderer   = sprintf('-%s',renderer);

special = [];
if exist('fn_ext','var') && ~isempty(fn_ext)
	switch lower(fn_ext)
		case 'pdf'
			outputType = '-dpdf';
		case 'svg'
			outputType = '-dsvg';
		case 'eps'
			outputType = '-depsc2';
         special    = '-loose';
		case 'ps'
			outputType = '-dpsc2';
		case 'png'
			outputType = '-dpng';
		case 'tiff'
			outputType = '-dtiff';
		otherwise
			outputType = '-dpdf';
			fn_ext     = 'pdf';
	end
else
	      outputType = '-dpdf';
	      fn_ext     = 'pdf';
end

if ~exist('EXPORTFIG','var') || isempty(EXPORTFIG)
   EXPORTFIG = false(numel(fighandle),1);
end


% Print!
fprintf('Printing all figures to: \n%s\n',outputDir);
for f = 1:numel(fighandle)
	figname = get(fighandle(f),'Name');
   fignum  = double(fighandle(f));
	fprintf('Printing figure %u with %s:\n%s ...',fignum,upper(outputType),figname);
	filename = sprintf('%s%s.%s',outputDir,figname,fn_ext);
%    DEBUG
%    get(gcf)
%    keyboard;
   if (numel(EXPORTFIG)>1 && ~EXPORTFIG(f)) || ~EXPORTFIG
      print(fighandle(f),outputType,filename,resolution,renderer,special);
   else
      % force load library for export_fig
      DIRS;
      if strcmpi(DIRS.host,'ZahootsBook')
         setenv('DYLD_INSERT_LIBRARIES','/opt/local/lib/libtiff.5.dylib:/opt/local/lib/libcurl.4.dylib');
      end
      export_fig(filename,['-' outputType(3:5)],renderer,fighandle(f));
   end
   % delete extraneous bounding boxes
   if strcmpi(fn_ext,'eps')
      system(sprintf('perl -0777 -i -pe ''s/[0-9\\.]{1,8}\\s[0-9\\.]{1,8}\\s[0-9\\.]{1,8}\\s[0-9\\.]{1,8}\\sre\\n[Wf]\\n//gd'' %s',filename));
   end

	fprintf('done.\n');
end

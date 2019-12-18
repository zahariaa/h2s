% figsetup.m
% Sets up MATLAB for nice figures, dumps params into struct fp
% 
% See also PRINTFIG.M BATCHPLOTREFINE.M TICKSCALE.M
% 
% 2011-10-06 AZ Created
% 2012-01-20 AZ Added renderer field
% 2012-01-30 AZ Added quarters fields
% 2012-02-15 AZ Default TickDir Out, added number of ticks (nticks) field
% 2012-03-14 AZ Added PaperPositionMode=manual
% 2012-06-15 AZ Initialize figname, fh
% 2012-07-18 AZ added sixths.pos
% 2013-03-15 AZ added onlineanalysis condition

fh          = struct('a',[],'f',[]); % get this out of the way

v = version;
fp.HG2 = str2double(v(1:3))>=8.4;
clear v

% flags
if ~exist('fp','var') || ~isfield(fp,'TONY')
   fp.TONY         = true;
end

ASPECTRATIO = [];
if ~isfield(fp,'DESTINATION')   || isempty(fp.DESTINATION)
   fp.DESTINATION = {'presentation'};
elseif iscell(fp.DESTINATION) && numel(fp.DESTINATION)>1
   ASPECTRATIO    = fp.DESTINATION{2};
elseif ischar(fp.DESTINATION)
   fp.DESTINATION = {fp.DESTINATION};
end
if ASPECTRATIO;   ASPECTRATIO = flip(ASPECTRATIO)/ASPECTRATIO(2);   end

if isfield(fp,'pap') && isfield(fp.pap,'sz') && ~isempty(fp.pap.sz)
   fp.FORCEPAPSIZE = true;
else
   fp.FORCEPAPSIZE = false;
end

if ~isfield(fp,'KEEPLINEWIDTH') || isempty(fp.KEEPLINEWIDTH)
   fp.KEEPLINEWIDTH = false;
end
if ~isfield(fp,'KEEPTICKPOS')   || isempty(fp.KEEPTICKPOS)
   fp.KEEPTICKPOS   = false;
end
if ~isfield(fp,'KEEPDOTSIZE')   || isempty(fp.KEEPDOTSIZE)
   fp.KEEPDOTSIZE = false;
end

fp.MINORTICKS = false;
fp.signif = 0;

% setup
set(0,'DefaultFigureColor','w');
set(0,'DefaultAxesTickDir','out');
% next line doesn't really work, batchPlotRefine explicitly enforces no boxes
set(0,'DefaultAxesBox','off');
fp.scr.sz = get(0,'ScreenSize');
fp.scr.sz = fp.scr.sz(3:4);
if exist('DIRS','var') && strcmpi(DIRS.host,'gunther')
    fp.scr.sz = [1440 900];
end

% % OpenGL needed for lighting/shading, painters higher res, less buggy (?)
% fp.renderer = 'painters'; %painters, zbuffer, OpenGL

% quarters.pos: TL; TR; BL; BR
fp.quarters.pos = [0   1/2 1/2 1/2
                   1/2 1/2 1/2 1/2
                   0     0 1/2 1/2
                   1/2   0 1/2 1/2];
% sixths.pos: TL; TC; TR; BL; BC; BR
fp.sixths.pos   = [0   1/2 1/3 1/2
                   1/3 1/2 1/3 1/2
                   2/3 1/2 1/3 1/2
                   0     0 1/3 1/2
                   1/3   0 1/3 1/2
                   2/3   0 1/3 1/2];

fp.quarters.campos = [ 8  0  0;
                       0  0  8;
                       1 -1  8;
                       0 -8  0];

fp.pap.opts = {'MenuBar','none','PaperUnits','inches','PaperSize',[],'PaperPosition',[],...
               'PaperPositionMode', 'manual','PaperOrientation','portrait','Renderer','painters'};
% PARAMETERS
fp.line = {'LineWidth',5};
fp.dots = {'MarkerSize',14.4,'MarkerEdgeColor',[1 1 1],'LineWidth',1};

switch lower(fp.DESTINATION{1})
case {'paper';'portrait';'landscape';'1col';'1.5col';'2col';'elife';'elifefull';'elifehalf'}
   fp.pap.u   = 'inches';
   if ~fp.FORCEPAPSIZE
      switch lower(fp.DESTINATION{1})                         %  (inches)  [width x height]
         case 'landscape'; fp.pap.opts{12}=fp.DESTINATION{1}; fp.pap.sz  = [11        8.5 ];
         case 'portrait';  fp.pap.opts{12}=fp.DESTINATION{1}; fp.pap.sz  = [ 8.5     11   ];
         case '1col';                                         fp.pap.sz  = [ 3.35     3.35];
         case '1.5col';                                       fp.pap.sz  = [ 4.57     4.57];
         case '2col';                                         fp.pap.sz  = [ 6.93     6.93];
         case 'elife';                                        fp.pap.sz  = [ 5.5      5.5 ];
         case 'elifehalf';                                    fp.pap.sz  = [ 5.5      5.5 ]*0.46;
         case 'elifefull';                                    fp.pap.sz  = [ 7.28     7.28];
         otherwise;                                           fp.pap.sz  = [ 1.5      1.5 ];
      end
      if  ASPECTRATIO;    fp.pap.sz  = fp.pap.sz.*ASPECTRATIO;   end
   end
	fp.scr.pos = round([fp.scr.sz*.1 (fp.pap.sz(1)*fp.scr.sz(2)*.8/fp.pap.sz(2)) fp.scr.sz(2)*.8]);
	fp.txtsz   = 7.5;
   fp.dots{2} = 6;
   fp.dots{6} = 0.5;
   fp.line{2} = 1;
case 'poster'
   fp.pap.u   = 'inches';
   if ~fp.FORCEPAPSIZE;   fp.pap.sz  = [4 4];                    end
   if     ASPECTRATIO;    fp.pap.sz  = fp.pap.sz.*ASPECTRATIO;   end
	fp.scr.pos = round([fp.scr.sz*.1 (fp.pap.sz(1)*fp.scr.sz(2)*.8/fp.pap.sz(2)) fp.scr.sz(2)*.8]);
	fp.txtsz   = 24;
   fp.dots{6} = 1;
case 'presentation'
	fp.pap.u   = 'points';
   if ~fp.FORCEPAPSIZE;   fp.pap.sz  = [1024 1024];              end
   if     ASPECTRATIO;    fp.pap.sz  = fp.pap.sz.*ASPECTRATIO;   end
	fp.scr.pos = round([fp.scr.sz*.1 (fp.pap.sz(1)*fp.scr.sz(2)*.8/fp.pap.sz(2)) fp.scr.sz(2)*.8]);
	fp.txtsz   = 24;
   
case 'onlineanalysis'
    fp.pap.u   = 'points';
    if ~fp.FORCEPAPSIZE
	 fp.pap.sz  = fp.scr.sz/6;
%  	fp.scr.pos = round([fp.scr.sz*.1 (fp.pap.sz(1)*fp.scr.sz(2)*.8/fp.pap.sz(2)) fp.scr.sz(2)*.8]);
%     fp.scr.pos = repmat(fp.scr.sz,[4 2]).*fp.quarters.pos;
    end
    if     ASPECTRATIO;   fp.pap.sz  = fp.pap.sz.*ASPECTRATIO;   end
    fp.scr.pos = repmat(fp.scr.sz,[6 2]).*fp.sixths.pos;
    if strcmpi(DIRS.host,'gunther')
        fp.scr.pos = fp.scr.pos + repmat([1600+40 1200 0 -30],[6 1]);
    end
   fp.line    = {'LineWidth',3};
   fp.dots{2} = 9;
	fp.txtsz   = 16;
otherwise
   if ~isfield(fp.pap,'u')
      fp.pap.u   = 'points';
   end
   if ~isfield(fp.pap,'sz') || ~fp.FORCEPAPSIZE
	fp.pap.sz  = [1024 1024];
   end
   if     ASPECTRATIO;    fp.pap.sz  = fp.pap.sz.*ASPECTRATIO;   end
	fp.scr.pos = round([fp.scr.sz*.1 (fp.pap.sz(1)*fp.scr.sz(2)*.8/fp.pap.sz(2)) fp.scr.sz(2)*.8]);
   fp.line    = {'LineWidth',3};
   fp.dots{2} = 9;
   if ~isfield(fp,'txtsz')
      fp.txtsz   = 20;
   end
end
fp.phyaxwidth = fp.dots{6};

if ~isfield(fp.pap,'u') || isempty(fp.pap.u)
   fp.pap.u   = fp.pap.opts{4};
end

fp.pap.pos = [fp.pap.sz(1:2)*0.01                fp.pap.sz];
fp.nticks = [4 4];%[5 5 5]
fp.RESCALE = false;

fp.pap.opts{4} = fp.pap.u;
fp.pap.opts{6} = fp.pap.sz;
fp.pap.opts{8} = fp.pap.pos;

% EXECUTE
if fp.TONY
	set(0,'DefaultAxesFontName' , 'Helvetica');
   if strcmpi(fp.DESTINATION{1},'paper') || ~isempty(strfind(fp.DESTINATION{1},'col'  )) ...
                                         || ~isempty(strfind(fp.DESTINATION{1},'elife'))
   set(0,'DefaultAxesFontAngle', 'Normal'   );
   else
   set(0,'DefaultAxesFontAngle', 'Oblique'  );
   end
%	set(0,'DefaultFontWeight'   , 'Normal'   );
	% add fp.ax{:} to axis properties, and xlabel/ylabel
	fp.ax  = {'FontName'     , get(0,'DefaultAxesFontName') ,...
             'AxesFontAngle', get(0,'DefaultAxesFontAngle'),'FontWeight','Normal'};
	fp.txt = {'FontName'     , get(0,'DefaultAxesFontName') ,...
                 'FontAngle', get(0,'DefaultAxesFontAngle'),'FontWeight','Normal',...
             'FontSize',fp.txtsz};
end

fp.pap.opts = [fp.pap.opts,{'Position',fp.scr.pos(1,:)}];

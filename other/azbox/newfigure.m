function [fh,f] = newfigure(NSUBPLOTSalso,fh,fp,figname,NSUBPLOTS,holdFlag)
% newfigure.m: Set up plots, AZ-style
% 
% See also FIGSETUP, PRINTFIG, BATCHPLOTREFINE, AXTIVATE
% 
% 2012-11-14 AZ Created
% 2013-02-15 AZ Added DeleteFcn

if ~exist('fh'      ,'var') || isempty(fh);        figsetup;
elseif ischar(fh) && ( ~exist('figname' ,'var') || isempty(figname ) )
   figname = fh;   figsetup;
elseif isnumeric(fh) && ( ~exist('figname' ,'var') || isempty(figname ) )
   NSUBPLOTS = fh; figsetup;
end
ftmp = figure;
if verLessThan('matlab','8.4'), fh.f=[fh.f ftmp];
else                            fh.f=[fh.f ftmp.Number];
end
f = numel(fh.f);

% NSUBPLOTSalso is now a dummy variable (for backwards compatibility).
% NEEDS TESTING
if exist('NSUBPLOTSalso','var') && ~isempty(NSUBPLOTSalso)
   if isnumeric(NSUBPLOTSalso) && ...
      ( ~exist('NSUBPLOTS','var') || isempty(NSUBPLOTS) )
      NSUBPLOTS = NSUBPLOTSalso;
   elseif ischar(NSUBPLOTSalso) && ( ~exist('figname' ,'var') || isempty(figname ) )
      figname = NSUBPLOTSalso;
   end
end

if ~exist('figname' ,'var') || isempty(figname );   figname  = char(32);  end
if ~exist('holdFlag','var') || isempty(holdFlag)
                 holdFlag = 'on';
elseif islogical(holdFlag)
   if holdFlag;  holdFlag = 'on';
   else          holdFlag = 'off';
   end
end

fh.n{f} = figname;
set(fh.f(f),'Name',fh.n{f},'DeleteFcn',@closeFH);
if exist('fp','var') && ~isempty(fp)
   if (~exist('NSUBPLOTS','var') || isempty(NSUBPLOTS)) && ...
         iscell(fp.DESTINATION) && numel(fp.DESTINATION)>1
      NSUBPLOTS = fp.DESTINATION{2};
   end
   set(fh.f(f),fp.pap.opts{:},'Position',fp.scr.pos);
end

% Make axes if desired
if ~exist('NSUBPLOTS','var')
   fh.a(f).h  = gca; hold(holdFlag);
   fh.a(f).ex = [];
elseif isnumeric(NSUBPLOTS) && numel(NSUBPLOTS)==2
   a = 0;
   for i = 1:NSUBPLOTS(1)
      for j = 1:NSUBPLOTS(2)
         a = a + 1;
         fh.a(f).h(a) = subplot(NSUBPLOTS(1),NSUBPLOTS(2),a); hold(holdFlag);
         fh.a(f).ex = [];
      end
   end
else
   fh.a(f).h  = [];
   fh.a(f).ex = [];
end

% Get calling workspace
[wstack,i] = dbstack;
if size(wstack) == 1
   wstack = 'base';
else
   wstack = wstack(i+1).name;
end
% Embed calling workspace in figure
set(fh.f(f),'UserData',wstack);

end



function closeFH(src,~)

varsInBase = evalin('base',...
   'exist(''fh'',''var'') && ~exist(''fTMP'',''var'') && isfield(fh,''f'') && ~isempty(fh.f)');

if varsInBase
   assignin('base','fTMP',src);
   try
      evalin('base','fh.f(fTMP) = [];fh.a(fTMP) = [];');
   catch
   end
   evalin('base','clear fTMP');
end

end

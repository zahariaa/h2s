function [fh,f] = newfigure(varargin)
% newfigure.m: Set up plots, AZ-style
% Inputs:
%    figname (DEFAULT = ' ')
%    NSUBPLOTS (DEFAULT = [1 1])
%    holdFlag (DEFAULT = true)
%    fh, fp (DEFAULT = figsetup)
% 
% See also FIGSETUP, PRINTFIG, BATCHPLOTREFINE, AXTIVATE
% 
% 2012-11-14 AZ Created
% 2013-02-15 AZ Added DeleteFcn
% 2020-06-09 AZ Simplified and made more flexible with varargin

for v = 1:nargin
   if     ischar(varargin{v}),        figname = varargin{v};
   elseif isstruct(varargin{v})
      if     isfield(varargin{v}, 'a' ), fhin = varargin{v};
      elseif isfield(varargin{v},'pap'),   fp = varargin{v};
      end
   elseif isnumeric(varargin{v}),   NSUBPLOTS = varargin{v};
   elseif islogical(varargin{v}),    holdFlag = varargin{v};
   end
end

if ~exist('fp','var') || isempty(fp);      figsetup; end
if exist('fhin','var') && ~isempty(fhin); fh = fhin; end

ftmp = figure;
if verLessThan('matlab','8.4'), fh.f=[fh.f ftmp];
else                            fh.f=[fh.f ftmp.Number];
end
f = numel(fh.f);

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

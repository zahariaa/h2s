function ax = tickscale(nticks,ax,RESCALE,signif)
% TICKSCALE.M
% Forces axis ax to have nticks number of ticks, spanning max & min of axis
% ax     is the axis handle
% nticks is a 2- or 3-vector defining how many ticks should be in the [x y]
%        or [x y t] axes
%
% 2012-02-15 AZ Created
% 2013-05-15 AZ Disabled warning

if ~exist('nticks' ,'var') || isempty(nticks);   nticks  = [5 5];  end
if ~exist('ax'     ,'var') || isempty(ax);       ax      = gca;    end
if ~exist('RESCALE','var') || isempty(RESCALE);  RESCALE = false;  end
if ~exist('signif' ,'var') || isempty(signif);   signif  = 0;      end

xyz = 'xyz';

for i = 1:3
   if numel(nticks)>(i-1) && nticks(min(i,2)) >= 0 && ...
      strcmpi(get(ax,[xyz(i) 'Scale']),'linear') && ~isempty(get(ax,[xyz(i) 'Tick']))
      axlim = get(ax,[xyz(i) 'Lim']);

      roundoff = 10^(signif-floor(log10(max(abs(axlim)))));
      
      % rescale axes to have nice bounds
      if RESCALE
         axlim = round( roundoff*axlim ) / roundoff;
         set(ax,[xyz(i) 'Lim'],axlim);
      end
      roundoff = roundoff * 10;
      
      % if 0 is stictly inside axis limits include, narrow options to ticks with 0
      if sum(sign(axlim)) == 0
         ix = arrayfun(@(NTICKS) ...
               sum(round( roundoff*linspace(axlim(1),axlim(2),NTICKS) )/roundoff==0),...
               3:nticks(i));
         if RESCALE && sum(ix)>0
            newticks = round( roundoff*linspace(axlim(1),axlim(2),nticks(min(i,2))) )/roundoff;
            set(ax,[xyz(i) 'Tick'],newticks);
            newticks = [];
            continue
         elseif RESCALE || sum(ix) == 0
            % if no nice ones found, force axlim to be {min,0,max}
            newticks = round( roundoff*[axlim(1) 0 axlim(2)])/roundoff;
            % if not symmetric, force symmetry
            switch sign(sum(newticks))
               case  1; newticks(3) = -newticks(1);
               case -1; newticks(1) = -newticks(3);
               case  0; % do nothing, already symmetric
            end
            set(ax,[xyz(i) 'Tick'],newticks);
            newticks = [];
            continue
         end
      end
      
      % override given nticks if integer valued ticks are possible
      if RESCALE
            newticks = round( roundoff*linspace(axlim(1),axlim(2),nticks(min(i,2))) )/roundoff;
      else
         if ~exist('newticks','var') || isempty(newticks)
            % Try all number of ticks from 2 to nticks(i)
            newticks = arrayfun( @(x)  linspace(axlim(1),axlim(2),x),...
                                       2:nticks(min(i,2)),'UniformOutput',false);
            % Delete ticks that would have been rounded
            newticks = cellfun(@(x) x(mod(roundoff*x,1)==0),newticks,'UniformOutput',false);
            % Choose sequence with maximum number of rounded ticks surviving
            newticks = newticks{maxix(cellfun(@numel,newticks))};
         end
      end
      
      % Getting desperate: sample Matlab's chosen ticks
      if isempty(newticks) || numel(newticks)==1
         currticks = get(ax,[xyz(i) 'Tick']);
         ix = arrayfun( @(x)  1:x:numel(currticks),2:numel(currticks),'UniformOutput',false);
         % choose 
         while numel(newticks) > nticks(min(i,2)) || numel(newticks)<2 && ~isempty([ix{:}])
            found     = maxix(cellfun(@numel,ix));
            newticks  = currticks(ix{found});
            ix{found} = [];
         end
      end
      
      if sum(sign(axlim)) == 0 && (~any(newticks==0) || numel(newticks) < nticks(i))
         % if no nice ones found, force axlim to be {min,0,max}
         newticks = round( roundoff*[axlim(1) 0 axlim(2)])/roundoff;
         set(ax,[xyz(i) 'Tick'],unique(newticks));
      end
      
      if numel(unique(newticks)) == numel(newticks) && ~isempty(newticks)
         set(ax,[xyz(i) 'Tick'],newticks);
      end
      newticks = [];
   end
end

% if axis is square, force ticks to be the same. take lower number of ticks
d = daspect; pb = pbaspect;
if all(abs(get(ax,'XLim')-get(ax,'YLim'))<eps) && d(1)==d(2) && pb(1)==pb(2)
   xticks = get(ax,'XTick');
   yticks = get(ax,'YTick');
   if     numel(xticks) < numel(yticks); set(ax,'YTick',xticks);
   elseif numel(xticks) > numel(yticks); set(ax,'XTick',yticks);
   elseif ~all(xticks==yticks)           % do nothing?
   end
end
% keyboard

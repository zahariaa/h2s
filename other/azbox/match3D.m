function match3D(sourceAx,targetAx,includeOP,felds)
% match3D(sourceAx,targetAx)
% match views of 3D plots
% 
% 2017-12-20 AZ Created

if ~exist('felds','var') || isempty(felds)
   felds = {'XLim','YLim','ZLim','DataAspectRatio'};
end
if exist('includeOP','var') && ~isempty(includeOP) && includeOP
   felds = [felds,{'OuterPosition'}];
else includeOP = false;
end

% Recurse
if numel(targetAx)>1
   if numel(sourceAx)>1;  arrayfun(@(s,t) match3D(s,t,includeOP,felds),sourceAx,targetAx)
   else                   arrayfun(@(  t) match3D(sourceAx,t,includeOP,felds)  ,targetAx)
   end
   return
end

for fn=felds
   set(targetAx,fn,get(sourceAx,fn))
end


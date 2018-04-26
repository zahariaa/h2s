function s = h2s_center(n,d,type,PLOT)

%% Preliminaries
if ~exist('PLOT','var') || isempty(PLOT),   PLOT = false;   end
colors = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1],[0 1 1]};
ne = 4;
N = 100; % bootstraps

%% Recurse if more than one dimension or number of samples input argument given
nd = numel(d);
nn = numel(n);
if nd > 1
   for ix = 1:nd
      s{ix} = h2s_center(n,d(ix),type,PLOT);
      if PLOT && ix>1,   plotCenterEsts(d(1:ix),s,colors,type);   end
      xlabel('dimensions')
   end
   return
end
if nn > 1
   for ix = 1:nn
      s{ix} = h2s_center(n(ix),d,type,PLOT);
      if PLOT && ix>1,   plotCenterEsts(n(1:ix),s,colors,type);   end
      xlabel('samples')
   end
   return
end

s = NaN(d,N,ne);

%% Sample & measure
X = sampleSpheres(n,d,N,type);
%% Estimators
s(:,:,1) = mean(X);
s(:,:,2) = median(X);
s(:,:,3) = mean(s(:,:,1:2),3);
cv = cvindex(n,10);
for i = 1:N
   s(:,i,4) = mean(cell2mat_concat(cv.crossvalidate(@mean,X(:,:,i))));
end
% try
% for i = 1:N
%    K = convhulln(X(:,:,i));
% %   [F,ix] = getFurthestPoints(X(:,:,i),d+1);
%    s(:,i,4) = mean(X(K(:,1),:,i));
% %   s(:,i,5) = mean(X(ix    ,:,i));
% end
% catch
% end

if PLOT,   figure(98);clf;plotEstimators(X,zeros(d,N),s,colors);   end

return

%% PLOT function
function plotCenterEsts(ds,s,colors,type)

figure(99);clf;hold on;
for e = 1:ne
   for ix = 1:numel(s)
      err(:,ix) = mean(squeeze(s{ix}(:,:,e)).^2);
   end
   plotErrorPatch(log2(ds),log10(err),colors{e});
end
ylabel('log_{10}(error^2)')
title([type ' distribution, center estimation'])
legend('','mean','','median',...
       '','mean of mean & median','','cv mean',...
       'Location','EastOutside')
set(legend,'Box','off')
return
end

%% DEMOS/DEBUG
h2s_center(200,log2space(1,12,12),'Uniform',true);
set(gca,'XLim',[1 12],'XTick',1:3:12,'XTickLabel',num2str(2.^(1:3:12)'))
set(99,'Name','CenterUniform','renderer','painters');
delete(findobj(get(legend,'Children'),'type','patch'))
axesSeparate;printFig(99,'~/Desktop/','eps');

h2s_center(200,log2space(1,12,12),'Gaussian',true);
set(gca,'XLim',[1 12],'XTick',1:3:12,'XTickLabel',num2str(2.^(1:3:12)'))
set(99,'Name','CenterGaussian','renderer','painters');
delete(findobj(get(legend,'Children'),'type','patch'))
axesSeparate;printFig(99,'~/Desktop/','eps');



end

function plotEstimators(data,varToHist,estimates,colors,fig)
% plotEstimators(data,varToHist,estimates,colors,fig=gcf)
% plots histogram of estimated variable varToHist,
% and compares {estimates} by plotting them in {colors}
% 
% 2018-04-13 AZ Created

if ~exist('fig','var') || isempty(fig),   fig = gcf;   end
if exist('estimates','var') && ~isempty(estimates) && isnumeric(estimates)
   estimates = {estimates};
end
if exist('colors','var') && ~isempty(colors) && isnumeric(colors)
   colors = {colors};
end

% make CDF
varToHist = sort(varToHist,1,'ascend');
cdf = cumsum(~~varToHist,1)/size(data,1);
% PCA, for visualization
P = simplepca(data(:,:,1));

%% PLOT
subplot(1,3,2:3);hold on;
[h,c]=hist(varToHist(:));bar(c,h);
%[h,c]=hist(varToHist(:),numel(varToHist)/100);barh(c,h);
plot(squeeze(varToHist),1.1*max(h)*squeeze(cdf),'k-')

for e = 1:numel(estimates)
   plotErrorPatch([estimates{e}(:) estimates{e}(:)] ,[0 1.1*max(h)],colors{e});
end
xlabel('L2 norm')
%title({type;
%       sprintf('%u points in %u dimensions',n,d)});

subplot(2,3,1);hold on;
plot(data(:,1,1),data(:,2,1),'ko');
for e = 1:numel(estimates)
   plot([0 mean(estimates{e}(:))]./sqrt(2),[0 mean(estimates{e}(:))]./sqrt(2),'Color',colors{e},'LineWidth',16);
end
axis equal square
drawCircle;
title({'r=expDistPerRad method,';
       'g=median*2^{1/d}, b=median' })

subplot(2,3,4);hold on;
plot(P(:,1),P(:,2),'ko');
for e = 1:numel(estimates)
   plot([0 mean(estimates{e}(:))]./sqrt(2),[0 mean(estimates{e}(:))]./sqrt(2),'Color',colors{e},'LineWidth',16);
end
axis equal square
drawCircle;
title('PCA')


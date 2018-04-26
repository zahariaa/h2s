function plotEstimators(data,varToHist,estimates,colors,type,fig)
% plotEstimators(data,varToHist,estimates,colors,fig=gcf)
% plots histogram of estimated variable varToHist,
% and compares {estimates} by plotting them in {colors}
% 
% 2018-04-13 AZ Created
if ndims(data)==3,   [n,d,N] = size(data);
else                 [n,N]   = size(data);
end

if ~exist('fig','var') || isempty(fig),   fig = gcf;   end
if exist('estimates','var') && ~isempty(estimates) && isnumeric(estimates)
   ne = last(size(estimates));
   if ndims(estimates)==3
      type = 'center'; % Determine estimate type by dimensionality of estimates
      estimates = squeeze(mat2cell(estimates,d,N,ones(1,ne)));
      % Estimates mean of center estimates
      meanCenters = cellfun(@(x) mean(x,2),estimates,'UniformOutput',false);
      % Quick and dirty radius estimate of center estimates, calculated from means above
      % NOTE: these are essentially "100% confidence intervals"
      varCenters  = cellfun(@(x,m) median(sqrt(2*mean((x-repmat(m,[1 N])).^2)))*2^(1/d),estimates,meanCenters);
   else
      type =  'radii'; % Determine estimate type by dimensionality of estimates
      estimates = mat2cell(estimates,N,ones(1,ne));
   end
end
if exist('colors','var') && ~isempty(colors) && isnumeric(colors)
   colors = {colors};
end

H = hyperdist2dist(data(:,:,1)); % for visualization

%% PLOT L2 norms
subplot(1,3,2:3);hold on;
if strcmpi(type,'radii')
   % make CDF
   varToHist = sort(varToHist,1,'ascend');
   cdf = cumsum(~~varToHist,1)/size(data,1);
   %[h,c]=hist(varToHist(:));bar(c,h);
   [h,c]=hist(varToHist(:),numel(varToHist)/100);bar(c,h);
   plot(squeeze(varToHist),1.1*max(h)*squeeze(cdf),'r-')
   for e = 1:numel(estimates)
      plotErrorPatch([estimates{e}(:) estimates{e}(:)] ,[0 1.1*max(h)],colors{e});
   end
else
   for e = 1:numel(estimates)
      plotErrorPatch(repmat(log10(sum((varToHist - estimates{e}).^2)),[2 1])',[0 1],colors{e});
   end
end
xlabel('L2 norm')
title({type; sprintf('%u points in %u dimensions',n,d)});

%% PLOT first 2 dimensions and estimates
subplot(2,3,1);hold on;
plot(data(:,1,1),data(:,2,1),'ko');
for e = 1:numel(estimates)
   if strcmpi(type,'radii')
      plot([0 mean(estimates{e}(:))]./sqrt(2),...
           [0 mean(estimates{e}(:))]./sqrt(2),'Color',colors{e},'LineWidth',16);
   elseif strcmpi(type,'center')
      plot(0,0,'kx')
      plot(meanCenters{e}(1),meanCenters{e}(2),'x','Color',colors{e});
      drawCircle(meanCenters{e}(1:2),varCenters(e),colors{e});
   end
end
axis equal square off
if strcmpi(type,'radii');   drawCircle;   end
title({'r=expDistPerRad method,';
       'g=median*2^{1/d}, b=median' })

%% PLOT hyperdist2dist and estimates
subplot(2,3,4);hold on;
plot(H(:,1),H(:,2),'ko');
for e = 1:numel(estimates)
   plot([0 mean(estimates{e}(:))]./sqrt(2),...
        [0 mean(estimates{e}(:))]./sqrt(2),'Color',colors{e},'LineWidth',16);
end
axis equal square
drawCircle;
title('hyperdist2dist')


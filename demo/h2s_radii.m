function [s,target] = h2s_radii(n,d,type,PLOT)

%% Preliminaries
if ~exist('PLOT','var') || isempty(PLOT),   PLOT = false;   end
colors = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1],[0 1 1]};
ne = numel(colors);
N = 100; % bootstraps
s = NaN(N,ne); % container for sampled estimators

%% Recurse if more than one dimension or number of samples input argument given
nd = numel(d);
nn = numel(n);
if nd > 1
   s = repmat(s,[1 1 nd]);
   target = NaN(1,nd);
   for ix = 1:nd
      [s(:,:,ix),target(ix)] = h2s_radii(n,d(ix),type,PLOT);
      if PLOT,   plotRadiusEsts(d,s,target,colors,type);   end
      xlabel('dimensions')
   end
   return
end
if nn > 1
   s = repmat(s,[1 1 nn]);
   target = NaN(1,nn);
   for ix = 1:nn
      [s(:,:,ix),target(ix)] = h2s_radii(n(ix),d,type,PLOT);
      if PLOT,   plotRadiusEsts(n,s,target,colors,type);   end
      xlabel('samples')
   end
   return
end

%% Sample & measure radii
X  = sampleSpheres(n,d,N,type);
for i = 1:N;   meanDists(i) = mean(pdist(X(:,:,i),'Euclidean'));   end

% Center X
Xc = X - repmat(mean(X,1),[n 1 1]);
radii = sqrt(sum(Xc.^2,2));
maxradii = max(radii,[],1);
% Target to measure estimates against
if     strcmpi(type,'gaussian'),   target = median(radii(:));
elseif strcmpi(type,'uniform' ),   target = 1;
end

%% Estimators
s(:,1) = meanDists/expDistPerRad(d);
s(:,3) = median(radii,1);
s(:,2) = s(:,3).*2^(1/d);   % median*2^(1/d)
s(:,4) = maxradii + maxradii*(n^-d);                   % MVUE for uniform distribution
v = std(reshape(Xc,[n*d N])); % squeeze(mean(std(X,[],2))); %v(i) = max(diag(cov(X(:,:,i))));
s(:,5) = v*sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2)); % MVUE for gaussian distribution
[~,~,tmp] = arrayfun(@(i) estimateHypersphere(X(:,:,i)),1:N,'UniformOutput',false);
s(:,6) = cell2mat(tmp);


if PLOT,   figure(98);clf;plotEstimators(Xc,radii/target,s/target,colors);   end

return


%% PLOT function
function plotRadiusEsts(ds,s,target,colors,type)

figure(99);clf;hold on;
for e = 1:ne
   err = (repmat(target,[N 1])-squeeze(s(:,e,:)))./repmat(target,[N 1]);
   plotErrorPatch(log2(ds),log10(err.^2),colors{e});
end
ylabel('log_{10}(error^2)')
title([type ' distribution, radius estimation'])
legend('','expDistPerRad','','median2^{1/d}',...
       '','median','','MVUE-Unif','','MVUE-Gauss',...
       'Location','EastOutside')
set(legend,'Box','off')
return
end


%% expDistPerRad function
function r = expDistPerRad(d)
% by interpolation among estimates precomputed in QT_expectedDistBetweenTwoPointsWithinHypersphere
ds= [1 2 3 5 8 15 30 60 100 1000];
expectedDists = [0.6582    0.9073    1.0261    1.1490    1.2343    1.3181    1.3651    1.3877    1.3993  1.4127];
if d<1000
    r = interp1(ds,expectedDists,d);
else
    r = sqrt(2);
end
return
end


%% DEMOS/DEBUG
h2s_radii(200,log2space(1,12,12),'Uniform',true);
set(gca,'XLim',[1 12],'XTick',1:3:12,'XTickLabel',num2str(2.^(1:3:12)'))
set(99,'Name','RadiusUniform','renderer','painters');
delete(findobj(get(legend,'Children'),'type','patch'))
axesSeparate;printFig(99,'~/Desktop/','eps');

h2s_radii(200,log2space(1,12,12),'Gaussian',true);
set(gca,'XLim',[1 12],'XTick',1:3:12,'XTickLabel',num2str(2.^(1:3:12)'))
set(99,'Name','RadiusGaussian','renderer','painters');
delete(findobj(get(legend,'Children'),'type','patch'))
axesSeparate;printFig(99,'~/Desktop/','eps');

end



function [s,target,ev] = h2s_radii(n,d,type,PLOT,SAVE)

%% Preliminaries
if ~exist('PLOT','var') || isempty(PLOT),   PLOT = false;   end
if ~exist('SAVE','var') || isempty(SAVE),   SAVE = false;   end
colors = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1],[0 1 1]};
ne = numel(colors);
N = 100; % bootstraps

%% Recurse if more than one dimension or number of samples input argument given
nd = numel(d);
nn = numel(n);
if nd > 1
   exptType = 'dimensions';
   s = NaN(N,ne,nd); % container for sampled estimators
   target = NaN(1,nd);
   for ix = 1:nd
      [s(:,:,ix),target(ix),ev(:,ix)] = h2s_radii(n,d(ix),type,PLOT);
      if PLOT,   plotRadiusEsts(d,s,target,colors,type);   end
      if SAVE,   saveRadiusEsts(type,exptType);   end
      xlabel(exptType)
   end
   return
end
if nn > 1
   exptType = 'samples';
   s = NaN(N,ne,nn); % container for sampled estimators
   target = NaN(1,nn);
   for ix = 1:nn
      [s(:,:,ix),target(ix)] = h2s_radii(n(ix),d,type,PLOT);
      if PLOT,   plotRadiusEsts(n,s,target,colors,type);   end
      if SAVE,   saveRadiusEsts(type,exptType);   end
      xlabel(exptType)
   end
   return
end
s = NaN(N,ne); % container for sampled estimators

%% Sample & measure radii
X  = sampleSpheres(n,d,N,type);
%X = X + 1*randn(size(X))/(sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2)));
for i = 1:N;   meanDists(i) = mean(pdist(X(:,:,i),'Euclidean'));   end

% Center X
Xc = X - repmat(mean(X,1),[n 1 1]);
radii = sqrt(sum(Xc.^2,2));
maxradii = max(radii,[],1);
% Target to measure estimates against
if     strcmpi(type,'gaussian'),   target = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));%median(radii(:));
elseif strcmpi(type,'uniform' ),   target = 1;
else                               target = hypercubeRadius(d);
end

%% Estimators
s(:,1) = meanDists/expDistPerRad(d);
s(:,3) = median(radii,1);
s(:,2) = s(:,3).*2^(1/d);   % median*2^(1/d)
% MVUEs
if     strcmpi(type,'gaussian'), s(:,4) = std(reshape(Xc,[n*d N]))*target;
elseif strcmpi(type,'uniform' ), s(:,4) = maxradii + maxradii*(n^-d);
else                             s(:,4) = s(:,3);
end

%% Testing ML joint center & radius estimator
tol = 10^(-7-log2(d));
opts = optimoptions('fminunc','TolX',tol,'TolFun',tol,'Algorithm','quasi-newton','Display','off','GradObj','on');%,'DerivativeCheck','on');
mest = arrayfun(@(i) fminunc(@(m) maxRadiusGivenCenter(m,Xc(:,:,i)),mean(Xc(:,:,i)),opts),1:N,'UniformOutput',false);
maxradii = arrayfun(@(i) maxRadiusGivenCenter(mest{i},Xc(:,:,i)),1:N)';
s(:,3) = maxradii;

% s(:,2) = arrayfun(@(i) mean(undeal(2,@() inferHyperspherePosterior(X(:,:,i),1000,false))),1:N);
% s(:,5) = arrayfun(@(i) mean(undeal(2,@() inferHyperspherePosterior(X(:,:,i),100000,false))),1:N);
% s(:,6) = arrayfun(@(i) undeal(3,@() estimateHypersphere(X(:,:,i))),1:N);

parfor i = 1:N
   stationarycounter(i,N)
   [~,~,tmp2(i)] = estimateHypersphere(X(:,:,i));
   tmp5(i)   = mean(undeal(2,@() inferHyperspherePosterior(X(:,:,i),1000,false)));
   tmp6(i)   = mean(undeal(2,@() inferHyperspherePosterior(X(:,:,i),100000,false)));
%   maxradii(i) = fminunc(@(m) maxRadiusGivenCenter(m,Xc(:,:,i)),mean(Xc(:,:,i)),opts);
end
s(:,2) = tmp2; s(:,5) = tmp5; s(:,6) = tmp6;
%s(:,6) = maxradii + maxradii*(n^-d);
ev = maxradii/target;

if PLOT,   figure(98);clf;plotEstimators(Xc,radii/target,s/target,colors);   end
if SAVE,   saveRadiusEsts(type);   end

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
legend('','expDistPerRad','','Normalized Std',...
       '','Joint ML','','MVUE','','MCMC 1k','','MCMC 100k',...
       'Location','EastOutside')
set(legend,'Box','off')
ylim([-10 0])
return
end

%% SAVE function
function saveRadiusEsts(distType,exptType)
   if ~exist('exptType','var'), exptType = []; end
   tickscale;logAxis(2);axesSeparate;
   set(99,'Name',['Radius_' distType exptType],'renderer','painters');
   printFig(99,[],'eps');
end

%% Objective function to optimize, with gradient
function [fval,grad] = maxRadiusGivenCenter(m,X)

grad = X - repmat(m(:)',[size(X,1) 1]);
fval = sqrt(sum(grad.^2,2));

% Apply max function
[fval,ix] = max(fval);
grad = -grad(ix,:)/fval;
return
end


%% stdPer
function v = stdPer(d)
expectedStds = [1.2733    1.0115    0.8796    0.8107    0.8384    0.8638    0.9579    1.0403    1.1938  1.4268    1.8384    2.4485];
v = interp1(2.^(1:12),expectedStds,d);
% d = log2(d);
% v = -0.1*d+1.2;
% v = 0.031*d^2 - 0.32*d + 1.6;
% v = 0.00049*d^4 - 0.013*d^3 + 0.12*d^2 - 0.51*d + 1.6;
% v = 0.00063*d^4 - 0.015*d^3 + 0.14*d^2 - 0.61*d + 1.8;
return
end


%% expDistPerRad function
function r = expDistPerRad(d)
% by interpolation among estimates precomputed in QT_expectedDistBetweenTwoPointsWithinHypersphere
expectedDists = [0.9067    1.1037    1.2402    1.3215    1.3660    1.3902    1.4018    1.4081    1.4112  1.4126    1.4134    1.4138];

if d<4097,   r = interp1(2.^(1:12),expectedDists,d);
else         r = sqrt(2);
end
return
end


%% DEMOS/DEBUG
h2s_radii(200,log2space(1,12,12),'Uniform',true);
h2s_radii(200,log2space(1,12,12),'Gaussian',true);
h2s_radii(200,log2space(1,12,12),'Hypercubic',true);

end



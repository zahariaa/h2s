function model = hypersphere2sphere(pointsOrDistMat,categories,estimator,dimLow)

% or should it be nball2ball?


%% preparations

% categories.vecs
% categories.labels
% categories.names

% criterion = {'MAP','minPercentileDeviationFromMarginalMedians', 'minRelativeDeviationFromMarginalMedians')
%
% MAP:      use the maximum a posteriori estimates as visualisation targets
%
% maxMarginalPosterior: use the maximum of the posterior marginal for each 
%           parameter as the visualisation target 
%
% minPercentileDeviationFromMarginalMedians (minPDMM): minimize the maximum
%           percentile-point deviation from the medians of the marginals of
%           the descriptive statistics.
%
% minRelativeDeviationFromMarginalMedians (minRDMM): minimize the maximum
%           relative deviation from the medians of the marginals of the
%           descriptive statistics.
%
% descriptive statistics:
%           (1) radii (sphere radii)
%           (2) dists (sphere-center dists)
%           (3) margins (sphere margins)
%
% infer post over these three summary statistics
%
% PAPER TITLE IDEAS
% Hypersphere to sphere: Visualizing the relationships among
% compact high-dimensional distributions
%
%
%
% TASKS
%
% write this function
%
% model comparisons using crossvalidation (hold out 10% of the points)
%       - hypersphere
%       - isotropic gaussian
%       - general gaussian
%       - optimal-kernel gaussian kernel smoother
%
% hs2s tests: comparison to MDS
%       - two touching equal-radius hyperspheres in 1-10 dimensions
%       - two concentric hyperspheres of different radius in 1-10 dimensions
%       - larger hypersphere enclosing smaller one, touching in one surface point
%       - two intersecting hypershperes
%       - 3 spheres (two intersecting) -> 2D circles
%       - 5 hyperspheres in 1-100 dimensions
%       - as previous, but illustrate effect of sample size changes


%% control variables
nPostSamples = 1000; % number of samples from the post over the parameters of each hypershpere
monitor = false;

%% preparations
[~,nCats] = size(categories.vectors);

%if ~exist('estimator','var'), estimator = 'MCMC'; end
if ~exist('estimator','var') || isempty(estimator), estimator = 'meanDist'; end
if ~exist('dimLow'   ,'var'),                       dimLow    =          3; end

%% create point data from distances if necessary

% under construction
% for now assuming points have been passed

points = pointsOrDistMat;

[~,nDim] = size(points);


%% infer posts over hypersphere parameters
switch estimator
    case 'MCMC'
        postSampleLocs = nan(nPostSamples,nDim,nCats);
        postSampleRadii = nan(nPostSamples,1,nCats);
        postSampleLogLikelihoods = nan(nPostSamples,1,nCats);

        for catI = 1:nCats
            [postSampleLocs(:,:,catI),postSampleRadii(:,:,catI),postSampleLogLikelihoods(:,:,catI)] =...
                inferHyperspherePosterior(points(logical(categories.vectors(:,catI)),:), nPostSamples, monitor);
        end

    case 'meanDist'
        estLocs = nan(nCats,nDim);
        estRadii = nan(1,nCats);
        
        nBootstrapSamples = 0;
        for catI = 1:nCats
            [estLocs(catI,:),locCI,estRadii(catI),radCI] = estimateHypersphere(points(logical(categories.vectors(:,catI)),:),nBootstrapSamples);
        end
end

%% compute posteriors for distances and margins
% margin = distance - r1 - r2

switch estimator
    case 'MCMC'
        nCatPairs = (nCats^2-nCats)/2;
        postSampleDists = nan(nPostSamples,nCatPairs);
        postSampleMargins = nan(nPostSamples,nCatPairs);

        for sampleI = 1: nPostSamples
            catPairI = 1;
            for cat1I = 1:nCats-1
                for cat2I = cat1I+1:nCats
                    postSampleDists(sampleI,catPairI) = sqrt(sum((postSampleLocs(sampleI,:,cat1I) - postSampleLocs(sampleI,:,cat2I)).^2,2));
                    postSampleMargins(sampleI,catPairI) = postSampleDists(sampleI,catPairI) - sum(postSampleRadii(sampleI,1,[cat1I cat2I]));
                    catPairI = catPairI+1;
                end
            end
        end

        % This could be sped up by using arrays instead of the for loops.
        % However, it would require a lot of memory (2*nCats times more).
    
    case 'meanDist'
        estDists = pdist(estLocs,'Euclidean');
        estMargins = nan(size(estDists));

        catPairI = 1;
        for cat1I = 1:nCats-1
            for cat2I = cat1I+1:nCats
                estMargins(catPairI) = estDists(1,catPairI) - sum(estRadii([cat1I cat2I]));
                catPairI = catPairI+1;
            end
        end
end

%% compute target point estimates to be visualized
if isequal(estimator,'MCMC')
    % use posterior medians as estimates
    estRadii   = median(postSampleRadii,1);   estRadii =   estRadii(:)';
    estDists   = median(postSampleDists,1);   estDists =   estDists(:)';
    estMargins = median(postSampleMargins,1); estMargins = estMargins(:)';
    disp(any2str('median radii: ',estRadii));
    disp(any2str('median distances: ',estDists));
    disp(any2str('median margins: ',estMargins));
end

%% initialize 3D model
mdsCenters_3d = mdscale(estDists,dimLow,'criterion','metricstress');
model.centers = mdsCenters_3d;
model.radii = estRadii;
model.categories = categories;
model.error = NaN;

%showSpheresModel(model, 200, '\bfinitial model based on MDS')


%% optimize 3D model
% TASK: further optimise the model to minimise various specific cost
% functions


%% return remaining error
model = measureError(model, estDists, estRadii,  estMargins);

%%%diagnostics(model, estDists, estRadii, estMargins);


%% measure error
function model = measureError(model,dataDists,dataRadii,dataMargins)

nCats = numel(model.radii);

model.dists = pdist(model.centers);
model.margins = squareform(model.dists) - ...
    ( repmat(model.radii,[nCats 1])+repmat(model.radii',[1 nCats]) )
model.margins(logical(eye(nCats)))=0;
model.margins = squareform(model.margins);

modelParams = [model.dists model.radii model.margins];
dataParams = [dataDists dataRadii dataMargins];
mx = max(modelParams);
errors = dataParams/mx - modelParams/mx;
model.error = max(abs(errors));



%% show summary statistics
function diagnostics(model,dataDists,dataRadii,dataMargins)

% compare model and data statistics
h=figure(505); set(h,'Color','w'); clf;
showSpheresModel(model, [505 2 2 1], {'\bfhs2s\rm',any2str('error <= ',ceil(model.error*100),'%')});

% summary stats for data and model as bar graphs
dataCol =  [0 0 0];
modelCol = [0.5 0.5 0.5];

subplot(2,3,4);
radiusBars = [dataRadii' model.radii'];
h=bar(radiusBars,'grouped','EdgeColor','none');
colormap([dataCol; modelCol]);

mx = max(model.radii);
radErrors = dataRadii/mx - model.radii/mx;
model.radiiError = max(abs(radErrors));

title({'\bfradii',['\rmerror: ',num2str(model.radiiError,2)]});
legend({'data','model'},'Location','SouthEast');

subplot(4,3,[8 9]);
distBars = [dataDists' model.dists'];
if numel(model.dists)==1, distBars = [distBars; [0 0]]; end
h=bar(distBars,'group','EdgeColor','none');
colormap([dataCol; modelCol]);

mx = max(model.dists);
distErrors = dataDists/mx - model.dists/mx;
model.distError = max(abs(distErrors));

title({'\bfdists',['\rmerror: ',num2str(model.distError,2)]});
%legend({'data','model'},'Location','SouthEast');
nanaxis([0.5 numel(model.dists)+0.5 nan nan]);

subplot(4,3,[11 12]); cla;
marginBars = [dataMargins' model.margins'];
if numel(model.margins)==1, marginBars = [marginBars; [0 0]]; end
h=bar(marginBars,'group','EdgeColor','none');
colormap([dataCol; modelCol]);

mx = max(model.margins);
marginErrors = dataMargins/mx - model.margins/mx;
model.marginError = max(abs(marginErrors));

title({'\bfmargins',['\rmerror: ',num2str(model.marginError,2)]});
%legend({'data','model'},'Location','SouthEast');
nanaxis([0.5 numel(model.margins)+0.5 nan nan]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colour-coded matrix comparison
% marginCols = colourScale([ 22 31 40; 255 100 100]/255,100);
% distRadCols = colourScale([240 240 240; 22 31 40]/255,100);
subplot(2,2,2);

marginCols = colourScale([1 0 0; 1 1 1; 0 0.5 1],99);
distRadCols = colourScale([0 0 0; 1 1 1],100);

nCats = numel(model.radii);
marginWidth = nCats^(1/2);
axis([0 nCats 0 (2+marginWidth)*nCats]); axis equal;
set(gca,'YDir','reverse');

maxDistRad = max([dataRadii(:);model.radii(:);dataDists(:);model.dists(:)]);
maxAbsMargin = max(abs([dataMargins(:);model.margins(:)]));


dists = squareform(dataDists);
radii = dataRadii;
margins = squareform(dataMargins);

shift = 0;
w = 0.9;
lw = 2;
for displayI = 1:2
    for cat1i = 1:nCats
        for cat2i = cat1i:nCats
            if cat1i==cat2i
                % colour patch for radius
                col = distRadCols(1+round(radii(cat1i)/maxDistRad*99),:);
                rectangle('Position',[cat1i-1,shift+cat2i-1,w,w],'Curvature',0.2,'FaceColor',col,'EdgeColor',model.categories.colors(cat1i,:),'LineWidth',lw);
            else
                % colour patch for dist
                if maxDistRad == 0
                    col = marginCols(round(size(distRadCols,1)/2),:);
                else
                    col = distRadCols(1+round(dists(cat1i,cat2i)/maxDistRad*99),:);
                end
                rectangle('Position',[cat1i-1,shift+cat2i-1,w,w],'Curvature',0.2,'FaceColor',col,'EdgeColor','none');

                % colour patch for margin
                if maxAbsMargin == 0
                    col = marginCols(round(size(marginCols,1)/2),:);
                else
                    col = marginCols(50+round(margins(cat1i,cat2i)/maxAbsMargin*49),:);
                end
                rectangle('Position',[cat2i-1,shift+cat1i-1,w,w],'Curvature',0.2,'FaceColor',col,'EdgeColor','none');
            end
        end % cat2i
    end % cat1i
    
    shift = nCats*(1+marginWidth);
    radii = model.radii;
    dists = squareform(model.dists);
    margins = squareform(model.margins);
    
end % displayI

% labels
txtMarg = 0.3;
text(0, 0-2*txtMarg, '\bfdata');
text(0, shift-2*txtMarg, '\bfmodel');

for catI = 1:nCats
    text(nCats+0.3, catI-1+0.5, model.categories.labels{catI});
    text(nCats+txtMarg, shift+catI-1+0.5, model.categories.labels{catI});
    text(catI-1+0.5, shift+nCats+txtMarg, model.categories.labels{catI},'Rotation',90,'HorizontalAlignment','right');      
end


% colour scales
hold on;
h = 1.0*txtMarg; % height of horizontal colourbars
s = 3*txtMarg; % vertical offset between label and colourbar
s2 = 5*txtMarg; % vertical offset between colourbars

y = shift+nCats+18*txtMarg; % vertical starting point
text(0, y, '\bfradii and dists \rm[min-max]');
image([0 nCats],[y+s y+s+h],reshape(distRadCols,[1 size(distRadCols)]));

y = y+s+h+s2;
text(0, y, '\bfmargins \rm[min-max]');
image([0 nCats],[y+s y+s+h],reshape(marginCols,[1 size(marginCols)]));

y = y+s+h+s2;
axis([0 nCats*1.5 -nCats-4 y]);
axis equal off;

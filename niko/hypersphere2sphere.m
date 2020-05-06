function [model,high] = hypersphere2sphere(pointsOrDistMat,categories,estimator,dimLow)

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
% monitor = false;

%% preparations
if ~exist('dimLow'   ,'var'),                       dimLow    =          3; end
if ~exist('estimator','var') || isempty(estimator), estimator = 'meanDist'; nPostSamples = 1;
end

high  = Hypersphere('estimate',pointsOrDistMat,categories,estimator,nPostSamples);
model = SetOfHyps(high).h2s(dimLow);

% model.show(200, '\bfinitial model based on MDS')
% diagnostics(model, high);
return



%% show summary statistics
function diagnostics(model,data)

% compare model and data statistics
h=figure(505); set(h,'Color','w'); clf;
showModel(model, [505 2 2 1], {'\bfhs2s\rm',any2str('error <= ',ceil(model.error*100),'%')});

% summary stats for data and model as bar graphs
dataCol =  [0 0 0];
modelCol = [0.5 0.5 0.5];

subplot(2,3,4);
radiusBars = [data.radii' model.radii'];
h=bar(radiusBars,'grouped','EdgeColor','none');
colormap([dataCol; modelCol]);

mx = max(model.radii);
radErrors = data.radii/mx - model.radii/mx;
radiiError = max(abs(radErrors));

title({'\bfradii',['\rmerror: ',num2str(radiiError,2)]});
legend({'data','model'},'Location','SouthEast');

subplot(4,3,[8 9]);
distBars = [data.dists' model.dists'];
if numel(model.dists)==1, distBars = [distBars; [0 0]]; end
h=bar(distBars,'group','EdgeColor','none');
colormap([dataCol; modelCol]);

mx = max(model.dists);
distErrors = data.dists/mx - model.dists/mx;
distError = max(abs(distErrors));

title({'\bfdists',['\rmerror: ',num2str(distError,2)]});
%legend({'data','model'},'Location','SouthEast');
nanaxis([0.5 numel(model.dists)+0.5 nan nan]);

subplot(4,3,[11 12]); cla;
marginBars = [data.margins' model.margins'];
if numel(model.margins)==1, marginBars = [marginBars; [0 0]]; end
h=bar(marginBars,'group','EdgeColor','none');
colormap([dataCol; modelCol]);

mx = max(model.margins);
marginErrors = data.margins/mx - model.margins/mx;
marginError = max(abs(marginErrors));

title({'\bfmargins',['\rmerror: ',num2str(marginError,2)]});
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

maxDistRad = max([data.radii(:);model.radii(:);data.dists(:);model.dists(:)]);
maxAbsMargin = max(abs([data.margins(:);model.margins(:)]));


dists = squareform(data.dists);
radii = data.radii;
margins = squareform(data.margins);

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

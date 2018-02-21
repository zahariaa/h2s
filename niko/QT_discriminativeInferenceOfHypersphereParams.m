% QT_discriminativeInferenceOfHypersphereParams

% PROBLEM: infer the radius of a hypersphere from points drawn uniformly
% within it.
%
% solution 1: bayesian inference by MCMC: inferHyperspherePosterior
% solution 2 (explored here): measure std (pooled across dimensions) and
% map approximately to the maximum posterior radius of the generating
% uniform-within-hypersphere distribution and its posterior variance
%
% OBSERVATIONS
% - draws from within a unit hypersphere multiplied by s are draws from a
% radius-s hypersphere.
% - std of s times the coords of a set of points equals s times the std of
% the coordinates.
% - simulations only needed for unit hypersphere centered at the origin.
% shifts and scalings of the generative parameters result in shifts and
% scalings of the summary stats (mean and standard deviation).
% - radius 0, implies std 0.
% - std is proportional to radius.
% - for >=0 points, the radius to std ratio depends only on the
% dimensionality of the space.
% - the expected distance to center of a point grows with
% dimensionality, but in a saturating fashion.
% - the std (sqrt(exp. dist. to center^2 / dimensionality) ) decays with
% the dimensionality (figs 10, 20, 30).
% - the radius can be accurately estimated as std (pooled across
% dimensions) times radius to std ratio (which is a function of the
% dimensionality) (fig 40).
% - the approximate width of the distribution of radius estimates (given
% the true radius is 1) is roughly linearly related to the reciprocal of
% the square root of the number of scalars in the data matrix (figs 50, 
% 60).
% - for the locations parameter p(loc|mean) = p(mean|loc), because the
% joint density of the data mean and the data-generating model's location
% is a diagonal convolved with a 1-d gaussian.
% - for the radius, it may be more complicated...

% TASKS
% estimate p(radius|std), rather than p(std|radius) or p(estimated
% radius|radius). 
% p(estimated radius | radius) = p(estimated radius | radius = 1) stretched
% by factor radius along the estimated radius dimension.


%% control variables
nsDims = [1 2 5 10 40 200 1000]
nsPoints = [2 3 5 10 20 40 100 1000]
nSims = 1000;
radius = 1;


%% simulate uniform draws within unit spheres in different dimensions
stds=nan(numel(nsPoints),numel(nsDims),nSims);

h=figure(10); set(h,'Color','w'); clf;

for nPointsI = 1:numel(nsPoints)
    nPoints = nsPoints(nPointsI)

    for nDimsI = 1: numel(nsDims)
        nDims = nsDims(nDimsI)

        for simI = 1:nSims
            points = randsphere(nPoints,nDims,radius);
            %points = randn(2*nPointsPerCat,nDim);

            stds(nPointsI,nDimsI,simI) = sqrt(mean(std(points,0,1).^2,2));
        end

        subplot(numel(nsPoints),numel(nsDims),(nPointsI-1)*numel(nsDims)+nDimsI);
        hist(stds(nPointsI,nDimsI,:)); title(any2str(nPoints,' points in ',nDims,' dims'));
        xlabel('dim-pooled std');
        
    end
end


%% show summary stats
meanStds = mean(stds,3);

h=figure(20); set(h,'Color','w'); clf;
imagesc(meanStds);
colormap('bone');
colorbar;

set(gca,'YTick',1:numel(nsPoints),'YTickLabel',nsPoints); ylabel('number of points');
set(gca,'XTick',1:numel(nsDims),'XTickLabel',nsDims); xlabel('number of dimensions');
title('mean across simulations of dim-pooled std');


h=figure(30); set(h,'Color','w'); clf;
plot(nsDims,meanStds(end,:),'k');
xlabel('number of dimensions');
ylabel('mean of simulated dim-pooled stds');


h=figure(40); set(h,'Color','w'); clf;

radToStdRatio_dimI = 1./meanStds(end,:);
estimatedRadii = stds.*repmat(radToStdRatio_dimI,[numel(nsPoints) 1 nSims]);

for nPointsI = 1:numel(nsPoints)
    nPoints = nsPoints(nPointsI)

    for nDimsI = 1: numel(nsDims)
        nDims = nsDims(nDimsI)

        subplot(numel(nsPoints),numel(nsDims),(nPointsI-1)*numel(nsDims)+nDimsI);
        hist(estimatedRadii(nPointsI,nDimsI,:)); title(any2str(nPoints,' points in ',nDims,' dims'));
        xlabel('estimated radius');
        
    end
end


h=figure(50); set(h,'Color','w'); clf;
estRadiiStds = std(estimatedRadii,0,3);
imagesc(estRadiiStds);
colormap('bone');
colorbar;

set(gca,'YTick',1:numel(nsPoints),'YTickLabel',nsPoints); ylabel('number of points');
set(gca,'XTick',1:numel(nsDims),'XTickLabel',nsDims); xlabel('number of dimensions');
title('std of radius estimate (true radius = 1)');


h=figure(60); set(h,'Color','w'); clf;

nDataScalars = nsPoints'*nsDims;

subplot(2,2,1); plot(nDataScalars(:),estRadiiStds(:),'k.');
xlabel('number of data scalars'); ylabel('std of radius estimate')
axis equal square;
subplot(2,2,2); plot(sqrt(nDataScalars(:)),estRadiiStds(:),'k.');
xlabel('sqrt of number of data scalars'); ylabel('std of radius estimate')
axis equal square;
subplot(2,2,3); plot(1./sqrt(nDataScalars(:)),estRadiiStds(:),'k.');
xlabel('1/sqrt of number of data scalars'); ylabel('std of radius estimate')
axis equal square;

nDF = nsPoints'*(nsDims-1);

subplot(2,2,4); plot(1./sqrt(nDF(:)),estRadiiStds(:),'k.');
xlabel('1/sqrt of number of degrees of freedom'); ylabel('std of radius estimate')
axis equal square;






function [posteriorSamplesLoc,posteriorSamplesRad,logLikelihoods] = inferHyperspherePosterior(points, nSamples, monitor,VERBOSE)

% FUNCTION
% Given data points in argument 'points', this function infers the
% posterior over the parameters (center, radius) of a hypersphere-uniform
% hypothesis space of data-generative models. The prior is an (improper)
% flat infinite distribution over all locations and nonnegative radii. The
% posterior is inferred using the Metropolis-Hastings algorithm for
% Markov-Chain Monte Carlo sampling and returned in the form of nSamples
% samples of the parameters. 
%
% ARGUMENTS
% points    number of data points by number of dimensions matrix of data
%           points
%
% nSamples  number of samples to draw (inlcuding burn in)
%
% RETURNS
% posteriorSamplesLoc location samples (nSamples by number of dimensions)
%           from the inferred posterior
%
% posteriorSamplesRad radius samples (nSamples by 1) from the inferred
%           posterior
%
% Note: In topology, hypersphere refers to the surface, n-ball to the
% volume. Hypersphere-uniform, here, refers to the uniform distribution
% within a hypersphere. This could also be called the uniform n-ball 
% distribution. The visualisation method this function is for could be 
% called hypersphere2sphere of n-ball2ball.


%% preparations
if ~exist('nSamples','var'), nSamples = 1000; end
if ~exist('monitor','var'), monitor = true; end
if ~exist('VERBOSE','var') || isempty(VERBOSE), VERBOSE = false; end

[nPoints,nDim] = size(points);

nSamples = max(400,nSamples); % forbid less than 400 samples
posteriorSamplesLoc = nan(nSamples,nDim);
posteriorSamplesRad = nan(nSamples,1);
logLikelihoods = nan(nSamples,1);


%% initialise with all-inclusive hypershpere
cLoc = mean(points);
dists2center = sqrt(sum((points - repmat(cLoc,[nPoints 1])).^2,2));

cRad = 1.5*max(dists2center);

if cRad>0
    proposalStdLoc = cRad/1;
    proposalStdRad = cRad/2;
else
    cRad = 1;
    proposalStdLoc = 100;
    proposalStdRad = 10;
end

cLogL = logLikelihood(cLoc, cRad, points);

sampleI = 1;
nAccepted = 0;


%% sample from the posterior
burnInAndStepSizeTuningMode = true;
nStableCycles = 0;

while true

    posteriorSamplesLoc(sampleI,:) = cLoc;
    posteriorSamplesRad(sampleI) = cRad;
    logLikelihoods(sampleI) = cLogL;
    sampleI = sampleI + 1;
    
    if burnInAndStepSizeTuningMode
        % burn-in and step-size tuning period

        if sampleI==400
            acceptRate = nAccepted/sampleI;
            
            f=1+rand/2;
            if acceptRate<0.2
                nStableCycles = 0;
                proposalStdLoc = proposalStdLoc/f;
                proposalStdRad = proposalStdRad/f;
                if VERBOSE, disp(any2str('MCMC: accept rate = ',nAccepted/sampleI,'. Reducing step size.')); end
            elseif 0.5<acceptRate
                nStableCycles = 0;
                proposalStdLoc = proposalStdLoc*f;
                proposalStdRad = proposalStdRad*f;
                if VERBOSE, disp(any2str('MCMC: accept rate = ',nAccepted/sampleI,'. Increasing step size.')); end
            else
                nStableCycles = nStableCycles + 1;
                if VERBOSE, disp(any2str('MCMC: step-size tuning ',nStableCycles,'/5 successful.')); end
            end

            if nStableCycles == 5
                % end burn-in and tuning mode
                burnInAndStepSizeTuningMode = false;
            end
            nAccepted = 0;
            sampleI = 1;
        end
    else
        % freeze step size         
        if sampleI > nSamples, break; end
        if mod(sampleI,10000)==0 && VERBOSE,
            disp([num2str(round(sampleI/nSamples*100)),'% done']);
        end
    end
        
    % propose step
    candLoc = cLoc + proposalStdLoc * randn(1,nDim);
    candRad = cRad + proposalStdRad * randn(1,1);
    candLogL = logLikelihood(candLoc, candRad, points);
    
    if rand < exp(candLogL - cLogL)
        cLoc = candLoc;
        cRad = candRad;
        cLogL = candLogL;
        nAccepted = nAccepted+1;
        % Iain Murray: "If you reject at least about half your samples,
        % this thing is probably working well."
    end
end     


%% visualise results
if ~monitor, return; end

% show joint posterior for radius and location (dimension 1)
h=figure(100); clf; set(h,'Color','w');
if VERBOSE, disp(['accepted ',num2str(nAccepted),' proposals']); end

subplot(2,1,1); 
plot(posteriorSamplesLoc(:,1),posteriorSamplesRad,'k.','MarkerSize',1); hold on;
plot(posteriorSamplesLoc(1,1),posteriorSamplesRad(1),'b.','MarkerSize',20);
plot(mean(posteriorSamplesLoc(:,1)),mean(posteriorSamplesRad),'r.','MarkerSize',30);
title({['\bf',num2str(round(nSamples/1000)),'K samples from the posterior'],['\rm(burn-in excluded; accept rate ',num2str(nAccepted/nSamples,2),')']});
xlabel('location (dimension 1)');
ylabel('radius'); axis equal;

subplot(2,1,2); 
plot(posteriorSamplesLoc(floor(end/2):end,1),posteriorSamplesRad(floor(end/2):end),'k.','MarkerSize',1); hold on;
plot(mean(posteriorSamplesLoc(floor(end/2):end,1)),mean(posteriorSamplesRad(floor(end/2):end)),'r.','MarkerSize',30);
title({'\bfsamples from the posterior','\rm(second half)'});
xlabel('location (dimension 1)');
ylabel('radius'); axis equal;


% show marginals posteriors (histograms) for radius and location (dimension 1) 
h=figure(110); clf; set(h,'Color','w');
nBins=50;

subplot(2,1,1); 
hist(posteriorSamplesLoc(floor(end/2):end,1),nBins);
title({'\bfhistogram of location samples from the posterior','\rm(second half)'});
xlabel('location (dimension 1)');

subplot(2,1,2); 
hist(posteriorSamplesRad(floor(end/2):end,1),nBins);
title({'\bfhistogram of radius samples from the posterior','\rm(second half)'});
xlabel('radius');


% draw interval samples from the posterior (for dimension 1) 
h=figure(120); clf; set(h,'Color','w');

plot(points(:,1),0,'k.','MarkerSize',20); hold on;

subSample = round(linspace(nSamples/2,nSamples,500));
for subSampleI = 1:numel(subSample)
    sampleI = subSample(subSampleI);
    loc = posteriorSamplesLoc(sampleI);
    rad = posteriorSamplesRad(sampleI);
    plot([loc-rad loc+rad],20+[subSampleI subSampleI],'k');
end
    
title({'\bfdata and interval samples from the posterior','\rm(second half)'});
xlabel('location (dimension 1)');
set(gca,'YTick',[]);


% probability of new point being inside, given the data
% (plotted for dimension 1, other dimensions assumed to match location)
h=figure(130); clf; set(h,'Color','w');

newPoints_dim1 = -50:0.1:50;
leftBounds = posteriorSamplesLoc(:,1)-posteriorSamplesRad;
rightBounds = posteriorSamplesLoc(:,1)+posteriorSamplesRad;

inInterval = (repmat(leftBounds,[1 numel(newPoints_dim1)]) < repmat(newPoints_dim1,[nSamples 1]) )...
                                                         & ( repmat(newPoints_dim1,[nSamples 1]) < repmat(rightBounds,[1 numel(newPoints_dim1)]) );                                                    
plot(points(:,1),0,'k.','MarkerSize',20); hold on;
plot(newPoints_dim1,mean(double(inInterval),1),'k');
                                                    
figure(100);


%% log likelihood function
function logL = logLikelihood(center, radius, points)

nDim = numel(center);
nPoints = size(points,1);

dists2center = sqrt(sum((points - repmat(center,[nPoints 1])).^2,2));
if any(dists2center>radius)
    logL = -inf;
else
    %vol = nBallVolume(nDim,radius);
    % The volume of the n-ball equals a nDim-dependent factor *
    % radius^nDim. Since the nDim-dependent factor is the same for all
    % hypotheses considered, it just scales all likelihoods by the same
    % factor and we can neglect it. Computationally, we need to neglect it
    % because the volume will approach 0 for nDim>25. 
    volFac = radius^nDim;
    probDensity = 1/volFac;

    %logL = log(probDensity^nPoints);
    logL = nPoints*log(probDensity);
end



function [patterns,errorSSQ]=rdm2patternEnsemble(rdm,nDim,distanceMeasure)


%% preparations
rdm=squareRDMs(unwrapRDMs(rdm));
[~,nCond]=size(rdm);
if ~exist('distanceMeasure','var'), distanceMeasure='euclidean'; end

stepSize=10;


%% try MDS
% patterns_mdscale=mdscale(rdm,nDim,'Criterion','metricstress');
% showRDMs(concatRDMs(rdm,pdist(patterns_mdscale,distanceMeasure)));

% conclusion:
% mdscale does not work on the easy case of finding an arrangement in a
% high-dimensional space


%% simple gradient descent

% initialise with random patterns
patterns=randn(nCond,nDim)+100;

% scale patterns such that the root mean square distance matches that of
% the target RDM
patterns = patterns./ sqrt(sum(pdist(patterns,distanceMeasure).^2)) .* sqrt(sum(vectorizeRDMs(rdm).^2));
% sqrt(sum(pdist(patterns,distanceMeasure).^2))
% sqrt(sum(vectorizeRDMs(rdm).^2))

errorSSQ=inf;




% drive distances to match rdm
iterationCount=0;
nSmallProgressSteps = 0;

while true;
    
    dists=squareRDMs(pdist(patterns,distanceMeasure));
    
    distDiffs=dists-rdm;
    prevErrorSSQ=errorSSQ;
    errorSSQ=nansum(distDiffs(:).^2);
    disp(any2str('errorSSQ=',errorSSQ,', stepSize=',stepSize,', nSmallProgressSteps=',nSmallProgressSteps));
    
    if errorSSQ<1e-1 || 100<nSmallProgressSteps; 
        break; 
    end
    
    if errorSSQ>prevErrorSSQ
        patterns=prevPatterns;
        stepSize = stepSize/2;
    end

    if abs(errorSSQ-prevErrorSSQ)<0.001
        patterns=prevPatterns;
        nSmallProgressSteps = nSmallProgressSteps+1;
    end

    if mod(iterationCount,100)==0 && false,
        showRDMs(cat(3,rdm,dists),100); 
        subplot(2,2,1); title('target RDM');
        subplot(2,2,2); title('RDM of re-embedded points');
        subplot(2,2,4); cla; hold on;
        plot([zeros(100,1) 100*randn(100,1)]', [zeros(100,1) 100*randn(100,1)]','Color',[0.7 0.7 0.7]);
        plot(rdm(:),dists(:),'k.');
        mx=max([rdm(:);dists(:)])
        axis square; axis([0 mx 0 mx]);
        xlabel('target RDM'); ylabel('pattern RDM');
        errorSSQ
        stepSize
    end
    
    distFacs=rdm./dists; 
    
    diffVecs=repmat(reshape(patterns,[nCond 1 nDim]),[1 nCond 1])-repmat(reshape(patterns,[1 nCond nDim]),[nCond 1 1]);
    stepVecs=stepSize*diffVecs.*repmat(distFacs-1,[1 1 nDim]);
    
    prevPatterns=patterns;
    patterns=patterns-reshape(nanmean(stepVecs,1),[nCond nDim]);
    
    %     if strcmp(distanceMeasure,'correlation')
    %         patterns=patterns-repmat(mean(patterns,2),[1 nDim]); % remove the overall activation component
    %         patterns=patterns./repmat(std(patterns,1,2),[1 nDim]); % normalise the variance across simulated voxels for each pattern
    %     end
    
    iterationCount=iterationCount+1;
    
end

showRDMs(cat(3,rdm,dists),100);
subplot(2,2,1); title('target RDM');
subplot(2,2,2); title('RDM of re-embedded points');
subplot(2,2,4); cla; hold on;
plot([zeros(100,1) 100*randn(100,1)]', [zeros(100,1) 100*randn(100,1)]','Color',[0.7 0.7 0.7]);
plot(rdm(:),dists(:),'k.');
mx=max([rdm(:);dists(:)])
axis square; axis([0 mx 0 mx]);
xlabel('target RDM'); ylabel('pattern RDM');




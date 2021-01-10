% QT_rdm2patternEnsemble

load(['92imageData',filesep,'92_brainRDMs.mat']);
rdm_hIT=mean(unwrapRDMs(reshape(RDMs,[1 8])),3);

showRDMs(rdm_hIT);

%distanceMeasure='euclidean';
distanceMeasure='correlation';
nDim=100;
[patterns,errorSSQ]=rdm2patternEnsemble(rdm_hIT,nDim,distanceMeasure);


dists=squareRDMs(pdist(patterns,distanceMeasure));

patterns2=patterns./repmat(std(patterns,1,2),[1 nDim]);
patterns2=patterns2-1.2*min(patterns2(:));
dists2=squareRDMs(pdist(patterns2,distanceMeasure));

showRDMs(concatRDMs(rdm_hIT,dists,dists2),200);

figurePanel(210);
subplot(2,1,1); imagesc(patterns); colorbar;
subplot(2,1,2); imagesc(patterns2); colorbar;


save('reconstructedPats92','patterns','patterns2');


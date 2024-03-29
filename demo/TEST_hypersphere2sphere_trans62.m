% TEST_hypersphere2sphere_transRSA62


%% control variables
DEBUG  = false;
dimLow = 3;


%% define categories
controlledHumanFaces =  false(62,1); controlledHumanFaces(1:14)=true;
humanFaces =  false(62,1); humanFaces(15:20)=true;
animalFaces =  false(62,1); animalFaces(21:26)=true;
humanBodies =  false(62,1); humanBodies(27:32)=true;
animalBodies =  false(62,1); animalBodies(33:38)=true;
artificialObjects =  false(62,1); artificialObjects(39:44)=true;
naturalObjects =  false(62,1); naturalObjects(45:50)=true;
artificialScenes =  false(62,1); artificialScenes(51:56)=true;
naturalScenes =  false(62,1); naturalScenes(57:62)=true;



%% roll all subsets and supersets into one
allcategories.vectors = [humanFaces | animalFaces | humanBodies | animalBodies ... % animate
                         artificialObjects | naturalObjects | artificialScenes | naturalScenes ... % inanimate
                         humanFaces | animalFaces humanBodies | animalBodies ... % faces & bodies
                         artificialObjects | naturalObjects artificialScenes | naturalScenes ... % objects & scenes
                         humanFaces animalFaces humanBodies animalBodies artificialObjects naturalObjects artificialScenes naturalScenes]; % all
allcategories.labels  = {'animate' 'inanimate' 'faces' 'bodies' 'objects' 'scenes' ...
                         'humanFaces' 'animalFaces' 'humanBodies' 'animalBodies' 'artificialObjects' 'naturalObjects' 'artificialScenes' 'naturalScenes'};
allcategories.colors = [215   148   196% animate   = magenta
                        135   184    98% inanimate = yellow
                        246   136   159% faces     = red
                        100   180   210% bodies    = blue
                         88   189   139% objects   = green
                        223   157    79% scenes    = orange
                        ]/255; 
ac = Categories(allcategories);

%% define ROIs
roi.name = {'V1','V2','V3','LOC','OFA','FFA','PPA','aIT'};
roi.hemisphere = {'left','right'};
roi.size = [20 40 80 120 160 320];

% choose fixed hemisphere and ROI size (for the moment)
roiSizeI=6;
hemisphereI = 2;


%% make hypersphere2sphere visualization (and MDS) for RDMs from walther et al. (in preparation)
load('crossvalidatedRDMs_leaveOneRunOut_concatYtest_alldata.mat');
avgRDMs=nan(1,1891,2);

selectedROIs = [1 4 2 6 3 7];
nROIs        = numel(selectedROIs);
fh = newfigure([nROIs/2 3*2],'trans62');
figrowI = 0;

for roiI = selectedROIs
    skipRoi = false;
    rdms = [];
    for sessI=1:2
        for subjI=1:24
            rdmStruct = getrow(D, D.subj==subjI & D.sess==sessI & D.hemi==hemisphereI & D.roi==roiI & D.roisize==4);
            try
                rdms = cat(3,rdms,rdmStruct.rdm);
            catch
                skipRoi=true;
            end
        end
    end % sessI
    
    if skipRoi, continue; end
    
    avgRdm = mean(rdms,3); % averaged across sessions and subjects
    avgRdm_nonNeg = avgRdm;
    avgRdm_nonNeg(avgRdm<0) = 0;
    
    % DEBUG: rand_rdm = pdist(randn(30,500),'Euclidean');
   if DEBUG
    showRDMs(cat(3,avgRdm,avgRdm_nonNeg));
    subplot(2,2,1); title(any2str(roi.hemisphere{hemisphereI},' ',roi.name{roiI},' (',roi.size(roiSizeI),' vox)'));
    subplot(2,2,2); title('non-negative version');
    subplot(2,2,4); plot(avgRdm(:),avgRdm_nonNeg(:),'k.');
    xlabel('rdm');
    ylabel('rdm-nonNeg');
   end
    % hs2s & MDS
    titleStr = any2str(roi.hemisphere{hemisphereI},' \bf',roi.name{roiI},'\rm (',roi.size(roiSizeI),' vox)');
    figrowI = figrowI + 1;
    roirow = ceil(figrowI/2);
    roicol = mod(figrowI-1,2)*3;
    roiplt = (roirow-1)*2*3 + roicol + 1;
    
    % use embedded points
    % re-embed as points in 100-dimensional space to reflect the representational geometry
    [patterns,errorSSQ] = rdm2patternEnsemble(avgRdm_nonNeg,100,'Euclidean');
    % hyp(roiI) = HS2SandMDS(patterns,ac.select(3:6),fh.a.h(figrowI*2+(-1:0)),...
    %                        titleStr,dimLow,100);
    hyp(roiI) = Hypersphere.estimate(patterns,ac.select(3:6),5000).meanAndMerge;
    hyp(roiI).h2s(dimLow).show(fh.a.h(roiplt))
    hyp(roiI).showSig(fh.a.h(roiplt+(1:2)))

    % % direct to distance-based estimator
    % hyp(roiI) = HS2SandMDS(avgRdm_nonNeg,ac.select(3:6),fh.a.h(figrowI*2+(-1:0)),...
    %                        titleStr,dimLow,'distance',1000);
    drawnow;
   % keyboard
   % figure(2);subplot(3,2,figrowI);hyp(roiI).show;
end % roiI


[fh,f] = newfigure([nROIs/2 2*2],'trans62_stats',fh);
for i = 1:nROIs
    % hyp(selectedROIs(i)).showValues(fh.a(f).h(i*2-1));
    % hyp(selectedROIs(i)).showSig(fh.a(f).h(i*2),'diff');
    hyp(selectedROIs(i)).showSig(fh.a(f).h(i*2+(-1:0)));
end

figure;plotErrorPatch(1:4,hyp(7).radii,hyp(7).ci.radii)
figure;hist(reshape([hyp(7).ci.bootstraps.radii],4,[])')

%% Finish up comparison figure
papsz = [4 3];
set(fh.f,'PaperUnits','inches','PaperSize',papsz,'PaperPosition',[0.01*papsz papsz],'PaperPositionMode','manual');
%set(findobj(fh.f,'type','line'),'MarkerSize',3); % make dots slightly smaller
for f = 1:2; subplotResize(fh.f(f),0.01,0.01); end
set(fh.a(1).h(16)   ,'CameraViewAngle',5  );
set(fh.a(1).h([4 9]),'CameraViewAngle',3.5);
set(fh.f(2),'Renderer','painters');   printFig;
delete(fh.a(1).h(setdiff(1:nROIs*3,1:3:nROIs*3)));
set(fh.f(1),'Renderer','openGL');     printFig(fh.f(1),[],'png',1200);




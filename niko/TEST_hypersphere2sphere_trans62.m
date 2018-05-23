% TEST_hypersphere2sphere_transRSA62


%% control variables
DEBUG = false;



%% define categories
controlledHumanFaces =  zeros(62,1); controlledHumanFaces(1:14)=1;
humanFaces =  zeros(62,1); humanFaces(15:20)=1;
animalFaces =  zeros(62,1); animalFaces(21:26)=1;
humanBodies =  zeros(62,1); humanBodies(27:32)=1;
animalBodies =  zeros(62,1); animalBodies(33:38)=1;
artificialObjects =  zeros(62,1); artificialObjects(39:44)=1;
naturalObjects =  zeros(62,1); naturalObjects(45:50)=1;
artificialScenes =  zeros(62,1); artificialScenes(51:56)=1;
naturalScenes =  zeros(62,1); naturalScenes(57:62)=1;

categoryVectors=[controlledHumanFaces+humanFaces+animalFaces humanBodies+animalBodies artificialObjects+naturalObjects artificialScenes+naturalScenes];

categories.vectors =[humanFaces+animalFaces humanBodies+animalBodies artificialObjects+naturalObjects artificialScenes+naturalScenes];
categories.labels = {'faces','bodies','objects','scenes'};
categories.colors = [0.8 0 0
                     0.2 0.8 0
                     0   0.5 0.8
                     0.5 0.5 0.5];

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
    % re-embed as points in high-dimensional space to reflect the representational geometry
    nDim = 100;
    distanceMeasure = 'Euclidean';
    [patterns,errorSSQ] = rdm2patternEnsemble(avgRdm_nonNeg,nDim,distanceMeasure);
    
    % hs2s & MDS
    figPanelSpec = [150 4 8 4*(roiI-1)+1];
    titleStr = any2str(roi.hemisphere{hemisphereI},' \bf',roi.name{roiI},'\rm (',roi.size(roiSizeI),' vox)');
    psOutput='';
    pdfOutput='';
    HS2SandMDS(patterns,categories,figPanelSpec,titleStr,3);%,psOutput,pdfOutput)

    %model = multidimensionalCategoryScaling(avgRDMs(1,:,1),categoryVectors,catLabels,[],'ellipsoid',3);
    
    
end % roiI
          


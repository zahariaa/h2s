% TEST_uniformWithinHypersphereParamEstimators

%% control variables
clear;

nsPoints = [4:10 15:5:50 75:25:200];
%nsPoints = [40 80 160 320 640];
%nsPoints = 10;
nsPoints = 4:20;

%prior
nsDim=1:4;
nsDim=[1 2 3 6 15 40 100];
nsDim=[1 2 3 6 10 20 40 500 1000];
nsDim=[1 2 3 10 40 80];
nsDim = 100;

% note that sampling radii or locations would be redundant:
% same situation shifted for any location and scaled for any radius.


nSim = 100;

lw=2;


%% simulated densities
for nDimI = 1:numel(nsDim)
    nDim = nsDim(nDimI);
    
    disp('_____________________');
    nDim
    
    for nI = 1:numel(nsPoints)
        n = nsPoints(nI);
        
        %points=randn(n,nDim);
        points=randsphere(n,nDim,1);
        
        monitor = false;
        nSamples = 5000;
        [posteriorSamplesLoc posteriorSamplesRad logLikelihoods] = inferHyperspherePosterior(points, nSamples, monitor);

        MCMCerror_loc(nDimI,nI) = mean(median(posteriorSamplesLoc,1).^2); % MSE
        MCMCerror_rad(nDimI,nI) = (median(posteriorSamplesRad,1)-1)^2;    % MSE
        
        nBootstrapSamples = 0;
        [loc locCI rad radCI] = estimateHypersphere(points,nBootstrapSamples);
        
        fastEstError_loc(nDimI,nI) = mean(loc.^2); % MSE
        fastEstError_rad(nDimI,nI) = (rad-1)^2; % MSE
        
        i=i+1;
    end
end


%% plot MSE against number of points for each number of dimensions
h=figure(10); clf; set(h,'Color','w');
cols=colorScale([.85 .85 .85; 0 0 0],numel(nsDim));
clear legendStrings;

for nDimI = 1:numel(nsDim)
    nDim = nsDim(nDimI);
    
    subplot(2,2,1); hold on;
    plot(nsPoints,fastEstError_loc(nDimI,:),'Color',cols(nDimI,:), 'LineWidth',lw);

    subplot(2,2,2); hold on;
    plot(nsPoints,fastEstError_rad(nDimI,:),'Color',cols(nDimI,:), 'LineWidth',lw);
    legendStrings{nDimI}=[num2str(nDim),' dim.'];

    subplot(2,2,3); hold on;
    plot(nsPoints, MCMCerror_loc(nDimI,:), 'Color',cols(nDimI,:), 'LineWidth',lw);

    subplot(2,2,4); hold on;
    plot(nsPoints, MCMCerror_rad(nDimI,:), 'Color',cols(nDimI,:), 'LineWidth',lw);
    legendStrings{nDimI}=[num2str(nDim),' dim.'];
end

subplot(2,2,1);
xlabel('number of points'); ylabel('MSE');
title({'\bflocation','\rmmean squared error(# dimensions)'});
legend(legendStrings,'Location','NorthEast');

subplot(2,2,2);
xlabel('number of points'); ylabel('MSE');
title({'\bfradius','\rmmean squared error(# dimensions)'});

subplot(2,2,3);
xlabel('number of points'); ylabel('MSE');
title({'\bfMCMC: location','\rmmean squared error(# dimensions)'});
legend(legendStrings,'Location','NorthEast');

subplot(2,2,4);
xlabel('number of points'); ylabel('MSE');
title({'\bfMCMC: radius','\rmmean squared error(# dimensions)'});


%% histogram of MSE
h=figure(20); clf; set(h,'Color','w');
nBins = 10;

subplot(2,1,1);
hist([MCMCerror_loc(:) fastEstError_loc(:)],nBins);
legend({'MCMC','fast est.'});
xlabel('location mean squared error');
ylabel('frequency');

subplot(2,1,2);
hist([MCMCerror_rad(:) fastEstError_rad(:)],nBins);
legend({'MCMC','fast est.'});
xlabel('radius mean squared error');
ylabel('frequency');




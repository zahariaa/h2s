% QT_inferHyperspherePosterior

nSamples = 10000;

%% 1-d scenario
%points = 0;
%points = [-1 1]';
points = linspace(-2,2,5)';

%points = [-9 -8 -7 9]';

[posteriorSamplesLoc posteriorSamplesRad logLikelihoods] = inferHyperspherePosterior(points, nSamples);


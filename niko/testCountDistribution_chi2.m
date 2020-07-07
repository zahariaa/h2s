function [p,chi2]=testCountDistribution_chi2(counts,probabilities)
% USAGE
%       [p,chi2]=testCountDistribution_chi2(counts[,probabilities])
%
% FUNCTION
%       tests the null hypothesis that the mutually exclusive
%       event types counted in argument counts have the probabilities (summing to
%       one) defined in argument probabilities. if the second argument is
%       omitted, the null hypothesis of equal probabilities is tested.
%
% CAVEAT
%       "The approximation to the chi-square distribution breaks down if
%       expected frequencies are too low. It will normally be acceptable so
%       long as no more than 10% of the events have expected frequencies
%       below 5. Where there is only 1 degree of freedom, the approximation
%       is not reliable if expected frequencies are below 10."
%       (Wikipedia, July 2008)

nPossibleOutcomes=numel(counts);
df=nPossibleOutcomes-1; 

nEvents=sum(counts);

if ~exist('probabilities','var')
    probabilities=ones(1,nPossibleOutcomes)/nPossibleOutcomes;
end

expectedCounts=probabilities*nEvents;

chi2=sum((counts-expectedCounts).^2./expectedCounts);

p=1-cdf('chi2',chi2,df);



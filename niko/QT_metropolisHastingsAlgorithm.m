% QT_metropolisHastingsAlgorithm

%% control variables

domain = [1 2];
pdf = [2 1];

domain = [1 2 3 4 5];
pdf = [1 2 3 4 5];

domain = -3:0.1:3;
pdf = exp(-domain.^2);


nSamples = 100000;


%% markov-chain monte carlo by metropolis-hastings algorithm

acceptedSamples = nan(nSamples,1);
% Iain Murray: "If you are rejecting at least half your samples, this thing
% is probably working well."

currLocI = 1;

sampleI = 1;

while true
    acceptedSamples(sampleI) = currLocI;
    sampleI = sampleI + 1;
    if sampleI > nSamples, break; end
        
    cProbabilityDensity = pdf(currLocI);
    
    % propose step
    candLocI = ceil(rand*(numel(domain)-1));
    if candLocI>=currLocI, candLocI = candLocI+1; end
    
    if rand < pdf(candLocI)/pdf(currLocI)
        currLocI = candLocI;
    end
end     


%% visualise results
h=figure(100); 
subplot(2,1,1); bar(domain,pdf); title('true pdf');
subplot(2,1,2); hist(domain(acceptedSamples),nSamples/50); title('histogram of MCMC samples');



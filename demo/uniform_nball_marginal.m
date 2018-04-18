function uniform_nball_marginal(n,d,nbins)
% uniform_nball_marginal(n,d,nbins=100)
% Plots histogram of a uniform d-dimensional ball with n samples
% nbins: number of bins in histogram
% 
% 2018-04-18 AZ Created

if ~exist('nbins','var') || isempty(nbins),   nbins = 100;   end

%% Sample
X = randnball(n,d);
%% Analytic function
x = linspace(-1,1,nbins);
y = (1+log10(d))*(n/nbins)*sqrt(1-x.^2).^(d-1);

%% Plot
hold on;
hist(X(:,1),nbins);
plot(x,y,'r-')


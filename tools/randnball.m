function X = randnball(n,d,r)
% X = randnball(n,d,r=1)
% Generates n points uniformly distributed in d-dimensional 
% ball with radius r using Marsaglia (1972) method
%
% Marsaglia, G. (1972). "Choosing a Point from the Surface of a Sphere".
% Annals of Mathematical Statistics. 43 (2): 645–646.
% doi:10.1214/aoms/1177692644
%
% 2018-03-16 Andrew Zaharia

if ~exist('r','var') || isempty(r),   r = 1;   end

% Initialize
X = randn(n,d);
U = rand( n,1).^(1/d);
N = zeros(n,1);        % To be filled with norms of rows of X,

% ... to transform X onto unit (n-1)-sphere,
for i = 1:n
   N(i,:) = sqrt(X(i,:)*X(i,:)');
end

% ... and distribute uniformly from (n-1)-sphere to n-ball, radius r.
X = r*repmat(U./N,[1 d]).*X;

return

%% Debug
X=randnball(1000,2);figure;hold on;plot(X(:,1),X(:,2),'ko');drawCircle;

end


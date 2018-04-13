function h2s_radii(n,ds,type)

N = 100; % bootstraps

s = NaN(N,5,numel(ds)); % container for sampled estimators
colors = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1]};

for d = ds
   dindex = find(d==ds);
   % Sample & measure radii
   X  = sampleSpheres(n,d,N,type);
   radii = sqrt(sum(X.^2,2));
   % Estimators
   % expDistPerRad method
   [~,~,tmp] = arrayfun(@(i) estimateHypersphere(X(:,:,i)),1:N,'UniformOutput',false);
   s(:,1,dindex) = cell2mat(tmp);
   s(:,3,dindex) = median(radii,1);
   s(:,2,dindex) = s(:,3,dindex).*2^(1/d);   % median*2^(1/d)
   maxradii = max(radii,[],1);
   s(:,4,dindex) = maxradii + maxradii/n;                      % MVUE for uniform distribution
   s(:,5,dindex) = ones(1,N)*sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2)); % MVUE for gaussian distribution

   % Results
   figure(98);clf;plotEstimators(X,radii,mat2cell(s(:,:,dindex),100,ones(1,5)),colors);
%   skewness(radii(:))
end
1-[mn md mx]

return

end



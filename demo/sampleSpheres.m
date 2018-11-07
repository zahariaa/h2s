function X = sampleSpheres(n,d,N,type)
X = zeros(n,d,N);
switch lower(type)
   case 'uniform'
      for i = 1:N
         X(:,:,i) = randsphere(n,d,1);
      %   X(:,:,i) = randnball(n,d,1);
      end
   case 'gaussian'
      X = 1*randn(n,d,N);
   otherwise
      X = 1*rand( n,d,N)-0.5;
end

return
end

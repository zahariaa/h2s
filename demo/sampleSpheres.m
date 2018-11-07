function X = sampleSpheres(n,d,N,type)
X = zeros(n,d,N);
if strcmpi(type,'uniform')
   for i = 1:N
   %   X(:,:,i) = randsphere(n,d,1);
      X(:,:,i) = randnball(n,d,1);
   end
else
   for i = 1:N
      X(:,:,i) = 1*randn(n,d);
   end
end

return
end

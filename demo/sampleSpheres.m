function X = sampleSpheres(n,d,N,type,stream)

if ~exist('stream','var') || isempty(stream)
	stream = RandStream.getGlobalStream;
end

X = zeros(n,d,N);

switch lower(type)
   case 'uniform'
      for i = 1:N
         X(:,:,i) = randsphere(n,d,1,stream);
      %   X(:,:,i) = randnball(n,d,1,stream);
      end
   case 'gaussian'
      X = 1*randn(stream,n,d,N);
   otherwise
      X = 1*rand( stream,n,d,N)-0.5;
end

return
end

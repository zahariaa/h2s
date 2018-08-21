function v = nsimplex(n)
% vertices = nsimplex(number_of_dimensions)
% Generates (n+1) vertices for a unit n-simplex centered at the origin.
% The first vertex is always [1 zeroes...]
% E.g., a 3-simplex is a tetrahedron in 3D, so 4 vertices are generated
%
% 2018-06-11 AZ Created

% Initialize first row, and the rest as zeros
v = [ 1 -ones(1,n)/n
      zeros(n-1,n+1) ];
% Loop
for i = 2:n
   % all vertices are unit norm. get v_ii by solving pythagoeran theorem
   % using (only) already filled in values: assume other elements are zero
   v(i,i)       = sqrt(1 - sum(v(1:i-1,i).^2));
   % the rest of the row must sum to cancel out v_ii
   v(i,i+1:end) = -v(i,i)/(n-i+1);
end

return


%% DEBUG/DEMO
n=4;
hyp = SetOfHyps(nsimplex(n)',ones(n+1,1));
figure;hyp.h2s.show


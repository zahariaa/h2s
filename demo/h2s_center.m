function h2s_center

d = 3;
N = 100; % bootstraps

for n = 20
   X = sampleSpheres(n,d,N);
   mn = mean(mean(X,1),3);
   md = mean(median(X,1),3);
   
   figure;clf;hold on;
   plot(X(:,1,1),X(:,2,1),'ko');
   plot(mn(1),mn(2),'bx','MarkerSize',16);
   plot(md(1),md(2),'rx','MarkerSize',16);
   axis equal square
   drawCircle([0 0],1);
end
[norm(mn) norm(md)]

return

end


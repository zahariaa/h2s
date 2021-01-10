% QT_nBallVolume

clear;
dims = 1:1000;

for dimI = 1:numel(dims)
    dim = dims(dimI);
    vols(dim) = nBallVolume(dim,4);
end


figure(20); clf; set(20,'Color','w');
plot(dims,vols,'k.','MarkerSize',15);
xlabel('number of dimensions');
ylabel('volume');
title('\bfvolume of unit n-ball as a function of # dimensions');

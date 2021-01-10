% QT2_nBallVolume

clear;
dims = 1:20;
rads = [1 2 3];
rads = 0.1:0.1:0.5;
rads = [0.1:0.3:1, 2:3];
rads = [0.9 1.1];

cols = colourScale([0 0 0; .8 .8 .8],numel(rads));
figure(20); clf; set(20,'Color','w'); 

for radI = 1:numel(rads)
    rad = rads(radI)
    legendStrings{radI}=['radius=',num2str(rad)];
    for dimI = 1:numel(dims)
        dim = dims(dimI);
        vols(dimI) = nBallVolume(dim,rad);
    end
    semilogy(dims,vols,'.','Color',cols(radI,:),'MarkerSize',15); hold on;
    drawnow;
end

xlabel('number of dimensions');
ylabel('log(volume)');
title('\bfvolume of n-ball as a function of # dimensions');
legend(legendStrings,'Location','SouthWest');
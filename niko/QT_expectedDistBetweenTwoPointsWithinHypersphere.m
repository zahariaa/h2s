% QT_expectedDistBetweenTwoPointsWithinHypersphere

% conclusions
% the expected distance between two points drawn from a uniform
% hyperspherical distribution is proportional to the radius R of the
% hypersphere and grows (as a saturating function) with the dimensionality d.
%
% for d=2 (circular disk): E(dist) = (128*R)/(45*pi) = 0.9054*R
% for d=3 (3D sphere):     E(dist) = (36*R)/35 = 1.0286*R


%% control variables
Rs=(0:3)/3;
ds= [1 2 3 5 8 15 30 60 100 1000];
n=1000;

samplingStyle='polarAdjusted';
%samplingStyle='rejection';


%% loop over hypersphere radii Rs and dimensionalities ds
avgDists=nan(numel(ds),numel(Rs));
nSamplesAccepted=nan(numel(ds),numel(Rs));

for dI=1:numel(ds)
    for RI=1:numel(Rs)
        R=Rs(RI);
        d=ds(dI);


        % monte carlo
        switch samplingStyle
            case 'polarAdjusted';
                % polar sampling adjusted for uniformity
                directions = randn(n,d);
                radii = rand(n,1).^(1/d)*R;
                data=directions.*repmat(radii./sqrt(sum(directions.^2,2)),[1 d]);

            case 'rejection'
                % rejection sampling
                data = R*(rand(n,d)*2-1);
                r = sqrt(sum(data.^2,2));
                data=data(r<=R,:);
        end
        
        nSamplesAccepted(dI,RI) = size(data,1);
        avgDists(dI,RI) = mean(pdist(data,'Euclidean'));
        
        
    end
end


%% plot results
figure(100); clf; set(100,'color','w');
subplot(2,1,1);
cols=repmat(linspace(0,0.7,numel(ds))',[1 3]);

for dI=1:numel(ds)
    plot(Rs,avgDists(dI,:),'Color',cols(dI,:),'LineWidth',3); hold on;
end

xlabel('radius of the hypersphere');
ylabel('average distance between points');

legend({'d=1','d=2','d=3','...'},'Location','SouthEast');



%% plot expected distance between two points drawn uniformly from a hypersphere as a function of the dimension
n = 1000;
i=1;
expectedDists=[];

for d = ds
    
    points = randsphere(n,d,1);
    expectedDists(i) = mean(pdist(points,'Euclidean'));

    i=i+1;
end

subplot(2,1,2);
plot(ds,expectedDists,'k');
xlabel('dimensionality d');
ylabel({'expected distance between two points','uniformly drawn within a hypersphere'});




% TEST_hypersphere2sphere

% hs2s tests: comparison to MDS
%       - two touching equal-radius hyperspheres in 1-10 dimensions
%       - two concentric hyperspheres of different radius in 1-10 dimensions
%       - larger hypersphere enclosing smaller one, touching in one surface point
%       - two intersecting hypershperes
%       - 3 spheres (two intersecting) -> 2D circles
%       - 5 hyperspheres in 1-100 dimensions
%       - as previous, but illustrate effect of sample size changes



%% general control variables
scenarios = 1:6; %[1 2 3];
dimLow    = 2; % 2=showCircles, 3=showSpheres

psFilespec = '';
pdfFilespec = '';


%% scenario 1: two touching equal-radius hyperspheres in 1-200 dimensions
if find(scenarios==1)
    radius = 1;

    nsDim = [1 2 3 5 10 20 40 200]

    for nPointsPerCat = [10 40]
        h=figure(2001);clf;
        set(h,'Name','scenario 1: two touching equal-radius hyperspheres in 1-200 dimensions');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
            points = randsphere(2*nPointsPerCat,nDim,radius);
            %points = randn(2*nPointsPerCat,nDim);
            points(end/2+1:end,1) = points(end/2+1:end,1)+2; % shift second half by 2 along first dimension

            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end

%% scenario 2: two identical hyperspheres in 1-200 dimensions
if find(scenarios==2)
    radius1 = 1;
    radius2 = 1;

    nsDim = [1 2 3 5 10 20 40 200];

    for nPointsPerCat = [20 100]
        h=figure(2002);clf;
        set(h,'Name','scenario 2: two identical hyperspheres in 1-200 dimensions');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
            points = nan(2*nPointsPerCat,nDim);
            points(1:nPointsPerCat,:) = randsphere(nPointsPerCat,nDim,radius1);
            points(nPointsPerCat+1:end,:) = randsphere(nPointsPerCat,nDim,radius2);
            
            % gaussian instead: points = randn(2*nPointsPerCat,nDim);

            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end


%% scenario 3: two concentric hyperspheres of different radius in 1-200 dimensions
if find(scenarios==3)
    radius1 = 1;
    radius2 = 2;

    nsDim = [1 2 3 5 10 20 40 200];

    for nPointsPerCat = [50]
        h=figure(2003);clf;
        set(h,'Name','scenario 3: two concentric hyperspheres of different radius in 1-200 dimensions');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
            points = nan(2*nPointsPerCat,nDim);
            points(1:nPointsPerCat,:) = randsphere(nPointsPerCat,nDim,radius1);
            points(nPointsPerCat+1:end,:) = randsphere(nPointsPerCat,nDim,radius2);
            %points = randn(2*nPointsPerCat,nDim);

            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end

%% scenario 4: larger hypersphere enclosing smaller one, touching in one surface point
if find(scenarios==4)
    radius1 = 1;
    radius2 = 2;

    nsDim = [1 2 3 5 10 20 40 200];

    for nPointsPerCat = [100]
        h=figure(2004);clf;
        set(h,'Name','scenario 4: larger hypersphere enclosing smaller one, touching in one surface point');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
            points = nan(2*nPointsPerCat,nDim);
            points(1:nPointsPerCat,:) = randsphere(nPointsPerCat,nDim,radius1);
            points(nPointsPerCat+1:end,:) = randsphere(nPointsPerCat,nDim,radius2);
            points(nPointsPerCat+1:end,1) = points(nPointsPerCat+1:end,1) + 1;
            %points = randn(2*nPointsPerCat,nDim);

            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end



%% scenario 5: two intersecting hypershperes
if find(scenarios==5)
    radius1 = 2;
    radius2 = 2;

    nsDim = [1 2 3 5 10 20 40 200];

    for nPointsPerCat = [50]
       h=figure(2005);clf;
        set(h,'Name','scenario 5: two intersecting hypershperes');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
            points = nan(2*nPointsPerCat,nDim);
            points(1:nPointsPerCat,:) = randsphere(nPointsPerCat,nDim,radius1);
            points(nPointsPerCat+1:end,:) = randsphere(nPointsPerCat,nDim,radius2);
            points(nPointsPerCat+1:end,1) = points(nPointsPerCat+1:end,1) + 1;
            %points = randn(2*nPointsPerCat,nDim);

            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end



%% scenario 6: as previous, but illustrate effect of sample size changes
if find(scenarios==6)
    radius1 = 2;
    radius2 = 2;
    
    fac = 5;

    nsDim = [1 2 3 5 10 20 40 200];

    for nPointsPerCat = [20]
       h=figure(2006);clf;
        set(h,'Name','scenario 6: two intersecting hyperspheres with different numbers of points');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
            points = nan((1+fac)*nPointsPerCat,nDim);
            points(1:nPointsPerCat,:) = randsphere(nPointsPerCat,nDim,radius1);
            points(nPointsPerCat+1:end,:) = randsphere(fac*nPointsPerCat,nDim,radius2);
            points(nPointsPerCat+1:end,1) = points(nPointsPerCat+1:end,1) + 1;
            %points = randn(2*nPointsPerCat,nDim);

            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));
            categories.vectors = [categories.vectors; repmat(categories.vectors(end/2+1:end,:),[(fac-1) 1])];
            
            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end


%% scenario 7: gaussian ellipsoid version of scenario 2
%%           : two overlapping, random covariance hyperellipsoids in 1-200 dimensions
if find(scenarios==1)
    radius = 1;

    nsDim = [1 2 3 5 10 20 40 200]

    for nPointsPerCat = [10 40]
        h=figure(2007);clf;
        set(h,'Name','scenario 7: 2 overlap hyperellipsoids w rand covar 1-200 dims');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
	    SIG{1} = rand(nDim);  SIG{1} = SIG{1}*SIG{1}';
	    SIG{2} = rand(nDim);  SIG{2} = SIG{2}*SIG{2}';
            points = nan(nPointsPerCat,nDim);
	    points(1:nPointsPerCat,:) = mvnrnd(zeros(nDim,1),SIG{1},nPointsPerCat);
	    points(nPointsPerCat+1:2*nPointsPerCat,:) = ...
	                               mvnrnd(zeros(nDim,1),SIG{2},nPointsPerCat);
            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end


%% scenario 8: gaussian ellipsoid version of scenario 5
%%           : two intersecting equal-covariance hyperellipsoids in 1-200 dimensions
if find(scenarios==1)
    radius = 1;

    nsDim = [1 2 3 5 10 20 40 200]

    for nPointsPerCat = [10 40]
        h=figure(2008);clf;
        set(h,'Name','scenario 8: 2 adj hyperellipsoids w same rand covar 1-200 dims');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
	    SIG{1} = rand(nDim);  SIG{1} = SIG{1}*SIG{1}';
            points = nan(nPointsPerCat,nDim);
	    points(1:nPointsPerCat,:) = mvnrnd( ones(nDim,1),SIG{1},nPointsPerCat);
	    points(nPointsPerCat+1:2*nPointsPerCat,:) = ...
	                                mvnrnd(-ones(nDim,1),SIG{1},nPointsPerCat);
            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end


%% scenario 9: gaussian ellipsoid version of scenario 6
%%           : two overlapping, identical hyperellipsoids in 1-200 dimensions
%%             different sampling
if find(scenarios==1)
    radius = 1;
    fac    = 5;
    nsDim = [1 2 3 5 10 20 40 200]

    for nPointsPerCat = [20]
        h=figure(2009);clf;
        set(h,'Name','scenario 9: 2 overlap hyperellipsoids w rand covar 1-200 dims');
        for nDimI = 1: numel(nsDim)
            nDim = nsDim(nDimI);
	    SIGMA{1} = rand(nDim);  SIGMA{1} = SIGMA{1}*SIGMA{1}';
            points = nan((1+fac)*nPointsPerCat,nDim);
	    points(1:nPointsPerCat,:) = mvnrnd( ones(nDim,1),SIG,nPointsPerCat);
	    points(nPointsPerCat+1:(1+fac)*nPointsPerCat,:) = ...
	                                mvnrnd(-ones(nDim,1),SIG,fac*nPointsPerCat);
            categories.labels = {'category 1','category 2'};
            categories.colors = [0.8 0 0; 0 0 0];
            categories.vectors = logical(blockDiagonalMatrix(2*nPointsPerCat,2,2));

            figPanelSpec = [h 4 6 1+(nDimI-1)*3];
            titleStr = any2str(nPointsPerCat, ' points/cat. in ',nDim,' dim.');
            HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
        end
    end
end

%% scenario 7: 3 spheres (two intersecting) -> 2D circles
%% scenario 8: 5 hyperspheres in 1-100 dimensions



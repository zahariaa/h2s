function mat = blockDiagonalMatrix(height,width,nBlocks)

blockMat = eye(nBlocks);
mat = interp2(1:nBlocks,1:nBlocks,blockMat,linspace(1,nBlocks,width),linspace(1,nBlocks,height)','nearest');

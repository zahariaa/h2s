function vol = nBallVolume(nDim,radius)

vol = pi^(nDim/2) / gamma(nDim/2+1) * radius^nDim;

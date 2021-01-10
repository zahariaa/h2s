function V = hypersphereVolume(radius,dimensions)
% Minifunction calculating hypersphere volume in arbitrary dimensions
% 2018-06-05 AZ Created

V = (radius.^dimensions).*(pi.^(dimensions/2))./gamma((dimensions/2)+1);

return


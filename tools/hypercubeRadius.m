function r = hypercubeRadius(d)
% hypercubeRadius: calculates expected radius of a uniform hypercube


% Hypercube correction factor from:
% https://math.stackexchange.com/questions/2446084/distance-from-the-centre-of-a-n-cube-as-n-rightarrow-infty
% integral2(@(x,y)   sqrt(x.^2 + y.^2)       ,-0.5,0.5,-0.5,0.5)          % d=2
% integral3(@(x,y,z) sqrt(x.^2 + y.^2 + z.^2),-0.5,0.5,-0.5,0.5,-0.5,0.5) % d=3
% integral(@(w) integral3(@(x,y,z) sqrt(w.^2 + x.^2 + y.^2 + z.^2),-0.5,0.5,-0.5,0.5,-0.5,0.5),...
%          -0.5,0.5,'ArrayValued',true) % d=4
% integral2(@(v,w) arrayfun(@(v,w) integral3(@(x,y,z) sqrt(v.^2 + w.^2 + x.^2 + y.^2 + z.^2),...
%          -0.5,0.5,-0.5,0.5,-0.5,0.5),v,w),-0.5,0.5,-0.5,0.5) % d=5
% Don't run this one, it takes forever:
% integral2(@(s,t) arrayfun(@(s,t) integral3(@(u,v,w) arrayfun(@(u,v,w) ...
%           integral3(@(x,y,z) sqrt(s.^2 + t.^2 + u.^2 + v.^2 + w.^2 + x.^2 + y.^2 + z.^2),...
%          -0.5,0.5,-0.5,0.5,-0.5,0.5),u,v,w),-0.5,0.5,-0.5,0.5,-0.5,0.5),s,t),-0.5,0.5,-0.5,0.5) % d=8
expectedDists = [0.25 0.3826 0.4803 0.5609]; 

if d<5,   r = expectedDists(d);
else      r = sqrt(d/12);
end
return
end


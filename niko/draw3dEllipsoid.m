function h=draw3dEllipsoid(location,covariance,col,n,opacity)


%% defaults
if ~exist('location'  ,'var') || isempty(location  ), location = [0 0 0];  end
if ~exist('col'       ,'var') || isempty(col       ), col = [0.5 0.5 0.5]; end
if ~exist('covariance','var') || isempty(covariance), covariance = eye(3); end
if ~exist('n'         ,'var') || isempty(n         ), n = 50;              end
if ~exist('opacity'   ,'var') || isempty(opacity   ), opacity = 0.5;       end


%% find principal axes and corresponding radii
[eigenvec eigenval]=pcacov(covariance);
eigen = eigenvec.*repmat(sqrt(eigenval)',[3 1]);


%% contruct sphere and distort its coordinates
[X,Y,Z] = sphere(n);
xyz2 = eigen(:,1)*X(:)' + eigen(:,2)*Y(:)' + eigen(:,3)*Z(:)';
X2 = location(1) + reshape(xyz2(1,:),size(X));
Y2 = location(2) + reshape(xyz2(2,:),size(Y));
Z2 = location(3) + reshape(xyz2(3,:),size(Z));


%% draw 3d ellipsoid
% figure(1); clf;

try
    h=surf(X2,Y2,Z2,'FaceColor',col,'EdgeColor','none','FaceAlpha',opacity);
catch
    disp('Could not draw.')
end

    % h=surf(X2,Y2,Z2,'FaceColor',col,'EdgeColor','none');
% alpha(h,opacity);

% CHECK
% for i = 1:350;
%     x=[X2(i) Y2(i) Z2(i)];
%     c=location;
%     val(i)=(x-c)*inv(covariance)*(x-c)'
% end
% val

% lighting phong;
% camlight('left');
% set(gcf,'Renderer','OpenGL');
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal;

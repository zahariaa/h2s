function h = draw3dEllipsoid(location,covariance,col,n,opacity,RECURSED)

%% defaults
if ~exist('location'  ,'var') || isempty(location  ), location = [0 0 0];  end
if ~exist('col'       ,'var') || isempty(col       ), col = [0.5 0.5 0.5]; end
if ~exist('covariance','var') || isempty(covariance), covariance = eye(3); end
if ~exist('n'         ,'var') || isempty(n         ), n = 50;              end
if ~exist('opacity'   ,'var') || isempty(opacity   ), opacity = 0.5;       end

%% Recurse if given multiple inputs
nl = numel(location  )/3;
nc = numel(col       )/3;
if iscell(covariance), nr = numel(covariance);
else                   nr = numel(covariance)/9;
end
mn = max([nl nc nr]);
if mn > 1
   if nl>1,     location   = num2cell(location,2);
   else         location   = repmat({ location(:)'},[mn 1]);
   end
   if iscell(covariance) && nr==mn, covariance = covariance(:);
   elseif nr>1, covariance = num2cell(covariance,3);
   else         covariance = repmat({ covariance  },[mn 1]);
   end
   if nc>1,     col        = num2cell(col,2);
   else         col        = repmat({ col(:)'     },[mn 1]);
   end
   if numel(n      ) < mn;   n       = repmat({n      },[mn 1]);   end
   if numel(opacity) < mn;   opacity = repmat({opacity},[mn 1]);   end
   % Recursive call, with RECURSED flag set to true ...
   h = cellfun( @draw3dEllipsoid, location,covariance,col,n,opacity,true);
   camlight('left'); lighting phong; % ... so that lighting is only applied 1x
   return
elseif ~exist('RECURSED','var'), RECURSED = false;
end

%% find principal axes and corresponding radii
[V,D] = eig(covariance);
eigen = V.*repmat(sqrt(diag(D)),[1 3]);

%% contruct sphere and distort its coordinates
[X,Y,Z] = sphere(n);
xyz2 = eigen(:,1)*X(:)' + eigen(:,2)*Y(:)' + eigen(:,3)*Z(:)';
X2 = location(1) + reshape(xyz2(1,:),size(X));
Y2 = location(2) + reshape(xyz2(2,:),size(Y));
Z2 = location(3) + reshape(xyz2(3,:),size(Z));

%% draw 3d ellipsoid
hold on;
h = surf(X2,Y2,Z2,'FaceColor',col,'EdgeColor','none','FaceAlpha',opacity);
axis vis3d off

% Make sure lighting is only created after last call
if ~RECURSED,   lighting phong;   camlight('left');   end


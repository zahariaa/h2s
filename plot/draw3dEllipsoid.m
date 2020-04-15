function [h,XYZ] = draw3dEllipsoid(center,covariance,color,alpha,n,SETLIGHT)
% [h,XYZ] = draw3dEllipsoid(center=[0 0 0],covariance=eye(3),...
%                           color=[0.5 0.5 0.5],alpha=0.2,n=100,SETLIGHT=true)
% -  center must be Nx3
% -  accepts multiple inputs and draws multiple ellipsoids
% 
% e.g.:
% [h,XY] = draw3dEllipsoid;           % draws a black ellipsoid at [0 0 0] with
%                                     % covariance = eye(3), outputs figure
%                                     % handle and XYZ values of ellipsoid
% draw3dEllipsoid([1 2 4])
% draw3dEllipsoid([1 2 4],diag(1:3))            % ellipsoid with 1, 2, and 3
%                                               % for major/minor axes
% draw3dEllipsoid([1 2 4],3*eye(3))             % sphere with radius 3
% draw3dEllipsoid([1 2 4],3*eye(3),'r')         % red sphere
% draw3dEllipsoid([1 2 4],3*eye(3),[1 0 0])     % same as previous
% draw3dEllipsoid([1 2 4],3*eye(3),'r',0)       % opaque sphere
% draw3dEllipsoid([1 2 4],3*eye(3),'r',[],500)  % hi-res!
% draw3dEllipsoid([1 2 4;3 4 7;5 6 9],{3*eye(3);4*eye(3);5*eye(3)},'krg')
                                                % 3 ellipsoids, in black, red,
                                                % and green
% draw3dEllipsoid([1 2 4;3 4 7;5 6 9])          % 3 spheres, all black and radius
%                                               % 1, but with specified centers
% 
% All inputs are optional:
%    center (DEFAULT = [0 0]): [n x 3] matrix of center coordinate(s) for n
%       number of ellipsoid(s). n should match n in radius/color argument(s), or
%       equal 1, in which case, that center coordinate is used for all ellipsoids.
%    covariance (DEFAULT = eye(3)): 3x3 matrix (or n-cell of 3x3 matrices)
%       specifying the covariance matrices for n number of ellipsoids. n should
%       match  n should match n in radius argument, or equal 1,
%       in which case, that center coordinate is used for all ellipsoids.
%    color (DEFAULT = [0.5 0.5 0.5]): [n x 3] matrix or n-vector of strings for
%       the color(s) of the ellipsoid(s). n should match n in center/radius
%       argument(s), or equal 1, in which case, that center coordinate is used
%       for all ellipsoids.
%    alpha (DEFAULT = 0.2): The number, in [0,1], determining the ellipsoid(s)'
%       FaceAlpha/opacity. Same value applied to all ellipsoids.
%    n (DEFAULT = 100): Level of surface detail, as measured by the number of
%       vertices for the ellipsoids.
%    SETLIGHT (DEFAULT = true): sets lighting for all ellipsoids in the plot.
%       This is mainly a placeholder; it's forced to be false in recursive calls
%       to avoid multiple lighting sources piling up.
% 
% SEE ALSO DRAWCIRCLE, SHOWMODEL

%% defaults
if ~exist('center'    ,'var') || isempty(center    ), center = [0 0 0];      end
if ~exist('covariance','var') || isempty(covariance), covariance = eye(3);   end
if ~exist('color'     ,'var') || isempty(color     ), color = [0.5 0.5 0.5]; end
if ~exist('n'         ,'var') || isempty(n         ), n = 100;               end
if ~exist('alpha'     ,'var') || isempty(alpha     ), alpha = 0.5;           end

%% Recurse if given multiple inputs
nl = numel(center)/3;
nc = numel(color )/3;
if iscell(covariance),  nr = numel(covariance);
else                    nr = numel(covariance)/9;
end
mn = max([nl nc nr]);
if mn > 1
   if nl>1,     center     = num2cell(center,2);
   else         center     = repmat({center(:)'},[mn 1]);
   end
   if iscell(covariance) && nr==mn, covariance = covariance(:);
   elseif nr>1, covariance = vectify(num2cell(covariance,[1 2]));
   else         covariance = repmat({covariance},[mn 1]);
   end
   if nc>1,     color      = num2cell(color,2);
   else         color      = repmat({ color(:)'},[mn 1]);
   end
   if numel(n    ) < mn;   n     = repmat({n    },[mn 1]);   end
   if numel(alpha) < mn;   alpha = repmat({alpha},[mn 1]);   end
   % Recursive call, with SETLIGHT flag set to false ...
   [h,XYZ] = cellfun(@draw3dEllipsoid,...
                     center,covariance,color,alpha,n,repmat({false},[mn 1]),'UniformOutput',false);
   lighting phong; camlight('headlight'); view(180,90); % ... so that lighting is only applied 1x
   h  = vertcat (h{:});   XYZ = cat(1,XYZ{:});
   return
elseif ~exist('SETLIGHT','var'), SETLIGHT = true;
end

%% find principal axes and corresponding radii
[V,D] = eig(covariance);
eigen = V.*repmat(sqrt(diag(D)),[1 3]);

%% contruct sphere and distort its coordinates
[X,Y,Z] = sphere(n);
xyz2 = eigen(:,1)*X(:)' + eigen(:,2)*Y(:)' + eigen(:,3)*Z(:)';
X2 = center(1) + reshape(xyz2(1,:),size(X));
Y2 = center(2) + reshape(xyz2(2,:),size(Y));
Z2 = center(3) + reshape(xyz2(3,:),size(Z));

XYZ = cat(3,X2,Y2,Z2);
%% draw 3d ellipsoid
hold on;
h = surf(X2,Y2,Z2,'FaceColor',color,'EdgeColor','none','FaceAlpha',alpha);
axis vis3d off

% Make sure lighting is only created after last call
if SETLIGHT,   lighting phong;   camlight('headlight');   end


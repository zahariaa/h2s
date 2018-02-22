function h = drawCircle(center,radius,color,border,alpha,n)
% h = drawCircle(center,radius,<color='k'>,<border=2>,<alpha=0.2>,<n=100>)
% center must be Nx2
% 
% 2018-02-22 AZ Created

%% Initialization
if ~exist('color' ,'var') || isempty(color );   color  = [0 0 0];   end
if ~exist('n'     ,'var') || isempty(n     );   n      = 100;       end
if ~exist('alpha' ,'var') || isempty(alpha );   alpha  = 0.2;       end
if ~exist('border','var') || isempty(border);   border = 2;         end

%% Recurse
nr = numel(radius);
nc = numel(center)/2;
switch (nr>1)*10+(nc>1)
   case  1;   h = arrayfun(@(c1,c2  ) drawCircle([c1 c2],radius,color,border,alpha,n),...
                  center(:,1),center(:,2)); return;
   case 10;   h = arrayfun(@(      r) drawCircle(center,r,color,border,alpha,n),...
                  radius); return;
   case 11;   h = arrayfun(@(c1,c2,r) drawCircle([c1 c2],r,color,border,alpha,n),...
                  center(:,1),center(:,2),radius(:)); return;
end

%% Generate plot
if border==0;   opts = {'LineStyle','none','FaceAlpha',alpha};
else            opts = {'LineWidth',border,'FaceAlpha',alpha};
end

t = linspace(0,2*pi,n);
h = fill(radius*cos(t)+center(1),radius*sin(t)+center(2),color);
set(h,opts{:});
axis square off
hold on

return

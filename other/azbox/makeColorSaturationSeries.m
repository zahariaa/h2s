function colors = makeColorSaturationSeries(RGBin,n)
% Makes color series for plotting lines:
% e.g., for blue input, colors are logspaced from
% almost white to blue to almost black
% 
% 2013-02-13 AZ Created
% 2013-03-25 AZ Now uses CIELAB if image processing toolbox is present

% ColorOn  =           RGorB;
% ColorOff = [1 1 1] - RGorB;
% 
% nBright = ceil(2*n/3);
% nDark   = n - nBright - 1;
% 
% % cold = [0.75.^(2:n) 0];
% colsBright = logspace(-0.1,-1,nBright);
% colsDark   = logspace(-0.2,log10(0.4),nDark);
% 
% % Build color matrix
% colors = repmat(ColorOn,[nBright 1])+colsBright'*ColorOff;
% colors = [colors; ColorOn; ...
%          repmat([0 0 0],[nDark 1])+colsDark'*ColorOn ];
% 
% % Get rid of extra color if n is even
% colors = colors(1:n,:);


% % Convert RGBin to HSV
% H = rgb2hsv(RGBin);
% maxS = H(2); % get max distance from white
% H = repmat(H(1),[n 1]);
% 
% % Build parabola through SV slice of HSV cylinder
% buffer = 1/3;
% V = linspace(buffer,1-buffer,n)';
% 
% % Parabola steepness coefficient
% a = maxS/(0.5)^2;
% S = (maxS - a*(V-0.5).^2);
% 
% colors = hsv2rgb([H S V]);


% L1 path white->color->black

if n == 1
   % do nothing
   colors = RGBin;

else
   % Detect black or white input, to cut out opposite later on
   if ~any(RGBin) || all(RGBin==1);   BW = true;  n = n+1;
   else                               BW = false;
   end

   try
      % if have image processing toolbox
      rl = makecform('srgb2lab');
      lr = makecform('lab2srgb');

      LABin  = applycform(RGBin,rl);

      colors = applycform([linspace(0,100,n);LABin(2)*ones(1,n);LABin(3)*ones(1,n)]',lr);

   %    colors(:,RGBin~=0) = repmat(mean(colors(:,RGBin~=0),2),[1 sum(RGBin~=0)]);
   catch

      himax = 0.2;
      lomax = 0.4;

      H = rgb2hsv(RGBin);
      maxS = H(2); % get max distance from white
      H = repmat(H(1),[n 1]);


      nBright = ceil(2*n/3);
      nDark   = n - nBright - 1;

      colsBright = logspace(-1,log10(1-himax),nBright);
      colsDark   = logspace(-0.2,log10(lomax),nDark);

      S = [colsBright'; maxS*ones(nDark+1,1)];
      try
      V = [ones(nBright+1,1); colsDark'];
      catch
      V = [ones(nBright  ,1); colsDark'];
      end

      % 0.25->1->0.25

      colors = hsv2rgb(flipud([H S V]));
   end
   
   % Cut out pure white if input is black, and vice versa
   if BW
      toCut  = maxix(abs(colors - repmat(RGBin,[n 1])));
      colors = colors(setdiff(1:n,toCut),:);
   end
   
end
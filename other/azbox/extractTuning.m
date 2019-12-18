function varargout = extractTuning(p,subset,SEM,SPIKERATE)
% Extract tuning parameters from p struct
% 
% Ex. 1:
% [bins,respmean,respstd,spontmean,spontstd] = extractTuning(p);
% 
% Ex. 2: Use standard error instead of standard deviation
% [bins,respmean,respstd,spontmean,spontstd] = extractTuning(p,[],true);
% 
% Ex. 3: Output to a single struct
% tuning = extractTuning(p);
%
% 2011-08-05 AZ Created
% 2011-11-11 AZ added option to select only a subset of the data
% 2011-11-22 AZ added handling rfsize9c (SIZE TUNING ONLY, NO AMRF (?) )
% 2012-01-30 AZ added ori16 (same as sf & tf)
%               made plaid120 give 2 rows of tuning [grating;plaid]
% 2012-12-02 AZ bugfix
% 2012-12-10 AZ added input arg to use standard error instead of standard
%               deviation

if ~exist('subset','var') || isempty(subset);	subset = true(size(p.nspikes'));   end
if ~exist('SEM'   ,'var') || isempty(SEM   );      SEM = false;                    end
if  exist('SPIKERATE','var') && ~isempty(SPIKERATE) && SPIKERATE
       durcorr   = p.duration(p.stim);
       durcorr   = mean(durcorr(subset)/10000);
       p.nspikes = p.nspikes(subset)/durcorr;
else   p.nspikes = p.nspikes(subset);durcorr = 1;
end
% Select subset
for feld = {'con','ori','sf','tf','ph','size','special'};feld=feld{1}; %#ok<FXSET>
   if numel(p.(feld)) > numel(p.run); p.(feld) = p.(feld)(subset,:); end
end
% NOTE: p fields unchanged: startTime, endTime, duration, plaid
%% Identify errors
matrixBlockIDs = unique(p.blockID(p.stim))';
blocks         = p.blockID(p.stim);
blocks         = blocks(subset);
cons           = sum(p.con,2)'; % Sum contrasts for plaids

spikeRate = p.nspikes(cons==0);
spontmean = mean(spikeRate);
spontstd  =  std(spikeRate);
if SEM
   spontstd = spontstd/sqrt(sum(cons==0));
end
stimmean  = nanmean(p.nspikes(abs(cons)>0)); % don't divide by mduration

%% Protocol-specific settings
fx = 1;
nSubsets = 1;
feld = @(i) 'special';

switch p.protocol
	case {'sf11','sf11c','sf11_B','sf11c_B','tf10','tf10low','tf10high','tf10_B','tf10low_B','tf10high_B'}
      feld = @(i) p.protocol(1:2);
   case {'rvc10','rvc10_B','rvcSur','rvcOrtho'}
      feld = @(i) 'con';
   case {'rfsize9c','rfsize9c_B'}
      nSubsets = 2;
   case {'plaid120','plaid120SS','plaid120_B','plaid120SS_B','plaid120_16'}
      fx   = 2;
      nSubsets = numel(unique(p.specialDescription(regexp(p.specialDescription,'[0-9]\='))));
      nBins    = numel(matrixBlockIDs)/nSubsets;
   case {'ori16','ori16_B','ori16c','ori16c_B'}
      feld = @(i) 'ori';
   case {'ValSpace'}
      feld     = @(i) valspacefield(p,blocks,i);
      nSubsets = sum(sum(p.special,1)>0);
      nBins    = round(numel(matrixBlockIDs)/nSubsets);
      p.nan    = NaN(size(p.con));
      p.specialDescription = p.specialDescription';
end

%% Extract Responses
bins     = NaN(size(matrixBlockIDs));
respmean = NaN(size(matrixBlockIDs));
respstd  = NaN(size(matrixBlockIDs));
series   = NaN(size(matrixBlockIDs));%debug
oracle   = NaN(size(p.nspikes     ));
binBlks  = NaN(size(matrixBlockIDs));
null     = oracle;
oracle(cons==0) = spontmean;
  null(cons==0) = spontmean;
for i = matrixBlockIDs(:)'
   if ~any(isnan(p.(feld(i))(blocks==i,fx))) && ismember(i,unique(blocks))
          bins(i-min(matrixBlockIDs)+1) = unique(p.(feld(i))(blocks==i,fx));
      binBlks( i-min(matrixBlockIDs)+1) = i;
      spikeRate = p.nspikes(blocks==i & cons ~= 0);
      respmean(i-min(matrixBlockIDs)+1) =   mean(spikeRate);
      respstd( i-min(matrixBlockIDs)+1) =    std(spikeRate);
      if SEM
       respstd(i-min(matrixBlockIDs)+1) = respstd(i-min(matrixBlockIDs)+1) / ...
                                              sqrt(sum(blocks==i & cons ~= 0));
      end
      series(  i-min(matrixBlockIDs)+1) = unique(p.(feld(i))(blocks==i,1));%debug
      if ~isempty(spikeRate)
         oracle(blocks==i) = respmean(i-min(matrixBlockIDs)+1);
           null(blocks==i) = stimmean; % want this to be raw spike rate
      end
   end
end

%% Reshape
switch p.protocol
   case {'sf11','sf11c','sf11_B','sf11c_B',...
         'tf10','tf10low','tf10high','tf10_B','tf10low_B','tf10high_B',...
         'rvc10','rvc10_B','rvcSur','rvcOrtho','ori16','ori16_B','ori16c','ori16c_B'}
      respmean = vectify(respmean(bins~=0));
      respstd  = vectify( respstd(bins~=0));
%       series   = vectify(  series(bins~=0));%debug
      binBlks  = vectify( binBlks(bins~=0));
      bins     = vectify(    bins(bins~=0));
   case {'plaid120','plaid120SS','plaid120_B','plaid120SS_B','plaid120_16'}
      bins     = reshape(bins          ,nBins,nSubsets);
      respmean = reshape(respmean      ,nBins,nSubsets);
      respstd  = reshape(respstd       ,nBins,nSubsets);
      series   = reshape(series        ,nBins,nSubsets);%debug
      binBlks  = reshape(binBlks       ,nBins,nSubsets);
      %% cut blanks and rearrange
      nBins    = max(sum(~isnan(series)));
      bins     =     bins(1:nBins,fliplr(find(sum(~isnan(series))>0)));
      respmean = respmean(1:nBins,fliplr(find(sum(~isnan(series))>0)));
      respstd  = respstd( 1:nBins,fliplr(find(sum(~isnan(series))>0)));
   case {'rfsize9c','rfsize9c_B'}
      n        = numel(bins);
      bins     = reshape(bins(          1:n-1),[],nSubsets);
      respmean = reshape(respmean(      1:n-1),[],nSubsets);
      respstd  = reshape(respstd(       1:n-1),[],nSubsets);
%       series   = reshape(series(        1:n-1),[],nSubsets);%debug
      %% cut blanks
      bins(bins<0) = NaN; % dropping both blank blocks!
      binBlks  = reshape(       binBlks(1:n-1),[],nSubsets);
   case {'ValSpace'}
      nSubsets = numel(p.specialDescription);
      % arrange blocks according to experiment series
      blks     = arrayfun(@(i) unique(p.blockID([false(~p.stim(1));p.special(:,i)]))',1:nSubsets,'UniformOutput',false);
      blks     = cell2mat(cellfun(@(x) [x;NaN(13-numel(x),1)],blks,'UniformOutput',false));
      blksPres = unique(p.blockID([false(~p.stim(1));~p.con]))';
      if any(blksPres)
         blks(blks==0) = indexm(blksPres,1:sum(blks(:)==0));
         % reshape
         bins     = bins(    blks);
         respmean = respmean(blks);
         respstd  = respstd( blks);
         binBlks  = binBlks( blks);
      end
      % package specialDescription cell into a nicely readable string
      str=[];for i = 1:nSubsets;str = sprintf('%s%2.0f = %s\n',str,i,p.specialDescription{i});end
      p.specialDescription = str;
   otherwise % (including HyperPlaid*)
      oracle   = p.nspikes;
end

% ORACLE AND NULL MODELS SHOULD BE IN COUNTS, NOT RATE
oracle = oracle*durcorr;
null   =   null*durcorr;

% Combine outputs into struct if desired
if nargout == 1
   varargout{1} = struct('bins',bins,'respmean',respmean,'respstd',respstd,...
                 'spontmean',spontmean,'spontstd',spontstd,'blocks',blocks,...
                 'binBlks',binBlks,'oracle',oracle,'null',null,'protocol',p.protocol,...
                 'specialDescription',p.specialDescription);
else
   varargout    = {bins;respmean;respstd;spontmean;spontstd;blocks;oracle;null;binBlks};
end

end

%% special function to determine which 
function feld = valspacefield(p,blocks,i)
   felds = [{'ori','ori','sf','sf','tf','ori','ori','ori','ori'},...
      cellfun(@(x) lower(x(1:2)),p.specialDescription(10:end),'UniformOutput',false)'];
   
   feldnum = find(logical(sum(p.special(blocks==i,:))));
   if any(feldnum);   feld = felds{feldnum};
   else               feld = 'nan';
   end
end

function data = mattloader(subject,session)
% e.g., data = mattloader('Vortex','20180808');

%datadir = ['~' filesep 'Desktop' filesep];
datadir = '/Volumes/kriegeskorte-labshare/OxfordMonkeyData/RSVP_700_Exp2/';

%% Load stimuli
stimfile = dir(sprintf('%sRSVP_700-%s-%s-%s-%s*.mat',datadir,subject,...
             session(5:6),session(7:8),session(1:4)));
load([stimfile.folder filesep stimfile.name],'a');

% Image file key
imkey = cellfun(@(x) str2double(x(5:8)),a.TaskObject(:,2));
stims = imkey(a.ConditionNumber)';

%% Load spikes
spikefile = dir(sprintf('%s%s%s%s_%s_*trialised.mat',...
         datadir,session(7:8),session(5:6),session(1:4),subject));
load([spikefile.folder filesep spikefile.name],'SpikeData')

%% Stimulus start and end times
stimstarts = arrayfun(@(x) x.behavMatrix(...
                     find(x.behavMatrix(:,2)==1,1,'last'),1),...
                     SpikeData,'UniformOutput',false);
stimends   = arrayfun(@(x) x.behavMatrix(...
                     find(x.behavMatrix(:,2)==2,1,'last'),1),...
                     SpikeData,'UniformOutput',false);
stimstarts(cellfun(@isempty,stimstarts)) = {NaN};
stimends(  cellfun(@isempty,stimends  )) = {NaN};
stimstarts = [stimstarts{:}];
stimends   = [stimends{:}];

%% Valid trials
valid = arrayfun(@(x) ~any(  x.behavMatrix(:,2)==15 ...
                           | x.behavMatrix(:,2)==17) ...
                      && any(x.behavMatrix(:,2)==34), SpikeData);
valid(any([stimstarts;stimends]==0)) = false;

nt    = numel(SpikeData);
%% Spike times
% Take cells of NaN-padded matrices and convert them to
% NaN-stripped cells, while also subtracting the stimulus start time
spiketimes.v4 = arrayfun(@(i) ...
                 cellfun(@(x) x(~isnan(x))-double(stimstarts(i)),...
    num2cell(SpikeData(i).V4spikeData,2),'UniformOutput',false),...
    1:nt,'UniformOutput',false);
spiketimes.te = arrayfun(@(i) ...
                 cellfun(@(x) x(~isnan(x))-double(stimstarts(i)),...
    num2cell(SpikeData(i).TEpspikeData,2),'UniformOutput',false),...
    1:nt,'UniformOutput',false);
spiketimes  = structfun(@cellify,spiketimes,'UniformOutput',false);
spikecounts = structfun(@(x) cellfun(@numel,x),spiketimes,'UniformOutput',false);

%% FORMAT OUTPUT
data = struct('stims',stims,'valid',valid,'stimstarts',stimstarts,...
              'stimends',stimends,'spikecounts',spikecounts,...
              'protocol','trans92');
return


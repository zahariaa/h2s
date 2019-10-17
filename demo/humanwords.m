% humanwords.m : Loads Laura's data and runs h2s on it
% for different cuts/categories of words

FAST              = true;
LOBES_NOT_REGIONS = true;

%% LOAD AND RESHAPE DATA
fprintf('Loading...');
load('~/Desktop/humanwords/datan.mat')
load('~/Desktop/humanwords/elec.mat')
load('~/Desktop/humanwords/features.mat')
fprintf(' done.\n');

if LOBES_NOT_REGIONS
   [regions,~,ir] = unique({elec.Loc2});
else
   [regions,~,ir] = unique({elec.Loc3});
end
nr = numel(regions);

datan = permute(datan,[3 1 2]);

%% MANUALLY MODIFY LAURA'S SEMANTIC CATEGORIES
% nountype mod
nountypeAZ      = nountype;
nountypeAZnames = {'count' 'mass' 'countAndCollective' 'CountAndMass'};
nountypeAZ(:,4)=false;               % drop plural (only 3 words)
ix = sum(nountypeAZ(:, 1:2 ),2)==2;  % find count and mass words
nountypeAZ(ix,4)   = true;           % replace plural with countAndMass
nountypeAZ(ix,1:2) = false;          % drop countAndMass from inividual count and mass
ix = sum(nountypeAZ(:,[1 3]),2)==2;  % find count and collective words
nountypeAZ(ix,1)   = false;          % drop countAndCollective from individual count

% partofspeech mod
partspeechAZ = partofspeech;
partspeechAZnames = {'nounOnly' 'nounAndVerb' 'nounVerbAdj'};
% all words are nouns.
ix = sum(partspeechAZ,2)==3;        % find noun, verb, and adj words
partspeechAZ(ix,1:2) = false;       % drop noun and verb from nounVerbAdj words
ix = sum(partspeechAZ(:,1:2),2)==2; % find noun and verb words
partspeechAZ(ix,1)   = false;       % drop noun and verb from nouns
% dropping adjectives that are not also verbs (3 words)
ix = partspeechAZ(:,3) & sum(partspeechAZ(:,1:2),2)==1;
partspeechAZ(ix,:) = false;

% build category structures
categoryTypes = {'myCategories', 'myCategories1', 'nountypeAZ', 'partspeechAZ', 'wordnet'};
mycats(1)= Categories(mycategories,mycategorynames);
mycats(2)= Categories(mycategory1 ,mycategorynames);
mycats(3)= Categories(nountypeAZ  ,nountypeAZnames);
mycats(4)= Categories(partspeechAZ,partspeechAZnames);
mycats(5)= Categories(wordnetcategories,wordnetcategorynames);

% Deal with words belonging to multiple categories
for i = 1:numel(mycats)
   cix{i} = sum(mycats(i).vectors,2)==1;
   icat = sum(mycats(i).vectors(cix{i},:))>4;
   if ~all(icat)
      mycats(i) = mycats(i).select(icat);
   end
   cix{i} = sum(mycats(i).vectors,2)==1;
   mycats(i).vectors = mycats(i).vectors(cix{i},:);
end

%% LOOP OVER ALL TYPES OF CATEGORIES, COMPUTE h2s
fprintf('Running h2s on regions ');
for r = 1:nr
   stationarycounter(r,nr);
   for i = 1:numel(mycats)
      dataslice = reshape(datan(ir==r,:,cix{i}),[],sum(cix{i}));
      % Run H2S
      reg(r).hi(i)  = SetOfHyps('estimate',dataslice,mycats(i),1);
      reg(r).lo(i)  = reg(r).hi(i).h2s(FAST);
      reg(r).lo(i).sig = reg(r).lo(i).significance(dataslice);
   end
end
fprintf(' done.\n');

%% RENDER h2s RESULTS
for i = 1:numel(mycats)
   fh(i) = figure('Name',categoryTypes{i});
   for r = 1:nr
      ax = subplot(ceil(sqrt(nr)),floor(sqrt(nr)),r);
      reg(r).lo(i).show(false);
      title({regions{r};ax.Title.String})
   end
   ann = mycats(i).legend([0.01 0.01 1 0.2]);
   ann.VerticalAlignment = 'bottom';
end

printFig(fh,[],'png',450)


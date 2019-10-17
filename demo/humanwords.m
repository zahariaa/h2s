% humanwords.m : Loads Laura's data and runs h2s on it
% for different cuts/categories of words

%% LOAD AND RESHAPE DATA
fprintf('Loading...');
load('~/Desktop/humanwords/datan.mat')
load('~/Desktop/humanwords/elec.mat')
load('~/Desktop/humanwords/features.mat')
fprintf(' done.\n');

datan = reshape(permute(datan,[1 3 2]),[],300);

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
mycats(1)= Categories(mycategories,mycategorynames);
mycats(2)= Categories(mycategory1 ,mycategorynames);
mycats(3)= Categories(nountypeAZ  ,nountypeAZnames);
mycats(4)= Categories(partspeechAZ,partspeechAZnames);
mycats(5)= Categories(wordnetcategories,wordnetcategorynames);

%% LOOP OVER ALL TYPES OF CATEGORIES, COMPUTE h2s
for i = 1:5
   % Deal with words belonging to multiple categories
   ix = sum(mycats(i).vectors,2)==1;
   icat = sum(mycats(i).vectors(ix,:))>4;
   if ~all(icat)
      mycats(i) = mycats(i).select(icat);
   end
   ix = sum(mycats(i).vectors,2)==1;
   mycats(i).vectors = mycats(i).vectors(ix,:);

   % Run H2S
   hi(i)  = SetOfHyps('estimate',datan(:,ix),mycats(i),1);
   lo(i)  = hi(i).h2s(false);%(true);
   lo(i).sig = lo(i).significance(datan(:,ix));
end

%% RENDER h2s RESULTS
[legy,legx] = meshgrid([0.65 0.35 0.05],[0.01 0.51]);
figure;
for i = 1:5
   subplot(3,2,i);
   lo(i).show(false);
   mycats(i).legend([legx(i) legy(i) 0.33 0.3])
end


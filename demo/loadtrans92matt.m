% load all matt's trans92 single unit data

subjects = [repmat({'Vortex'},[1 12]) repmat({'Vulcan'},[1 10])];
sessions = {'20180420','20180425','20180426','20180430','20180606','20180607','20180608',...
            '20180618','20180619','20180704','20180807','20180808',...
            '20180420','20180606','20180607','20180608','20180618','20180619','20180625',...
            '20180704','20180705','20180710'};

%% LOAD ALL DATA, skipping intermediate data structures
tun = cellfun(@(subj,sess) utahTuning(mattloader(subj,sess)),subjects,sessions);


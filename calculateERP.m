%% calculateERP
% load data befor and after artifact attenuation
% LPF
% reject invalid button presses: too close to each other
% epoch: +-500ms
% baseline correct: -500 to -300ms
% define channel and time of interest for N1 & movement-related-cortical potentials (MRCP) by inspection of group average
% save maps (centered around peaks)
% estimate signal (average of ROI) and noise (standard deviation of baseline)
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020

PATHINall = {fullfile(PATH, 'derivates', 'rawPrepared'),...
    fullfile(PATH, 'derivates','artifactAttenuated')};

PATHOUTall = fullfile(PATH, 'derivates','specificity');

% parameters
EV = {{'LeftButtonPress'}, {'RightButtonPress'},{'LeftButtonPress','RightButtonPress'}};
COND = {'beforeCleaning', 'afterCleaning'};
REJ         = 3;    % SD for rejection

% define channel and time of interest for N1 & MRCP by inspection of group average
MRCP(1).chanName = {'C2', 'C4', 'CP2', 'CP4'};
MRCP(1).chanIdx  = find(ismember({chanlocs.labels}, MRCP(1).chanName));
MRCP(1).ev       = EV(1);
MRCP(1).time     = [-110 10];
MRCP(1).baseline = [-500 -300];

MRCP(2).chanName = {'C1','C3', 'CP1','CP3'};
MRCP(2).chanIdx  = find(ismember({chanlocs.labels}, MRCP(2).chanName));
MRCP(2).ev       = EV(2);
MRCP(2).time     = MRCP(1).time;
MRCP(2).baseline = MRCP(1).baseline;

N1(1).chanName = {'Fz', 'FC1', 'FC2','Cz'};
N1(1).chanIdx  =  find(ismember({chanlocs.labels}, N1(1).chanName));
N1(1).ev       = EV(3);
N1(1).time     = [150 170];
N1(1).baseline = [-500 -300];

%%
for c = 1:length(PATHINall) % loop through data before and after artifact attnuation
    
    % load rejection overview
    numSubj = readtable(fullfile(PATHINall{c}, 'group', 'subjStats.xls'));
    
    for s = 1:height(participants) % process each subject
        
        PATHOUT = fullfile(PATHOUTall, participants(s,1:end).participant_id);
        if ~exist(PATHOUT,'dir') % create output directories if necessary
            mkdir(PATHOUT)
        end
        
        % get file handles
        if c == 1       % befor artifact attenuation
            FILES = dir(fullfile(PATHINall{c}, participants(s,1:end).participant_id ,'*.set'));
        elseif c == 2   % after artifact attenuation
            FILES = dir(fullfile(PATHINall{c}, participants(s,1:end).participant_id ,'*clean.set'));
        end
        
        % load data _______________________________________________________
        EEG = pop_loadset(FILES.name, FILES.folder);
        
        % interpolate removed channels ____________________________________
        % (only necessary for data befor artifact attenuation,
        % channels were interpolated during artifact attenuation)
        EEG = pop_interp(EEG, chanlocs, 'spherical');
        
        % LPF _____________________________________________________________
        % zero-phase FIR, order: 84, cut-off 45 Hz
        EEG = pop_eegfiltnew(EEG, 'hicutoff',40);
        
        % reject invalid button presses ___________________________________
        % delete button presses that are less than 800ms apart
        idxButton = find(contains({EEG.event.type}, 'ButtonPress')); % index all Button presses (also the missed attempts)
        rejButton = find(diff([EEG.event(idxButton).latency])<.8*EEG.srate); % find those less than .8s apart
        idxRej = idxButton([rejButton,rejButton+1]); % get index of both button presses
        [EEG.event(idxRej).type] = deal('RejButt'); % rename those events -> won't  be used for epoching
        % store information
        numSubj.numButton(s) = length(idxButton); % number of total button presses
        numSubj.numValidButton(s) = length(idxButton)-length(rejButton); % number of button presses further than .8s apart
        
        
        for ev = 1:length(EV) % epoching __________________________________
            
            % epoch 
            EEGtmp = pop_epoch( EEG, EV{ev}, [-.5  .5], 'newname', [FILES.name(1:7), EV{ev}{:},'Epochs_' COND{c}], 'epochinfo', 'yes');
            
            % remove baseline 
            EEGtmp = pop_rmbase( EEGtmp, N1.baseline); % MRCP and N1 baseline are the same
            
            % epoch rejection 
            EEGtmp = pop_jointprob(EEGtmp,1,1:size(EEGtmp.data,1),REJ,REJ,1,1); % probability
            
            if ev ==1
                numSubj.numLeftButton(s) = EEGtmp.trials;
            elseif ev == 2
                numSubj.numRightButton(s) = EEGtmp.trials;
            end
            
            % store data 
            data(:,:,s,c,ev) = mean(EEGtmp.data,3); 
            % chans, pnts, event (lef/right button press), condition (before/after artifact attenuation)
            
            % create group folder, if needed
            if ~exist(fullfile(PATHOUTall, 'group'),'dir')
                mkdir(fullfile(PATHOUTall, 'group'))
            end
            
            % create a dataset for group averages
            if ev ==1 && s ==1 && c==1 % only do this once
                EEGav = pop_select(EEGtmp, 'trial', 1:height(participants));
                EEGav.setname = 'sub-all_';
                pop_saveset(EEGav, EEGav.setname, fullfile(PATHOUTall, 'group'));
                
                [~, from] = min(abs(EEGav.times-MRCP(1).time(1)));
                [~, to] = min(abs(EEGav.times-MRCP(1).time(2)));
                MRCP(1).pnts = from:to;
                MRCP(2).pnts = MRCP(1).pnts;
                
                [~, from] = min(abs(EEGav.times-N1(1).time(1)));
                [~, to] = min(abs(EEGav.times-N1(1).time(2)));
                N1(1).pnts = from:to;
            end
            
            % save %%%
            pop_saveset(EEGtmp, 'filename',EEGtmp.setname,'filepath',PATHOUT);
            
        end
    end
    
    % save
    writetable(numSubj,fullfile(PATHOUTall, 'group', ['subjStats_', COND{c} ,'.xls']));
    save(fullfile(PATHOUTall, 'group',['data_'  COND{c}]), 'data');
end


%% group averages _________________________________________________________
EEG = EEGav;
for c = 1:length(PATHINall)
    for ev = 1:length(EV)
        EEG.data = data(:,:,:,c,ev);                                    % extract data from matrix it was perv. stored in
        EEG.setname = [EEGav.setname, EV{ev}{:},'Epochs_' COND{c}];     % rename set
        pop_saveset(EEG, EEG.setname, fullfile(PATHOUTall, 'group'));   % save
        
        % find points for noise estimate
        [~, from] = min(abs(EEGav.times-MRCP(1).baseline(1)));
        [~, to] = min(abs(EEGav.times-MRCP(1).baseline(2)));
        
        if ev ==3 % N1 ____________________________________________________
            
            % map
            N1.topo(:,:,c) = squeeze(mean(EEG.data(:,N1.pnts,:),2));
            
            % noise estimate
            N1.noisePnts = from:to;
            N1.noise(:,c) = std(mean(EEG.data(N1.chanIdx, N1.noisePnts,:)));
            
            % signal estimate
            N1.signal(:,c) = mean(mean(EEG.data(N1.chanIdx, N1.pnts,:)));
            
        else % MRCP _______________________________________________________
            
            % map
            MRCP(ev).topo(:,:,c) = squeeze(mean(EEG.data(:,MRCP(ev).pnts,:),2));
            
            % noise estimate
            MRCP(ev).noisePnts = from:to;
            MRCP(ev).noise(:,c) = std(mean(EEG.data(MRCP(ev).chanIdx, MRCP(ev).noisePnts,:)));
            
            % signal estimate
            MRCP(ev).signal(:,c) = mean(mean(EEG.data(MRCP(ev).chanIdx, MRCP(ev).pnts,:)));
        end
    end
end
save(fullfile(PATHOUTall, 'group', 'ERPs' ), 'MRCP','N1', 'data');

%% housekeeping
clearvars -except PATH participants chanlocs
close all
clc
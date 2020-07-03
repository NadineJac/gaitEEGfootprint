%% prepareRawData
% 1 Hz HPF
% bad channel rejection
% average Ref
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 29-June-2020

% start EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% create table to store removed channels
numSubj = table('Size',[height(participants),1],'VariableTypes', {'string'},  'VariableNames', {'ID'});
PATHOUT = fullfile(PATH, 'derivates', 'rawPrepared');
PATHOUTacc = fullfile(PATH, 'derivates', 'acc');
if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
    mkdir(PATHOUTacc)
end

%%
for s = 1:height(participants)
    
    PATHIN = fullfile(PATH, participants(s,1:end).participant_id, 'eeg');
    
    FILES = dir([PATHIN filesep '*.set']);
    numSubj.ID(s) = participants(s,1:end).participant_id;
    
    % load data
    EEG = pop_loadset(FILES.name, FILES.folder);
    
    % 1 Hz HPF: 
    % order: 826, cutoff frequency: 1 Hz,(zero-phase, non-causal)
    EEG = pop_eegfiltnew(EEG, 'locutoff',2);
    
    % 135 Hz LPF: 
    % order: 56, cutoff: 135 Hz, (zero-phase, non-causal)
    EEG = pop_eegfiltnew(EEG, 'hicutoff',120);
    
    % downsample to 250 Hz
    EEG = pop_resample( EEG, 250);
    
    % save acceleration data
    Acc = pop_select(EEG, 'nochannel', [1:64]);
    tmp = fullfile(PATHOUTacc, participants(s,1:end).participant_id);
    mkdir(tmp)
    pop_saveset(Acc, [FILES.name(1:end-7), '_acc'], tmp);
    
    % only keep EEG data
    EEG = pop_select(EEG, 'channel', [1:64]);
    
    % bad channel rejection: default parameters of clean_rawdata toolbox
    nbChan = EEG.nbchan;
    EEG = clean_rawdata(EEG, 5, -1, 0.8, 4, -1, -1);
    
    % store number of removed channels
    numSubj.rmChan(s) = nbChan-EEG.nbchan;
    
    % common average re-refrencing
    EEG = pop_reref(EEG, []);
    
    % save
    tmp = fullfile(PATHOUT, participants(s,1:end).participant_id);
    mkdir(tmp)
    pop_saveset(EEG, FILES.name, tmp);
end

tmp = fullfile(PATHOUT, 'group');
mkdir(tmp)
writetable(numSubj,[tmp filesep 'subjStats.xls']);

% housekeeping
clearvars -except PATH participants chanlocs
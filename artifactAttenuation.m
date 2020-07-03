%% artifactAttenuation
% clean w/ ASR: flatlines, drifts, clean_asr w/ designated calib data
% only keep walking data
% AMICA
% fit dipoles
% reject components w/ IC label and dipole RV + location
%
% ATTENTION: adapt eeglab path to run on your machine (specified in main.m)
%
% Developed in MATLAB R2018b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020

% header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% directories
PATHIN = fullfile(PATH, 'derivates','rawPrepared');
PATHOUT = fullfile(PATH, 'derivates','artifactAttenuated');
if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
end

mkdir(fullfile(PATHOUT, 'group')) % create directory for group files
numSubj = readtable(fullfile(PATH, 'derivates', 'rawPrepared', 'subjStats.xls'));

% preprocessing parameters
REJ         = 3;    % SD for rejection of dummy epochs before ICA
CONDS.walk  = {'standing', 'difficult', 'difficult_button', 'easy', 'easy_button'}; % conditions of experiment

% ASR parameters: 
cutoff      = 7;    % repair bursts, rest kept at default (hard coded)

% ICA parameters
numprocs    = 1;    % # of nodes (default = 1)
max_threads = 2;    % # of threads per node
num_models  = 1;    % # of models of mixture ICA
max_iter    = 2000; % max number of learning steps


%% ASR correction and ICA decomposition   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1:height(participants)
    
    FILES = dir(fullfile(PATHIN, participants(s,1:end).participant_id ,'*.set')); 
    SETNAME = FILES.name(1:27);
    
    % load data
    EEG = pop_loadset(FILES.name, FILES.folder);
    
    %%% clean w/asr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract calibration data
    from = round(EEG.event(strcmp({EEG.event.type}, 'start_restEEG')).latency);
    mybaseline = pop_select(EEG, 'point', [from:from+60*EEG.srate]); % extract 1 min baseline recording
    
    % extract data to be cleaned: 
    % remaining data of experiment after snipped extracted for ASR (has to stay continous)
    to = length(EEG.data);
    EEG = pop_select(EEG, 'point', [from+60*EEG.srate:to]);
    data_old = EEG.data; % keep uncorrected data for visualization only
    
    % correct with ASR
    EEG = clean_asr(EEG,cutoff,[],[],[],mybaseline);
    
    % visualize:
    % eegplot(EEG.data, 'data2', data_old);
    
    %%% ICA decomposition and dipole fitting %%%%%%%%%%%%%%%%%%%%
    
    % preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract relevant data: walking data + rest of standing data
    % exclude data from inbetween conditions
    lat = [];
    for c = 1:length(CONDS.walk)
        lat(c,1) = EEG.event(strcmp({EEG.event.type},['start_',CONDS.walk{c}])).latency-1;
        lat(c,2) = EEG.event(strcmp({EEG.event.type},['end_',CONDS.walk{c}])).latency+1;
    end
    EEG = pop_select(EEG,'point', [1 60*EEG.srate;sort(lat)]);
    
    % create pseudo-epochs, 1 second long, unrelated to task structure
    EEGtmp = eeg_regepochs(EEG,1);
    
    % remove improbable epochs
    EEGtmp = pop_jointprob(EEGtmp,1,1:size(EEGtmp.data,1),REJ,REJ,1,1);
    
    % run AMICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mkdir(fullfile(PATHOUT, SETNAME(1:6))) % create output directory
    % attention: AMICA data has to be chans x frames, i.e. not epoched!
    [weights,sphere,mods] = runamica15(EEGtmp.data(:,:),...
        'num_models',num_models,'outdir',fullfile(PATHOUT,SETNAME(1:6),[SETNAME(1:7) 'weights']), ...
        'numprocs', numprocs, 'max_threads', max_threads, 'max_iter',max_iter);
    
    % load computed and stored weights
    EEG = pop_loadmodout(EEG,fullfile(PATHOUT,SETNAME(1:6),[SETNAME(1:7) 'weights']));
    
    % interpolate previously removed channels %%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_interp(EEG, chanlocs, 'spherical');   
       
    % fit dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up head model
     eeglBEM = fullfile(eegl, 'plugins', 'dipfit3.2', 'standard_BEM');
     EEG=pop_chanedit(EEG, 'lookup',fullfile(eeglBEM, 'elec', 'standard_1005.elc'));
     EEG = pop_dipfit_settings( EEG, ...
         'hdmfile', fullfile(eeglBEM, 'standard_vol.mat'),...
         'coordformat','MNI',...
         'mrifile',fullfile(eeglBEM, 'standard_mri.mat'),...
         'chanfile',fullfile(eeglBEM, 'elec','standard_1005.elc'),...
         'coord_transform',[0 0 0 0 0 -1.5708 1 1 1],...
         'chansel',[1:EEG.nbchan] );

    
    % estimate dipoles: 
    % only keep dipoles located insode the brain and w/ residual variance
    % <15 % (does not reject ICs, only removes entry from EEG.dipfit)
    EEG = pop_multifit(EEG, 1:EEG.nbchan ,'threshold',15,'rmout','on'); %did not work in R2019b but consistently in R2018b
    
    % visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot topographies w/ dipole locations
    pop_topoplot(EEG,0, [1:size(EEG.icawinv,2)],EEG.setname,[] ,1,'electrodes','on');
    sgtitle(SETNAME, 'interp', 'none')
    print('-dpng', fullfile(PATHOUT,SETNAME(1:6),[SETNAME(1:7) 'ICA_comps.png'])); % save
    close;
    
    % save
    EEG.setname = SETNAME;
    pop_saveset(EEG, [SETNAME '_decomp'], fullfile(PATHOUT, SETNAME(1:6)));
    

end

%% IC rejection
for s = 1:height(participants)      
    
      % load data
    FILES = dir(fullfile(PATHOUT, participants(s,1:end).participant_id ,'*decomp.set'));
    EEG = pop_loadset(FILES.name, FILES.folder);
    
    % IC label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % classify ICs W/IC label
    EEG = iclabel(EEG);
    
    % remove eye and muscle ICs (classification >.9)
    rmEye = EEG.etc.ic_classification.ICLabel.classifications(:,3)>.9;
    numSubj.IClabelEye(s) = sum(rmEye);
    
    rmMuscle = EEG.etc.ic_classification.ICLabel.classifications(:,2)>.9;
    numSubj.IClabelMuscle(s)= sum(rmMuscle);
    
    % Dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove ICs w/ dipole residual variance > 15 %
    rmDipole =[EEG.dipfit.model.rv]'>.15;
    numSubj.dipfit(s) = sum(rmDipole);
    
    % reject components
    EEG = pop_subcomp( EEG, find(any([rmEye, rmMuscle, rmDipole],2)), 0);
    numSubj.ICremain(s) = size(EEG.icawinv,2);
    
    % save dataset
    pop_saveset(EEG, [FILES.name(1:end-10) 'clean.set'], fullfile(PATHOUT, participants(s,1:end).participant_id));
    
end

writetable(numSubj,fullfile(PATHOUT, 'group','subjStats.xls'));

% housekeeping
clearvars -except PATH participants chanlocs
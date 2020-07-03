%% prepareFootprint
% IN: raw or artifact attenuated continous data
% Out: time-warped ERP/head accelaration and ERSP needed for footprint
% calculation
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 29-June-2020

PATHINall= {fullfile(PATH, 'derivates','rawPrepared');...% rawdata
    fullfile(PATH, 'derivates','artifactAttenuated');... % artifact attenuated data (ASR+ICA)
    fullfile(PATH, 'derivates','artifactAttenuated')}; % only ASR corrected data
PATHINacc= fullfile(PATH, 'derivates','acc');
PATHOUT = fullfile(PATH, 'derivates','footprint');
if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
end

CONDS = {'before', 'after', 'afterASR'};
stand = {'restEEG', 'standing'}; % conditions to extract data from

% epoch and wrping parameters
EP          = [-1 3];               % epoching in seconds
EV          = {'RightHS'};          % event to epoch around
newLat      = [1 16 50 66 100];     % new latencies (to warp to), in pnts
newLat2     = [0 160 500 660 1000]; % new latencies (to warp to), in ms
allGaitEV   = {'LeftTO', 'LeftHS', 'RightTO', 'RightHS', 'RightHS'}; % all gait events that will be used for warping
RESPONSEtime = [25 1500];           % define time window in which resonse (aka next RHS) can occur [ms]

%%
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start eeglab
for c = 1:length(CONDS)
    
    PATHIN = PATHINall{c}; % set input directory
    numEpochsGait = table('Size',[height(participants), 1],'VariableTypes',{'cell'}, 'VariableNames', {'ID'});
    
    for s = 1:height(participants)
        % create spbject's output directory
        PATHOUTtmp = [PATHOUT filesep participants(s,1:end).participant_id];
        mkdir(PATHOUTtmp)
        
        % load EEG data ___________________________________________________
        if c == 1
            FILES = dir(fullfile(PATHIN, participants(s,1:end).participant_id, '*.set')); % only preprocessed
        elseif c == 2
            FILES = dir(fullfile(PATHIN, participants(s,1:end).participant_id, '*clean.set')); % artifact attenuation 1
        elseif c == 3
            FILES = dir(fullfile(PATHIN, participants(s,1:end).participant_id, '*decomp.set')); % artifact attenuation 2
        end
        FILESacc =  dir(fullfile(PATHINacc, participants(s,1:end).participant_id, '*.set')); % accleration
        
        SubjectNames{s} = FILES.name(1:6);
        numEpochsGait.ID(s) = SubjectNames(s);
        
        EEG = pop_loadset(FILES.name, FILES.folder);
        
        % interpolate missing channels
        EEG = pop_interp(EEG, chanlocs, 'spherical');
        
        % extract standing baseline
        if c == 1
            lat=[];
            for ev = 1:length(stand)
                lat(ev,1) = EEG.event(strcmp({EEG.event.type},['start_',stand{ev}])).latency-1;
                lat(ev,2) = EEG.event(strcmp({EEG.event.type},['end_',stand{ev}])).latency+1;
            end
            
        else
            lat = [1, 5*60*EEG.srate]; % length of standing snippet after artifact attenuation
        end
        EEGstand = pop_select(EEG,'point', lat);
        EEGstand.setname =  [FILES.name(1:6), '_standingBL_' CONDS{c}];
        pop_saveset(EEGstand, EEGstand.setname, PATHOUTtmp);
        
        % add accelertion data to EEG structure ___________________________
        if c == 1
            Acc = pop_loadset(FILESacc.name, FILESacc.folder);
            EEG.data = [EEG.data; Acc.data];
            EEG.chanlocs = [EEG.chanlocs, Acc.chanlocs];
            EEG.nbchan = size(EEG.data,1);
            EEG = eeg_checkset(EEG);
        end
        
        
        % epoch ___________________________________________________________
        EEG = pop_epoch(EEG, EV, EP,  'epochinfo', 'yes');
        numEpochsGait.numCycle(s) = EEG.trials;
        
        % only keep valid gait cycles _____________________________________
        % get gait event latencies (ms)
        latGaitEV = zeros(EEG.trials, length(allGaitEV));
        for i = 1:length(allGaitEV)-1
            evalc('latGaitEV(:,i+1) = eeg_getepochevent(EEG ,allGaitEV(i), RESPONSEtime, ''latency'')');
        end
        
        % flag epochs w/o all events
        % contain nans
        rmEp1 =any(isnan(latGaitEV),2);
        
        % and events in wrong order
        % negtive difference between events
        rmEp2 = any(diff(latGaitEV,[],2)<0,2);
        
        % reject flagged epochs
        EEG = pop_select(EEG, 'notrial', find(rmEp1+rmEp2));
        numEpochsGait.numValidCycle(s) = EEG.trials;
        
        % reject epochs with extreme parameters ___________________________
        EEG = pop_jointprob(EEG,1,[1:64] ,3,3,0,1,0,[],0);
        numEpochsGait.numRemainCycle(s) = EEG.trials;
        EEGorg = EEG;
        
        % TF transform and time-normalization _____________________________
        latGaitEV  = nan(EEG.trials,length(allGaitEV)); % preallocate matrix
        for i = 1:length(allGaitEV)-1                   % get latencies of gait events
            evalc('latGaitEV (:,i+1) = eeg_getepochevent(EEG ,allGaitEV(i), RESPONSEtime, ''latency'');'); % in ms, suppress output
        end
        oldLat = round(EEG.srate/1000*latGaitEV);       % convert latencies from ms to pnts (for newtimef)
        
        % time-warp (EEG and accelerometer) -> ERP
        oldLat(:,1) = 1;                                % 1st latency cant be 0, has to be 1
        from        = find(EEG.times == 0);
        
        % ERP _________________________________________
        % warping (resampling would likly be fine too (young sample, regular gait), and quicker!)
        for e = 1:size(EEG.data,3)                      % loop through all epochs
            LAT = oldLat(e,end);                        % get next RHS
            warpmat = timewarp(oldLat(e,:), newLat);    % get warping matrix
            for ch = 1:size(EEG.data,1)                 % loop through channels
                data = squeeze(EEG.data(ch, from:from+LAT-1,e))'; % get data from each epoch
                warped_data(ch,:,e) = warpmat*data;     % warp data
            end
        end
        %         EEG         = pop_resample(EEG, 100); % prepare structure to store data
        %         EEG = pop_select(EEG, 'time', [0,1]);
        %         EEG.data = warped_data;
        %         pop_saveset(EEG, [FILES(s).name(1:11) '_warpedGaitEpochs_' CONDS{c}], PATHOUTtmp);
        
        % average over epochs
        Mdata(:,:,s) = mean(warped_data, 3);
        allOldLat(s,:) = mean(oldLat);                  % store avarage of old latencies
        
        % ERSP _________________________________________
        EEG = pop_select(EEGorg, 'channel', [1:64]);    % only keep EEG
        latGaitEV(:,1) = 0;
        
        for ch = 1:length(EEG.chanlocs) % loop through all channels
            
            % get TF of standing baseline
            [~, ~, ~, ~, ~, ~, ~, tfdatab] = newtimef...
                (EEGstand.data(ch,:), EEGstand.pnts, [EEGstand.xmin EEGstand.xmax]*1000, EEGstand.srate,...
                'cycles', [3 1], 'wletmethod','dftfilt3', 'timesout', 200, 'freqs',[4 60],...
                'plotitc', 'off', 'nfreqs', 57,'plotersp', 'off');
            baseline(ch,:,s) = squeeze(mean(abs(tfdatab).^2,2)); %store mean time
            
            % TF and warp all gait cycles
            [~, ~, ~, times,  freqs, ~, ~, tfdata] = newtimef...
                (EEG.data(ch,:,:), EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate,...
                'cycles', [3 1], 'wletmethod','dftfilt3', 'timesout', [0:10:1000], ...
                'freqs',[4 60],'timewarp', latGaitEV, 'timewarpms', newLat2, ...
                'plotitc', 'off', 'nfreqs', 57, 'plotersp', 'off'); % add 'alpha', 0.05, to show thresholded map, 'baseline',[EEG.xmax],
            
            % store
            TF.data(ch,:,s,:) = mean(abs(tfdata).^2,3);
        end
        
    end
    
    %% save averaged data _________________________________________________
    % store subject averages in one EEGlab dataset
    
    % create group's output directory
    PATHOUTtmp = [PATHOUT filesep 'group'];
    mkdir(PATHOUTtmp)
    
    % create new EEGlab structure w/correct size
    allEEG = pop_resample(EEGorg, 100);
    allEEG = pop_select(allEEG, 'time', [0 1], 'trial', [1:height(participants)]);
    
    % add data
    allEEG.data = Mdata;
    
    % add new events
    allEEG.event = [];allEEG.urevent = []; allEEG.epoch = [];
    for ep = 1:allEEG.trials
        allEEG.epoch(ep).event = [1:length(allGaitEV)];
        allEEG.epoch(ep).eventtype = allGaitEV;
        allEEG.epoch(ep).eventlatency = newLat(1,:);
        for ev = 1:length(allGaitEV)
            allEEG.event(end+1).latency = (ep-1)*allEEG.pnts+newLat(ev);
            allEEG.event(end).type = allGaitEV{ev};
            allEEG.event(end).epoch = ep;
        end
    end
    
    % save
    allEEG = eeg_checkset(allEEG, 'eventconsistency');
    allEEG.setname = ['sub-all_warpedGaitEpochs_', CONDS{c}];
    allEEG.oldLat = allOldLat;
    allEEG.subjectNames = SubjectNames;
    
    % save EEG and acceleration seperatly
    EEG = pop_select(allEEG, 'channel', [1:64]);% EEG
    pop_saveset(EEG, allEEG.setname, PATHOUTtmp);
    
    if c == 1 %acceleration
        Acc = pop_select(allEEG, 'nochannel', [1:64]);
        pop_saveset(Acc, 'sub_all_warpedGaitEpochs_acc', PATHOUTtmp);
    end
    
    %%% TF %%%%%
    TF.times            = times;
    TF.chanlocs         = EEG.chanlocs;
    TF.freqs            = freqs;
    TF.baseline         = baseline;
    TF.events.latency   = newLat2;
    TF.events.labels    = allGaitEV;
    TF.subjectNames     = SubjectNames;
    
    % save data
    save([PATHOUTtmp,filesep, 'TF_' CONDS{c}], 'TF');
    writetable(numEpochsGait, [PATHOUTtmp,filesep, 'numEpochsGait_' CONDS{c}]);
end

% housekeeping
clearvars -except PATH participants chanlocs
close all
clc
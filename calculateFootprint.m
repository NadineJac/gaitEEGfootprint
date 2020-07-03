%% calculateFootprint
% with the prepared data calculate the footprint and calculate the
% distances
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020

% directories
PATHIN = fullfile(PATH, 'derivates','footprint','group');
PATHOUT = fullfile(PATHIN, 'results');
if ~exist(PATHOUT,'dir') % create output directory if necessary
    mkdir(PATHOUT)
end
CONDS = {'before', 'after', 'afterASR'};

% params
% channel indices (dim 1 of the time-frequency decomposed data)
lateralChanIdx = [1,4,5,6,9,10,15,16,17,18,20,21,26,27,30,31,32,33,35,37,41,44,45,46,48,49,51,52,53,57,61,64]; % index of channels labelled as central
neckChanR = [49 52 18 51];  % index of channels located over the right side of the neck
neckChanL = [45 48 46 16];  % index of channels located over the left side of the neck

% sample indices 
% [dim 4 (dim 3 after extraction of single-subject data) of the time-frequency decomposed data]
pntsRHS    = 1:16;          % double support following right-heel strike
pntsLHS    = 51:66;         % double support following left-heel strike
pntsDouble = [pntsRHS, pntsRHS];

%%
for c = 1:length(CONDS)
    
    % set up table to store features (only 1st time)
    tmp =[PATHOUT filesep 'gait_footprint_' CONDS{c} '.mat'];
    if ~isfile(tmp)
        gaitFootprint = table(participants.participant_id, 'VariableNames',{ 'ID'});
        save(tmp, 'gaitFootprint');
    else
      load([PATHOUT filesep 'gait_footprint_' CONDS{c}], 'gaitFootprint');      % load features
    end
 
    % load data
    EEG  = pop_loadset(['sub-all_warpedGaitEpochs_', CONDS{c} '.set'], PATHIN); % gait ERPs
    Acc  = pop_loadset('sub_all_warpedGaitEpochs_acc.set', PATHIN);             % mean gait acceleration
    load([PATHIN filesep 'TF_' CONDS{c}]);                                      % time-frequency decomposed EEG
    load([PATHIN, filesep,'sourceActivations_' CONDS{c}], 'M1', 'wholeBrain')   % source activations
        
    for s = 1:height(participants) % calculate for each subject individually
        
        % extract single-subject data
        ERP = EEG.data(:,:,s);
        ACC = Acc.data([1:3],:,s);
        TFdata = squeeze(TF.data(:,:,s,:));
        TFbaseline = TF.baseline(:,:,s);
        
        % A) R² gfp rms (head acc) -------------------------------------
        gaitFootprint.R2GFPacc(s) = RgfpHeadAcc(ERP, ACC);
        
        % B) power ratio lateral/all channels -------------------------
        [gaitFootprint.lateralPowRatio(s), TF.StandBL] = lateralPowRatio(TFdata, TFbaseline, lateralChanIdx);
        
        % C) correlation across frequencies --------------------------
        gaitFootprint.freqCorr(s) = Rfreq(TFdata, TFbaseline);
        
        % D) power double support/whole gait cycle power -------------
        gaitFootprint.DoubleSuppRatio(s) = doubleSuppRatio(TFdata, TFbaseline, pntsDouble);
        
        % E) power contralateral to all neck electrodes --------------
        gaitFootprint.NeckChanRatio(s) = neckChanRatio(TFdata, TFbaseline, neckChanL, neckChanR, pntsLHS, pntsRHS);
        
        % F) 1-S/W power ratio --------------------------------------
        gaitFootprint.SW(s) = swRatio(TFdata, TFbaseline);
        
        % G) source based metric ------------------------------------
        gaitFootprint.M1(s) = sourceRatio(M1.Value(:,s), wholeBrain.Value(:,s));
    end
    
    % save
    save([PATHOUT,filesep, 'gait_footprint_' CONDS{c}], 'gaitFootprint');
    writetable(gaitFootprint, [PATHOUT,filesep, 'gait_footprint_' CONDS{c}]);
end

%% housekeeping
clearvars -except PATH participants chanlocs
close all
clc
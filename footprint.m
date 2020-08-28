function [feature, tfstand, tfwalk, EEGwalkEP, AccEP] = footprint(EEGstand, EEGwalk, Acc, params)
%% footprint.m
% prepare and calculate gait artifact related footprint features A to F
%
% Syntax:  [feature, tfstand, tfwalk, EEGwalkEP, AccEP] = footprint(EEGstand, EEGwalk, Acc, params)
%
% Inputs:
%EEGstand, EEGwalk, Acc, params
%   EEGstand        - continous EEG data of standing baseline (EEGLAB dataset)
%   EEGwalk         - continous EEG data of walking conditions(EEGLAB dataset)
%   Acc             - continous head acceleration data of walking conditions
%                   (EEGLAB dataset, same samples as EEGwalk)
%   params          - structure with the following parameters:
%         .EV           - event type to epoch around (cell)
%         .allGaitEV    - correct order of all gait event types that will be used for warping (cell)
%         .EP           - epoch start and end (seconds)[vector]
%         .RESPONSEtime - define time window (start, end) in which resonse (aka next RHS) can occur (ms)[vector]
%         .newLatpnts   - new latencies (to warp to)(pnts)[vector, same number of entries as allGaitEV)
%         .newLatms     - new latencies (to warp to)(ms)[vector, s.a.]
%         .newtimef.cycles  - newtimef() parameters
%         .newtimef.wletmethod
%         .newtimef.timesout
%         .newtimef.freqs
%         .newtimef.nfreqs
%         .B.lateralChanIdx - indices of channels labelled as lateral [vector]
%         .E.neckChanR      - indices of channels located over the right side of the neck [vector]
%         .E.neckChanL      - indices of channels located over the left side of the neck [vector]
%         .E.pntsRHS        - sample indices of double support following right-heel strike [vector]
%         .E.pntsLHS        - sample indices double support following left-heel strike [vector]
%         .D.pntsDouble     - sample indices double support [vector]
%
% Outputs:
% feature   - feature values of feature A to F [vector]
% tfstand   - time-frequency transformed EEG baseline (i.e standing), chan x freq
% tfwalk    - time-frequency transformed EEG (averaged and time-normalized to the gait cycle), chan x freq x pnts
% EEGwalkEP - EEGLAB dataset with time-warped gait epochs
% AccEP     - EEGLAB dataset with time-warped head accelaration of gait epochs
%
% Other m-files required: RgfpHeadAcc.m, lateralPowRatio.m, Rfreq.m, doubleSuppRatio.m, neckChanRatio.m, swRatio.m
% Subfunctions: epoch_gait, warp_ERP, gaitERSP
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% August 2020; Last revision: 30-August-2020

%------------- BEGIN CODE --------------
%% check inputs
% datasets
% EEGLAB datasets
eeg_checkset(EEGstand);
eeg_checkset(EEGwalk);
eeg_checkset(Acc),

% same number of channels
if ~(EEGstand.nbchan == EEGwalk.nbchan)
    error('Standing and walking EEG do not have the same number of channels'); end

% same number of points
if ~(Acc.pntsn == EEGwalk.pnts)
    error('Walking EEG and head acceleration do not have the same number of samples'); end

% params: all fields thete (fild themselves will be checked during feature
% calculation)
neededParams = {'EV','allGaitEV','EP','RESPONSEtime',...
    'newLatpnts','newLatms','newtimef',...
    'B', 'E', 'D'};
neededParamsE = {'neckChanR','neckChanL','pntsRHS','pntsLHS'};
neededParamsNewtimef = {'cycles','wletmethod','timesout','freqs', 'nfreqs'};

for i = 1:length(neededParams)
    if~isfield(params, neededParams{i})
        error(['Structure params does not contain the field ', neededParams{i}]); end; end
for i = 1:length(neededParamsNewtimef)
    if~isfield(params.newtimef, neededParamsNewtimef{i})
        error(['Structure params.newtimef does not contain the field ', neededParamsNewtimef{i}]); end; end
for i = 1:length(neededParamsE)
    if~isfield(params.E, neededParamsE{i})
        error(['Structure params.E does not contain the field ', neededParamsE{i}]); end; end
if~isfield(params.B, 'lateralChanIdx')
    error('Structure params.E does not contain the field lateralChanIdx'); end
if~isfield(params.D, 'pntsDouble')
    error('Structure params.E does not contain the field pntsDouble'); end


%% start processing -------------------------------------------------------
EEG = EEGwalk; % combine EEG and acceleration data --> only process once
if isempty(Acc)
    a = 0;
    disp('No head motion data available, skip calculation of feature A')
else
    a =1;
    EEG.data = [EEG.data; Acc.data];
    EEG.chanlocs = [EEG.chanlocs, Acc.chanlocs];
    EEG.nbchan = size(EEG.data,1);
    evalc('EEG = eeg_checkset(EEG)');
end

[EEGEP,latGaitEV] = epoch_gait(EEG, params);

% reject epochs with extreme parameters ___________________________
EEGEP = pop_jointprob(EEGEP,1,[1:EEGwalk.nbchan] ,3,3,0,1,0,[],0);

% get time-warped epochs from continous data
EEGwarped = warp_ERP(EEG, EEGEP,latGaitEV, params);
evalc('EEGwalkEP = pop_select(EEGwarped, ''channel'', [1:EEGwalk.nbchan])'); %delete acceleration cnahhels, return this set

if a % if present return the warped acceleration
    evalc('AccEP = pop_select(EEGwarped, ''nochannel'', [1:EEGwalk.nbchan])');
    evalc('EEGEP = pop_select(EEGEP, ''channel'', [1:EEGwalk.nbchan])');
end

% time frequency decomposition & warping
evalc('[tfstand, tfwalk] = gaitERSP(EEGstand, EEGEP, latGaitEV, params)');

%% extract data for feature calculation
ERP = mean(EEGwalkEP.data,3);
ACC = mean(AccEP.data([1:3],:,:),3);
TFdata = squeeze(tfwalk);
TFbaseline = tfstand;

% A) R² gfp rms (head acc) -------------------------------------
feature(1) = RgfpHeadAcc(ERP, ACC);

% B) power ratio lateral/all channels -------------------------
feature(2)= lateralPowRatio(TFdata, TFbaseline, params.B.lateralChanIdx);

% C) correlation across frequencies --------------------------
feature(3) = Rfreq(TFdata, TFbaseline);

% D) power double support/whole gait cycle power -------------
feature(4) = doubleSuppRatio(TFdata, TFbaseline, params.D.pntsDouble);

% E) power contralateral to all neck electrodes --------------
feature(5) = neckChanRatio(TFdata, TFbaseline, params.E.neckChanL, params.E.neckChanR, params.E.pntsLHS, params.E.pntsRHS);

% F) 1-S/W power ratio --------------------------------------
feature(6) = swRatio(TFdata, TFbaseline);
end
-------------------------------------------------------------------------
-------------------------------------------------------------------------
function [EEGEP, latGaitEV] = epoch_gait(EEG, params)
% epoch around params.EV and store latencies of allgaitEV for later warping
% no input check!

% epoch ___________________________________________________________
EEGEP = pop_epoch(EEG, params.EV, params.EP,  'epochinfo', 'yes');
% numEpochsGait.numCycle(s) = EEGwalk.trials;


% only keep valid gait cycles _____________________________________
% get gait event latencies (ms)
latGaitEV = zeros(EEGEP.trials, length(params.allGaitEV));
for i = 1:length(params.allGaitEV)-1
    evalc('latGaitEV(:,i+1) = eeg_getepochevent(EEGEP ,params.allGaitEV(i), params.RESPONSEtime, ''latency'')');
end

% flag epochs w/o all events
% contain nans
rmEp1 =any(isnan(latGaitEV),2);

% and events in wrong order
% negtive difference between events
rmEp2 = any(diff(latGaitEV,[],2)<0,2);

% reject flagged epochs
EEGEP = pop_select(EEGEP, 'notrial', find(rmEp1+rmEp2));

latGaitEV = latGaitEV(~(rmEp1+rmEp2),:);
end

function EEGwarped  = warp_ERP(EEG, EEGEP, latGaitEV, params)
% extract epochs from continous data and warp all trials to same length (params.newLatpnts) using events of
% params.allGaitEV
% no input check!

% create new EEGlab structure w/correct size
evalc('EEGwarped = pop_resample(EEGEP, 100)');
evalc('EEGwarped = pop_select(EEGwarped, ''time'', [0 1])');

% time-warp (EEG and accelerometer) -> ERP
oldLat = round(EEG.srate/1000*latGaitEV);       % convert latencies from ms to pnts (for newtimef)
oldLat(:,1) = 1;                                % 1st latency cant be 0, has to be 1
from        = find(EEG.times == 0);

% ERP _________________________________________
% warping (resampling is quicker and would work for regular striges, e.g. treadmill)
for e = 1:size(EEG.data,3)                      % loop through all epochs
    LAT = oldLat(e,end);                        % get next RHS
    warpmat = timewarp(oldLat(e,:), params.newLatpnts);    % get warping matrix
    for ch = 1:size(EEG.data,1)                 % loop through channels
        data = squeeze(EEG.data(ch, from:from+LAT-1,e))'; % get data from each epoch
        EEGwarped.data(ch,:,e) = warpmat*data;     % warp data
    end
end
%         EEG         = pop_resample(EEG, 100); % prepare structure to store data
%         EEG = pop_select(EEG, 'time', [0,1]);
%         EEG.data = warped_data;
%         pop_saveset(EEG, [FILES(s).name(1:11) '_warpedGaitEpochs_' CONDS{c}], PATHOUTtmp);

end

function [tfstand, tfwalk, times, freqs] = gaitERSP(EEGstand, EEGEP, latGaitEV, params)
% time-frequency decompose standing baselinge and gait epochs using newtimef
% no input check!
latGaitEV(:,1) = 0;
disp(['Starting time-frequency decomposition of ' num2str(EEGstand.nbchan), ' channels'])
for ch = 1:EEGstand.nbchan % loop through all channels
    
    % get TF of standing baseline
    [~, ~, ~, ~, ~, ~, ~, tfdatab] = newtimef...
        (EEGstand.data(ch,:), EEGstand.pnts, [EEGstand.xmin EEGstand.xmax]*1000, EEGstand.srate,...
        'cycles', params.newtimef.cycles,...
        'wletmethod',params.newtimef.wletmethod,...
        'timesout', params.newtimef.timesout, ...
        'freqs',params.newtimef.freqs,...
        'nfreqs', params.newtimef.nfreqs, ...
        'plotitc', 'off', 'plotersp', 'off');
    tfstand(ch,:) = squeeze(mean(abs(tfdatab).^2,2)); %store mean time
    
    % TF and warp all gait cycles
    [~, ~, ~, ~, ~, ~, ~, tfdata] = newtimef(EEGEP.data(ch,:,:),...
        EEGEP.pnts, [EEGEP.xmin EEGEP.xmax]*1000, EEGEP.srate,...
        'cycles', params.newtimef.cycles,...
        'wletmethod',params.newtimef.wletmethod,...
        'timesout', params.newtimef.timesout, ...
        'freqs',params.newtimef.freqs,...
        'nfreqs', params.newtimef.nfreqs, ...
        'timewarp', latGaitEV,...
        'timewarpms', params.newLatms, ...
        'plotitc', 'off', 'plotersp', 'off'); % add 'alpha', 0.05, to show thresholded map, 'baseline',[EEG.xmax],
    
    % store
    tfwalk(ch,:,1,:) = mean(abs(tfdata).^2,3);
end
end

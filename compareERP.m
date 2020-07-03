%% compareERP
% evaluation of button press MRCP and auditory N1 before and after artifact
% attenuation and store results as .csv for analysis in R, JASP etc
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020

PATHIN = fullfile(PATH, 'derivates','specificity', 'group');
PATHOUT = [PATHIN filesep 'results'];
if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
end

% load data
load([PATHIN,filesep, 'ERPs']);
EEG = pop_loadset('sub-all_.set', PATHIN); % use timepoints for extraction of data

%% ERP similarity _________________________________________________________
% correlate whole epoch before and after artifact attenuation

[~, from] = min(abs(EEG.times+100));
[~, to] = min(abs(EEG.times-300));

for s = 1:size(N1.signal,1)
    for ev = 1:2
        r = corr(mean(data(MRCP(ev).chanIdx,from:to,s,1))', mean(data(MRCP(ev).chanIdx,from:to,s,2))', 'type', 'Pearson');
        MRCP(ev).erpR(s,1) = r;
    end
    
    r = corr(mean(data(N1.chanIdx,from:to,s,1))',mean(data(N1.chanIdx,from:to,s,2))','type', 'Pearson');
    N1.erpR(s,1) = r;
end

% fisher z transform correlations
N1.erpZ =  log((1+N1.erpR)./(1-N1.erpR))/2; % use atanh() instead
MRCP(1).erpZ =  log((1+MRCP(1).erpR)./(1-MRCP(1).erpR))/2;
MRCP(2).erpZ =  log((1+MRCP(2).erpR)./(1-MRCP(2).erpR))/2;

% save results for analysis in JASP or R
erpZ = table(N1.erpZ, MRCP(1).erpZ, MRCP(2).erpZ, ...
    'VariableNames', {'N1', 'MRCP_left', 'MRCP_right'},...
    'RowNames', string(participants.participant_id(:,1:6)));
writetable(erpZ, [PATHOUT filesep 'erpZ']);

% add mean and SD and save as excel file
erpZ{end+1,:} = mean(erpZ.Variables,1);
erpZ{end+1,:} = std(erpZ.Variables);
erpZ.Properties.RowNames(end-1:end) = {'M', 'SD'};
writetable(erpZ, [PATHOUT filesep 'erpZ.xls'],'WriteRowNames',true);

%% SNR ____________________________________________________________________

for c=1:2
    for ev = 1:2
        MRCP(ev).SNRdB(:,c) = 10*log10(abs(MRCP(ev).signal(:,c)./MRCP(ev).noise(:,c)));
    end
    N1.SNRdB(:,c) = 10*log10(abs(N1.signal(:,c)./N1.noise(:,c)));
end

% save results for analysis in JASP or R
SNRdB = table(N1.SNRdB, MRCP(1).SNRdB, MRCP(2).SNRdB,  ...
    'VariableNames', {'N1', 'MRCP_left', 'MRCP_right'},...
    'RowNames', string(participants.participant_id(:,1:6)));
writetable(SNRdB, [PATHOUT filesep 'SNRdB']);

% add mean and SD and save as excel file
SNRdB{end+1,:} = mean(SNRdB.Variables,1);
SNRdB{end+1,:} = std(SNRdB.Variables);
SNRdB.Properties.RowNames(end-1:end) = {'M', 'SD'};
writetable(SNRdB, [PATHOUT filesep 'SNRdB.xls'],'WriteRowNames',true);
%% map similarity _________________________________________________________
% R² of peak ERP maps 

for s = 1:size(N1.signal,1)
    for ev = 1:2
        r = corr(MRCP(ev).topo(:,s,1),MRCP(ev).topo(:,s,2), 'type', 'Pearson');
        MRCP(ev).topoR(s,1) = r;
    end
    
    r = corr(N1.topo(:,s,1),N1.topo(:,s,2),'type', 'Pearson');
    N1.topoR(s,1) = r;
end

% fisher z transform correlations
N1.topoZ =  log((1+N1.topoR)./(1-N1.topoR))/2;
MRCP(1).topoZ =  log((1+MRCP(1).topoR)./(1-MRCP(1).topoR))/2;
MRCP(2).topoZ =  log((1+MRCP(2).topoR)./(1-MRCP(2).topoR))/2;

% save results for analysis in JASP or R
topoZ = table(N1.topoZ, MRCP(1).topoZ, MRCP(2).topoZ,  ...
    'VariableNames', {'N1', 'MRCP_left', 'MRCP_right'},...
    'RowNames', string(participants.participant_id(:,1:6)));
writetable(topoZ, [PATHOUT filesep 'topoZ']);

% add mean and SD and save as excel file
topoZ{end+1,:} = mean(topoZ.Variables,1);
topoZ{end+1,:} = std(topoZ.Variables);
topoZ.Properties.RowNames(end-1:end) = {'M', 'SD'};
writetable(topoZ, [PATHOUT filesep 'topoZ.xls'],'WriteRowNames',true);

%% housekeeping
save([PATHOUT filesep 'ERPresults'],'N1', 'MRCP');
clearvars -except PATH participants chanlocs
clc
%% footprint distances
% calculate eucidean distances between footprint feature vectors of
% different processing stages
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020

% directories
PATHIN = fullfile(PATH, 'derivates','footprint','group','results');
PATHOUT = PATHIN;
CONDS = {'before', 'after', 'afterASR'};

% load footprints
load([PATHIN, filesep, 'gait_footprint_before']);
fingerprintBefore = gaitFootprint;

load([PATHIN, filesep, 'gait_footprint_afterASR']);
fingerprintAfterASR = gaitFootprint;

load([PATHIN, filesep, 'gait_footprint_after']);
fingerprintAfter = gaitFootprint;


%% calculate distance of feature vectors
distFootprint = fingerprintBefore(:,1);
for s = 1:height(fingerprintBefore)
    distFootprint.raw2ASR(s) = norm(fingerprintAfterASR{s,2:end}-fingerprintBefore{s,2:end});
    distFootprint.raw2ICA(s) = norm(fingerprintAfter{s,2:end}-fingerprintBefore{s,2:end});
    distFootprint.ASR2ICA(s) = norm(fingerprintAfter{s,2:end}-fingerprintAfterASR{s,2:end});
end
writetable(distFootprint, [PATHOUT filesep 'footprintDistances']);

%% housekeeping
clearvars -except PATH participants chanlocs
clc
%%  main
% Analyze a mobile EEG dataset of free outdoor walking [LINK DATA here] 
% to obtain a gait-related artifact footprint and test sensitivity and 
% specificity of one artifact attenuation pipeline as suggested in [DOI here]
%
% Developed in MATLAB R2018b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020
%
% Dependencies (add to path previously):
% - eeglab (v14_1_2b) w/ following plugins (can be added via the eeglab GUI)
%   download from: https://sccn.ucsd.edu/eeglab/download.php
%   - AMICA (v1.5.1)
%   - IClabel (v1.1)
%   - clean_rawdata (v1.00)
%   - dipfit (v3.2)
% - brainstorm (vJan-2020)
%   download from https://neuroimage.usc.edu/bst/download.php
% other versions should work as well but were not tested


% clear workspace
clc; clear; close all;

% set directories depending on workstation

%%% data %%%
% PATH = fullfile('E:','nadine','test_footprint_BIDS'); % laptop
PATH = 'D:\DATA_D\test_footprint_scripts_PC\'; % PC

%%% scripts %%%
cd([PATH filesep 'code']);

%%% load info %%%
nSubj = 3;
participants = import_pb_info(fullfile(PATH, 'participants.tsv'), 2:nSubj+1); % skip header
load(fullfile(PATH, 'code', 'chanlocs')); % channel locations

%% preprocessing __________________________________________________________
% 1 Hz HPF
% bad channel rejection
% average Ref
run prepareRawData

%% artifact attenuation ___________________________________________________
% clean w/ ASR: flatlines, drifts, clean_asr w/ designated calib data
% only keep walking data
% AMICA
% fit dipoles
% reject components w/ IC label and dipole RV + location

% ATTENTION: adapt eeglab path to run on your machine (handed over to dipfit
% plugin during artifact attenuation)
%eegl = fullfile('C:','Users','nadine','Documents','MATLAB','eeglab14_1_2b'); %PC
eegl = fullfile('E:','nadine','MATLAB','eeglab14_1_2b');%Laptop

run artifactAttenuation

%% specificity analysis ___________________________________________________
% ERPs (MRCP and N1) before and after artifact attenuation
% (interpolate channels)
% 45 Hz LPF
% reject invalid button presses
% epoch around button presses
% reject epochs with (residual) artifacts
% save subject averages
run calculateERP

% compare them and save results for analysis in R, JASP etc
% similariity of ERP time-course, maps and SNR before and after artifact
% attenuation
run compareERP
% statistical analysis of results performed in R with "statsSpecificity.R"

%% sensitivity analysis: gait-artifact footprint __________________________
% preprocessing:
%     epoch around gait cycle
%     only keep valid steps
%     time-frequency transform & time-warp steps to same length
run prepareFootprint
run estimateSourceActivity 
% brainstorm w/correctly setup protocols needed
% script will prompt you once to load custom scouts (details in script)

% footprintprint calculation and distances between footprints
run calculateFootprint

% calculate euclidean distances between all footprints
run calcFootprintDistances
% statistical analysis of results performed in R with "statsSensitivity.R"

% plot footprint
run plotFootprint 
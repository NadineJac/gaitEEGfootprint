%% plotFootprint
%%% adapt from where data is loaded! %%%%___________________
% plot footprint before and after artifact attenuation
% Calculate distance between feature vectors and illustrate them
% ToDo: decide on metrics order and change it in the table --> plot again!
%
%
% Developed in MATLAB R2019b
% Nadine Jacobsen (nadine.jacobsen@uni-oldenburg.de), 
% March 2020, last revision: 30-June-2020


% directories
PATHIN = fullfile(PATH, 'derivates','footprint','group', 'results');
PATHOUT = PATHIN;
if ~exist(PATHOUT, 'dir')
    mkdir(PATHOUT)
end
CONDS = {'before', 'after', 'afterASR'};

% load footprints
load([PATHIN, filesep, 'gait_footprint_before']);
fingerprintBefore = gaitFootprint;

load([PATHIN, filesep, 'gait_footprint_afterASR']);
fingerprintAfterASR = gaitFootprint;

load([PATHIN, filesep, 'gait_footprint_after']);
fingerprintAfter = gaitFootprint;

% radar plot params
%  get radar plot from https://de.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
% radar plot parameters 
axes_interval = 5;% Axes properties
axes_precision = 2;
axes_display = 'one';
marker_type = 'none';
axes_font_size = 10;
label_font_size = 12;
axes_labels = {'A','B','C','D','E','F','G'}; % Axes labels
fill_option = 'on';
fill_transparency = 1;
orange = [.9 .6 0];
blue = [0 .45 .7];
lightorange = [251 240 217]/255;
lightblue = [217 234 244]/255;
pink = [.8 .6 .7];
lightpink = pink*1.25;
edgeColors = [orange;pink;blue];
colors = [lightorange;lightpink;lightblue];
line_width = 1;
line_style ={'--',':','-'};

% extract data: averace across subjects
P = [mean(fingerprintBefore{:,2:end});mean(fingerprintAfterASR{:,2:end}); mean(fingerprintAfter{:,2:end})];% extract data (obs x vars)
axes_limits = repmat([-.25; 1],1,size(P,2));

%% Spider plot (customized, script in "code")
%  get radar plot from https://de.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
spider_plot_nj(P,...
    'AxesLabels', axes_labels,...
    'AxesInterval', axes_interval,...
    'AxesPrecision', axes_precision,...
    'AxesDisplay', axes_display,...
    'AxesLimits', axes_limits,...
    'FillOption', fill_option,...
    'FillTransparency', fill_transparency,...%        'Color', colors,...
    'LineWidth', line_width,...
    'LineStyle', line_style,...
    'Marker', marker_type,...
    'Color', colors,...
    'edgecolor', edgeColors,...
    'AxesFontSize', axes_font_size,...
    'LabelFontSize', label_font_size);
legend1 = legend({'before', 'stage 1', 'stage 2'}, 'Fontsize', label_font_size);
set(legend1, 'Position',[0.45 0.059 0.18 0.06]);
legend('boxoff')
title(legend1, 'artifact attenuation')
set(gcf, 'units', 'centimeters', 'Position',[0 0 10 15])
sgtitle('Mean gait artifact footrint')

savefig([PATHOUT filesep 'footprint_sensitivity_radarPlot'])
% print([PATHOUT filesep 'footprint_sensitivity_radarPlot'],'-dpng')
close

%% housekeeping
clearvars -except PATH participants chanlocs
clc
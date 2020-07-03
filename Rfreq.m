function [R, standBL] = Rfreq(TFdata, TFbaseline)
%% RgfpHeadAcc - calculate variance of head acceleration that is explained by global field potential of gaitERP
% variance of any frequency time-course explained by other frequnency time
% courses
%
% Syntax:  [R, standBL] = Rfreq(TFdata, TFbaseline)
%
% Inputs:
%   TFdata          - time-frequency transformed EEG (averaged and time-normalized to the gait cycle), chan x freq x pnts
%   TFbaseline      - time-frequency transformed EEG baseline (i.e standing), chan x freq
%   
%
% Outputs:
%   R               - squared correlation
%   standBL         - baseline corrected time-frequency data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% May 2020; Last revision: 12-May-2020

%------------- BEGIN CODE --------------
% check inputs
% at least 2, numeric vectors, third a number, set default lag to 5
narginchk(2,2)

if ~all([size(TFdata, 1), size(TFdata,2)] == size(TFbaseline))
    error('data and baseline matrix do not have the same number of rows and/or columns');
end

% calculation
standBL = TFdata./repmat(TFbaseline,1,1,size(TFdata,3)); % standing baseline correction
dat = squeeze(mean(standBL));               %average over channels
Z = atanh(corr(dat','type', 'Pearson'));    % z transformed pearson correlation
triuZ = triu(Z,1);                          % only keep upper triangle of correlation matrix (w/o diagonal)
meanZ = sum(triuZ,'all')/sum(triuZ~=0, 'all');% average all nonzero elements
meanR = tanh(meanZ);                        % transform back to r
R = meanR.^2;                               % store squared mean correlation --> coefficient of Determination

%------------- END OF CODE --------------
end
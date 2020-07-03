function [suppRatio, standBL] = doubleSuppRatio(TFdata, TFbaseline, pntsDouble)
%% doubleSuppRatio - calculate power during double support to summed power over the whole gait cycle
% across all channels (compared to baseline power)
%
% Syntax:  [suppRatio, standBL] = doubleSuppRatio(TFdata, TFbaseline, pntsDouble)
%
% Inputs:
%   TFdata          - [matrix,chan x freq x pnts] time-frequency transformed EEG (averaged and time-normalized to the gait cycle)
%   TFbaseline      - [matrix, chan x freq] time-frequency transformed EEG baseline (i.e standing),
%   pntsDouble      - [vector] indices of samples during double support
%
%
% Outputs:
%   suppRatio       - [double] power during double support to summed power over the whole gait cycle
%   standBL         - [matrix,chan x freq x pnts] sbaseline corrected time-frequency data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% May 2020; Last revision: 12-May-2020

%------------- BEGIN CODE --------------
% check inputs
narginchk(3,3)

if ~all([size(TFdata, 1), size(TFdata,2)] == size(TFbaseline))
    error('data and baseline matrix do not have the same number of rows and/or columns');
elseif any(~ismember(pntsDouble, 1:size(TFdata,3)))
    error('Some sample indices exeed the numer time points in the EEG data (dimension 3)');
end

% calculation
standBL = TFdata./repmat(TFbaseline,1,1,size(TFdata,3)); % standing baseline correction
datM = squeeze(sum(sum(standBL,1),2)); % sum over chan and frequency -> subj x pnts

% extract data and sum over time
datDoubleSupp = sum(datM(pntsDouble));

suppRatio = datDoubleSupp./sum(datM);

%------------- END OF CODE --------------
end
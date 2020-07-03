function [powRatio, standBL] = lateralPowRatio(TFdata, TFbaseline, lateralChanIdx)
%% lateralPowRatio ADD INFO HERE!
%
% Syntax:  [powRatio, standBL] = lateralPowRatio(TFdata, TFbaseline, lateralChanIdx)
%
% Inputs:
%   TFdata          - time-frequency transformed EEG (averaged and time-normalized to the gait cycle), chan x freq x pnts
%   TFbaseline      - time-frequency transformed EEG baseline (i.e standing), chan x freq
%   lateralChanIdX  - index of channels labelled as lateral (remaining
%   channels will be treated as central)
% 
% Outputs:
%   powRatio       - power ratio of laterall to all chennels
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
narginchk(3,3)

if ~all([size(TFdata, 1), size(TFdata,2)] == size(TFbaseline))
    error('data and baseline matrix do not have the same number of rows and/or columns');
elseif any(~ismember(lateralChanIdx, 1:size(TFdata,1)))
    error('Some central channel indices exeed the numer of channels in the EEg data (dimension 1)');
end

% calculation
standBL = TFdata./repmat(TFbaseline,1,1,size(TFdata,3)); % standing baseline correction
datM = sum(sum(standBL,2),3); % sum over time and frequency -> chan

% calculate and save power ratio
powRatio = sum(datM(lateralChanIdx))./sum(datM);

%------------- END OF CODE --------------
end
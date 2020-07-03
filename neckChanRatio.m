function [neckChRatio, standBL] = neckChanRatio(TFdata, TFbaseline, neckChanL, neckChanR, pntsLHS, pntsRHS)
%% neckChanRatio - ADD HERE!
% across all channels (compared to baseline power)
%
% Syntax:  [neckChRatio, standBL] = neckChanRatio(TFdata, TFbaseline, neckChanL, neckChanR, pntsLHS, pntsRHS)
%
% Inputs:
%   TFdata          - [matrix,chan x freq x pnts] time-frequency transformed EEG (averaged and time-normalized to the gait cycle)
%   TFbaseline      - [matrix, chan x freq] time-frequency transformed EEG baseline (i.e standing),
%   neckChanL       - [vector] indices channels located over the left neck
%   neckChanR       - [vector] indices channels located over the left neck
%   pntsLHS         - [vector] indices of samples following LHS (double support)
%   pntsRHS         - [vector] indices of samples following RHS (double support)
%
%
% Outputs:
%   neckChRatio     - [double] ratio of power at channels contralateral to the previous heel strike to one ipsilateral
%   standBL         - [matrix,chan x freq x pnts] sbaseline corrected time-frequency data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% May 2020; Last revision: 12-May-2020

%------------- BEGIN CODE --------------
% check inputs
narginchk(6,6)

if ~all([size(TFdata, 1), size(TFdata,2)] == size(TFbaseline))
    error('data and baseline matrix do not have the same number of rows and/or columns');
elseif any([~ismember(pntsLHS, 1:size(TFdata,3)), ~ismember(pntsRHS, 1:size(TFdata,3))])
    error('Some sample indices exeed the numer time points in the EEG data (dimension 3)');
elseif any([~ismember(neckChanL, 1:size(TFdata,1)), ~ismember(neckChanR, 1:size(TFdata,1))])
    error('Some channel indices exeed the numer time points in the EEG data (dimension 3)');
end

% calculation
standBL = TFdata./repmat(TFbaseline,1,1,size(TFdata,3)); % standing baseline correction
datRHS = squeeze(sum(sum(standBL(:,:,pntsRHS),2),3)); % sum over time and frequency -> chan
datLHS = squeeze(sum(sum(standBL(:,:,pntsLHS),2),3));

datIpsi = sum(datRHS(neckChanR))+sum(datLHS(neckChanL));
datContra = sum(datRHS(neckChanL))+sum(datLHS(neckChanR));

% calculate and save ratio
neckChRatio = 1-(datIpsi./datContra);

%------------- END OF CODE --------------
end
function SW = swRatio(TFdata, TFbaseline)
%% swRatio- ADD ONE-LINER HERE!
% more info
%
% Syntax:  SW = swRatio(TFdata, TFbaseline)
%
% Inputs:
%   TFdata          - [matrix,chan x freq x pnts] time-frequency transformed EEG (averaged and time-normalized to the gait cycle)
%   TFbaseline      - [matrix, chan x freq] time-frequency transformed EEG baseline (i.e standing),
%
%
% Outputs:
%   SW              - [double] standing/walking ratio
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% May 2020; Last revision: 12-May-2020

%------------- BEGIN CODE --------------

% check inputs
narginchk(2,2)

if ~all([size(TFdata, 1), size(TFdata,2)] == size(TFbaseline))
    error('data and baseline matrix do not have the same number of rows and/or columns');
end

% calculate
SWratio = repmat(TFbaseline,1,1,size(TFdata,3))./TFdata;    % ratio of standing/walking power
SW = mean(1-SWratio,'all');                                 % subtract every value from 1, then average

%------------- END OF CODE --------------
end
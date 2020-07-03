function ROIratio = sourceRatio(ROI, wholeCortex)
%% ROIratio- ADD ONE-LINER HERE!
% more info
%
% Syntax:  ROIratio = sourceRatio(ROI, wholeCortex)
%
% Inputs:
%   ROI         - [vector, vertices] source activations over the gait cycle at a motor ROI 
%   wholeCortex - [vector, vertices] source activations over the gait cycle of the whole cortex
%
%
% Outputs:
%   ROIratio    - [double] M1 source activation ratio
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% May 2020; Last revision: 12-May-2020

%------------- BEGIN CODE --------------

% check inputs
narginchk(2,2)

if ~ismatrix(ROI)||~ismatrix(wholeCortex)
    error('At least one of your inputs is not a matrix');
end

% calculate activity per vertex
% sum activity divide by numver of vertices (to normalize)
actROI          = sum(ROI)./size(ROI,1); 
actWholeCortex  = (sum(wholeCortex)-sum(ROI))./(size(wholeCortex,1)- size(ROI,1)); 

% subtract ratio from 1
ROIratio        = (actWholeCortex./actROI)';

%------------- END OF CODE --------------
end
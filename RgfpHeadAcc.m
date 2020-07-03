function R = RgfpHeadAcc(ERP, ACC, varargin)
%% RgfpHeadAcc - calculate variance of head acceleration that is explained by global field potential of gaitERP
% cross-correlates the global field potential of a gait ERP with the root mean square of the head acceleration
%
% Syntax:  R = RgfpHeadAcc(ERP, ACC, lag)
%
% Inputs:
%    ERP - time series of gait ERP, chans x pnts
%    ACC - time series of head acceleration, chans x pnts
% Note: both vectors have to have the same length
%
% Optional inputs:
%   lag - maximal number of samples cross-correlation may be shifted
%
% Outputs:
%    R - squared correlation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Nadine Jacobsen, nadine.jacobsen@uni-oldenburg.de
% May 2020; Last revision: 12-May-2020

%------------- BEGIN CODE --------------
% set default params


% check inputs
% at least 2, numeric vectors, third a number, set default lag to 5
narginchk(2,3)

 if ~ismatrix(ERP)||~ismatrix(ACC)
    error('EEG and/or acceleration are not a matrix');
 elseif size(ERP,2)~=size(ACC,2)
    error('EEG and accereloration signal do not have the same length');
end

% check optional inputs
if isempty(varargin)
    lag = 5;
elseif ~isempty(varargin) 
    if ~isinteger(varargin)
        error('Lag has to be an integer');
    else
        lag = varargin;
    end
end

%%
% calculate global field potential of gait ERP
gfp = std(ERP);
gfp = gfp-mean(gfp);

% calculate root mean square of gait ERP
rmsHead = rms(ACC);
rmsHead = rmsHead-mean(rmsHead);

% cross-correlate both time series
% subtract mean of each signal before to not bias cross-correlation
[all_r, ~] = xcorr(gfp, rmsHead, lag, 'coeff'); %use normalized cross-correlation

r = max(all_r); % use strongest correlation
R = r^2;        % square to obtein explained variance

%------------- END OF CODE --------------
end
function [eta] = fInterpEta(xq, yq, tq, waveData, method, extrap)

% eta = fInterpEta(xq, yq, tq, waveData, method, extrap)
% ------------------------------------------------------------------------
% Interpolates the free surface at query points (xq, yq, tq) according to
% a pre-computed data file.
%
% inputs:
%   xq, yq [m], t [s] - spatial and time points.
%                       Acceptable input formats for x and t:
%                         1) Only one is a non-scalar.
%                         2) Arrays of the same size.
%   waveData - path to the data file, OR the data file loaded as a struct.
%   method - interpolation method for scatteredInterpolant, see MATLAB
%            documentation. Defaults to 'linear', higher orders unverified.
%   extrap - extrapolaton method for scatteredInterpolant, see MATLAB
%            documentation. Defaults to 'none', where out-of-range queires
%            will return NaN. Extrapolation results are unverified.
% outputs:
%   eta [m] - free surface elevation relative to SWL.
% ------------------------------------------------------------------------
% lm808, 10/2019.
% github.com/lm808, all rights reserved.

if ischar(waveData)
    waveData=load(waveData);
end

switch nargin
    case 4
        method = 'linear';
        extrap = 'none';
    case 5
        extrap = 'none';
end

%% Check query points
uniTq = unique(tq);
if isempty(uniTq)
    eta = [];
    return
end

for i = 1:numel(uniTq)
    if single(min(abs(waveData.t-uniTq(i))))>=1e-10%zeroTol
        error('Time-step mismatch between query and pre-computed wave data.')
    end
end

nxq = size(xq);
nyq = size(yq);
ntq = size(tq);

if nxq==nyq & nxq==ntq
    eta = zeros(size(xq+yq+tq));
else
    eta = zeros(size(xq+yq+tq));
    s = find([isscalar(xq),isscalar(yq),isscalar(tq)]==0);
    if length(s)>1
        error('Query points dimension mismatch.')
    end
    switch s
        case 1
            yq = repmat(yq,nxq);
            tq = repmat(tq,nxq);
        case 2
            xq = repmat(xq,nyq);
            tq = repmat(tq,nyq);
        case 3
            xq = repmat(xq,ntq);
            yq = repmat(yq,ntq);
    end
end

%% Establish interpolant
F = scatteredInterpolant(reshape(waveData.X,waveData.nx*waveData.ny,1),reshape(waveData.Y,waveData.nx*waveData.ny,1),reshape(waveData.ETA(:,:,1),waveData.nx*waveData.ny,1));
F.Method = method;
F.ExtrapolationMethod = extrap;

%% Interpolate
for i = 1:length(uniTq)
    % progress report
    dgt = length([num2str(i-1),' /',num2str(length(uniTq))]);
    if i>1
        fprintf(1,repmat('\b',1,dgt));
    end
    fprintf(1,' %u/%u',i,length(uniTq))

    s = find(tq==uniTq(i));
    p = find(abs(waveData.t-uniTq(i))==min(abs(waveData.t-uniTq(i))));
    if length(p)>1
        error('More than one time step is found to match query.');
    end
    F.Values = reshape(waveData.ETA(:,:,p),waveData.nx*waveData.ny,1);
    eta(s) = F(xq(s),yq(s));
end
% progress report
dgt = length([num2str(i),' /',num2str(length(uniTq))]);
fprintf(1,repmat('\b',1,dgt));

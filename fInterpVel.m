function [u, v, w] = fInterpVel(xq, yq, zq, tq, waveData, method, extrap)

% [u,v,w] = fInterpVel(xq, yq, zq, tq, waveData, method, extrap)
% ------------------------------------------------------------------------
% Interpolates the fluid velocities at query points (xq, yq, zq, tq)
% according to a pre-computed data file.
%
% inputs:
%   xq, yq, zq [m], t [s] - spatial and time points.
%                           Acceptable input formats for x and t:
%                             1) Only one is a non-scalar.
%                             2) Arrays of the same size.
%   waveData - path to the data file, OR the data file loaded as a struct.
%   method - interpolation method for scatteredInterpolant, see MATLAB
%            documentation. Defaults to 'linear', higher orders unverified.
%   extrap - extrapolaton method for scatteredInterpolant, see MATLAB
%            documentation. Defaults to 'none', where out-of-range queires
%            will return NaN. Extrapolation results are unverified.
% outputs:
%   [u, v, w] - fluid velocty components in 3 Cartesian directions.
% ------------------------------------------------------------------------
% lm808, 10/2019.
% github.com/lm808, all rights reserved.

if ischar(waveData)
    waveData = load(waveData);
end

switch nargin
    case 5
        method = 'linear';
        extrap = 'none';
    case 6
        extrap = 'none';
end

%% Check query points
uniTq = unique(tq);
if isempty(uniTq)
    u = [];
    v = [];
    w = [];
    return
end

for i = 1:numel(uniTq)
    if min(abs(waveData.t-uniTq(i)))>=1e-10
        error('Time-step mismatch between query and pre-computed wave data.')
    end
end

nxq = size(xq);
nyq = size(yq);
nzq = size(zq);
ntq = size(tq);

if nxq==nyq & nxq==ntq & nxq == nzq
    u = zeros(size(xq+yq+zq+tq));
    v = zeros(size(xq+yq+zq+tq));
    w = zeros(size(xq+yq+zq+tq));
else
    u = zeros(size(xq+yq+zq+tq));
    v = zeros(size(xq+yq+zq+tq));
    w = zeros(size(xq+yq+zq+tq));

    s = find([isscalar(xq),isscalar(yq),isscalar(zq),isscalar(tq)]==0);

    if length(s)>1
        error('Query points dimension mismatch.')
    else
        switch s(1)
            case 1
                yq = repmat(yq,nxq);
                zq = repmat(zq,nxq);
                tq = repmat(tq,nxq);
            case 2
                xq = repmat(xq,nyq);
                zq = repmat(zq,nyq);
                tq = repmat(tq,nyq);
            case 3
                xq = repmat(xq,nzq);
                yq = repmat(yq,nzq);
                tq = repmat(tq,nzq);
            case 4
                xq = repmat(xq,ntq);
                yq = repmat(yq,ntq);
                zq = repmat(zq,ntq);
        end
    end
end

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

    F = scatteredInterpolant(waveData.x{p},waveData.y{p},waveData.z{p},waveData.ux{p});
    F.Method = method;
    F.ExtrapolationMethod = extrap;
    u(s) = F(xq(s),yq(s),zq(s));

    F.Values = waveData.uy{p};
    v(s) = F(xq(s),yq(s),zq(s));

    F.Values = waveData.uz{p};
    w(s) = F(xq(s),yq(s),zq(s));
end

% progress report
dgt = length([num2str(i),' /',num2str(length(uniTq))]);
fprintf(1,repmat('\b',1,dgt));


%% Set all points above water level to zero
eta = fInterpEta(xq,yq,tq,waveData,method,extrap);
u(zq>eta) = 0;
v(zq>eta) = 0;
w(zq>eta) = 0;

% u(single(zq)>single(eta)) = 0;
% v(single(zq)>single(eta)) = 0;
% w(single(zq)>single(eta)) = 0;

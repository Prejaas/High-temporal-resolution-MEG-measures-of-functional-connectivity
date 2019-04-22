function [dyn_wc, mean_wc] = dyn_wavcoh(data, options)
% Computes the wavelet coherence (WC) without pairwise
% leakage correction/pairwise orthogonalisation
% Tewarie et al., 2019 Tracking dynamic brain networks using high temporal 
%                      resolution MEG measures of functional connectivity 
%
%       Input:  - data or timeseries (N x M), where N is channels or nodes and
%               M number of samples
%               - required 'options' are
%                 options.sample        sampling frequency
%                 options.freq_low      lower bound of frequency band
%                 options.freq_high     upper bound of frequency band
%
%       Output: 
%               - dyn_wc weighted connectivity tensor spanning the
%               experiment (N x N X M)
%               - mean_wc weighted connectivity matrix (N x N), average/collapse
%               over N
%
%
% code is based on and parameters follow from:
% 
% Grinsted A, Moore JC, Jevrejeva S. 2004. Application of the cross wavelet transform 
% and wavelet coherence to geophysical time series. Nonlinear Process Geophys. 11:561–566.
%
% Code requires the wavelet coherence toolbox which can be downloaded from 
% https://sites.google.com/a/glaciology.net/grinsted/wavelet-coherence
% or https://github.com/grinsted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('.../grinsted-wadatalet-coherence-d987ea4'))

% frequency to scale
Freq = scal2frq(1:options.sample,'morl',1/options.sample );
[~, index_low] = min(abs(Freq - options.freq_high));
[~, index_high] = min(abs(Freq - options.freq_low));
no_chan = size(data,1);
M = size(data,2);

% parameters
Pad = 0;
S0 = index_low;
MaxScale = index_high;
J1 = round((index_high - index_low)/3); 
Mother = 'Morlet';

dyn_wc = zeros(no_chan,no_chan,size(data,2));
for r1 = 1:no_chan
    x = data(r1,:);

    for r2 = 1:no_chan
        if r1 ~= r2
        y = data(r2,:);      
        wc = wtc(x,y,'Pad',Pad,'S0',S0,'MaxScale',MaxScale,'Mother',Mother,'J1',J1);
        dyn_wc(r1,r2,:) = mean(wc,1); % average across frequencies within the frequency band 
        end
    end
end

% delete last and first part of the data
dyn_wc(:,:,M-25:M-1) = [];
dyn_wc(:,:,1:25) = [];
mean_wc = mean(dyn_wc,3); % mean connectivity

fprintf('computed wavelet coherence (no leakage correction) \n')
    
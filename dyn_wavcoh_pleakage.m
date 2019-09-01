function [dyn_wc, mean_wc] = dyn_wavcoh_pleakage(data, options)
% Computes the wavelet coherence (WC) with pairwise
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
scales = frq2scale(linspace(options.freq_low,options.freq_high,10),'morl',options.sample);
no_chan = size(data,1);
M = size(data,2);

% parameters
Pad = 0;
S0 = min(scales);
MaxScale = max(scales);
J1 = round((MaxScale - S0)/3); 
Mother = 'Morlet';

dyn_wc = zeros(no_chan,no_chan,size(data,2));
for r1 = 1:no_chan
    x = data(r1,:);
    inv_x = pinv(x);
    for r2 = 1:no_chan
        if r1 ~= r2
        y = data(r2,:);
 
        % pairwise leakage correction/orthogonalisation
        y_cor = y - x*(inv_x'*y');    
        
        wc = wtc(x,y_cor,'Pad',Pad,'S0',S0,'MaxScale',MaxScale,'Mother',Mother,'J1',J1);
        dyn_wc(r1,r2,:) = mean(wc,1); % average across frequencies within the frequency band 
        end
    end
end

% delete last and first part of the data
dyn_wc(:,:,M-25:M-1) = [];
dyn_wc(:,:,1:25) = [];
mean_wc = mean(dyn_wc,3); % mean connectivity

fprintf('computed wavelet coherence (with leakage correction) \n')
end

function scales = frq2scale(freq, wname, fs)   
% freq : frequencies to convert to scales
% wname : name of the wavelet : for example 'morl'
% fs : frequency of the sampling in Hertz

    if isempty(freq) , f = freq; return; end
    if nargin == 2, fs = 1; end
    err = (min(size(freq))>1) | (min(freq)<eps);
    if err
        error(message('Wavelet:FunctionArgVal:Invalid_ScaVal'))
    end
    if fs <= 0
        error(message('Invalid Value for fs !'))
    end
    scales = centfrq(wname)./(freq*fs);
end

    
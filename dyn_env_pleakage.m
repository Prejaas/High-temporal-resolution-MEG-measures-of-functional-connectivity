function [dyn_IAC, mean_IAC] = dyn_env_pleakage(data)
% Computes the Instantaneous Amplitude Correlation (IAC) with pairwise
% leakage correction/pairwise orthogonalisation
% Tewarie et al., 2019 Tracking dynamic brain networks using high temporal 
%                      resolution MEG measures of functional connectivity 
%
%       Input:  - data or timeseries (N x M), where N is channels or nodes and
%               M number of samples
%
%       Output: 
%               - dyn_IAC weighted connectivity tensor spanning the
%               experiment (N x N X M)
%               - mean_IAC weighted connectivity matrix (N x N), average/collapse
%               over N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_chan = size(data,1);
M = size(data,2);
data = zscore(data,0,2);

% loop over channels
dyn_IAC = zeros(no_chan,no_chan,size(data,2));
for r1 = 1:no_chan
    x = data(r1,:);
    inv_x = pinv(x);
    env_x = abs(hilbert(x));
            
    % loop over channels
    for r2 = 1:no_chan
        if r1 ~= r2
        y = data(r2,:);
        
        % pairwise leakage correction/orthogonalisation
        y_cor = y - x*(inv_x'*y');
       
        env_y = abs(hilbert(y_cor));
        dyn_IAC(r1,r2,:) = env_x.*env_y;
        
        end
    end
end

% delete last and first part of the data
dyn_IAC(:,:,M-25:M-1) = [];
dyn_IAC(:,:,1:25) = [];
mean_IAC = mean(dyn_IAC,3);

fprintf('computed instantaneous amplitude correlation (with leakage correction) \n')


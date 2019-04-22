function [pdd, pddnode, meanpdd] = dyn_pdd(data,sample,mindur,meanfreq)
% Computes the phase difference derivative (PDD) without pairwise
% leakage correction/pairwise orthogonalisation
%
% Tewarie et al., 2019 Tracking dynamic brain networks using high temporal 
%                      resolution MEG measures of functional connectivity 
%
% input:     - timeseries: N x M, with N no of channels and M no of samples
%            - mindur: minimum period of phase coupling in samples that is considered to be genuine coupling
%             (usually at least one cycle of an oscillation, duration of a cycle in samples)
%            - meanfreq: mean frequency of the frequency band 
%
% output:
%           - dphasematrix: N x N x M matrix of d_diff_phase/dt for all
%             pair combinations 
%           - dphasenode: average d_diff_phase/dt for each channel (N x M matrix)
%           - phasematrix: d_diff_phase/dt averaged over time (static connectivity N x N matrix)
%             N = number of nodes/channels
%             M = duration of recording in samples
%      
%
% original paper on PDD: Breakspear M, Williams LM, Stam CJ. 2004. A novel method for the topographic 
% analysis of neural activity reveals formation and dissolution of “dynamic cell assemblies.” J Comput Neurosci. 16:49–68.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(data,1);
M = size(data,2);
pdd = zeros(N,N,M-1);

% loop over channels
for j=1:N
    x = data(j,:);
    ht1 = hilbert(x);
    ht1 = bsxfun(@minus,ht1,sum(ht1,1)/M);
    phase1 = angle(ht1);   
    uphase1 = unwrap(phase1);

    % loop over channels
    for k=1:N
        if k~=j
            y = data(k,:);
            ht2 = hilbert(y);
            ht2 = bsxfun(@minus,ht2,sum(ht2,1)/M);
            phase2 = angle(ht2);
            uphase2 = unwrap(phase2);
            phasedif = uphase1 - uphase2;
            
            % compute pdd
            diftemp = exp(-(abs(diff(phasedif).*sample./(meanfreq))));
            diftemp_cut = filter_coupling(diftemp,mindur);
            pdd(j,k,:) = smooth(diftemp_cut,5,'moving');
        end
    end
end

pdd(:,:,M-25:M-1) = [];
pdd(:,:,1:25) = [];
pddnode = squeeze(mean(pdd,1));
meanpdd = squeeze(mean(pdd,3));

fprintf('computed phase difference derivative (no pairwise leakage correction) \n')
end
            
function k = filter_coupling(input,sync_cut)
% 'filter' phase derivative data take out random fluctuations. the PDD can
% be near zero if there is a switch between phase lag/lead. These near zero PDD values are not a sign of coupling and 
% can be removed by only taking into account episodes of d(diffpahse)/dt near zero with a mininum duration of one oscillation: mindur

k = zeros(size(input,1),size(input,2));
for no_chan = 1: size(input,1)
    c = input(no_chan,:);
m=1;
while m<numel(c)
    count = 0;
    if c(m)>0.4
            a=m; %start
            count = count + 1;
            m = m + 1;
        
        while c(m)>0.4 && m<numel(c)
            count = count+1;
            m = m + 1;
        end
        if count < sync_cut % 
             k(no_chan,a:a+count)=0;
        else k(no_chan,a:a+count)=c(a:a+count);
        end
    else
            k(no_chan,m) = 0;
            m=m+1;
    end
end
  
if numel(k) == numel(c)-1
    k(no_chan,m)=0;
end
end
end

function peak_loc = get_peakloc(data,win_width,sample,meanfreq,no_osc)
% Extracting data driven windows from recurrence plots 
%
% Tewarie et al., 2019 Tracking dynamic brain networks using high temporal 
%                      resolution MEG measures of functional connectivity 
%
% 
% Output:       - peak_loc refer to local maxima in the gradient matrix of
%                 the recurrence plot, and thus correspond to transitions
%                 between states in state space. output is the index of
%                 these local maxima in the data
%
%               - data refers to timeseries, can be frequency filtered (N x
%               M), where N is no channels and M number of samples
%               - sample is sampling frequency
%               - meanfreq refers to the mean frequency in the frequency
%               band
%               - win_width refers to a working memory limitation rather
%               than limitation of the method. For 32GB system, a
%               correlation matrix of M x M (sample x sample) of all the data is not feasible.
%               Therefore the recurrence plot/correlation matrix is
%               computed for subsequent window of size win_width (in
%               samples). 
%               - no_osc refers to the size of the blocks around the
%               diagonals, in which we are interested. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Hilbert transform/amplitude envelope
M = size(data,2);
env = abs(hilbert(data'));
env(M-25:M-1,:) = [];
env(1:25,:) = [];

% initialisation
no_win = floor(size(data,2)/win_width);
beg = 1;
eind = win_width;
peak_loc=[];

% create toeplitz like matrix
N = win_width;
K = N*no_osc;
CIJ = zeros(N);
CIJ1 = ones(N);
KK = 0;
cnt = 0;
seq = 1:N-1;
while (KK<K)
    cnt = cnt + 1;
    dCIJ = triu(CIJ1,seq(cnt))-triu(CIJ1,seq(cnt)+1);
    dCIJ = dCIJ+dCIJ';
    CIJ = CIJ + dCIJ;
    KK = sum(sum(CIJ));
end

% loop over windows
for win = 1:no_win
    if eind < M
        
        % compute recurrence plot
        recur_temp = corr(env(beg:eind,:)');
        vert_lns = sum(abs(gradient(recur_temp).*CIJ));
        vert_lns_sm = smooth(vert_lns,round(sample/meanfreq(1)*0.2));       
        [~, loc] = findpeaks(vert_lns_sm,'MinPeakDistance',round((sample/meanfreq)));
        peak_loc = cat(1,peak_loc,loc+((win-1).*win_width));
        beg = beg + win_width;
        eind = eind + win_width;
    end
    if eind > M
        win_width = size(data,2) - beg;
        % create toeplitz like matrix
        N = win_width;
        K = N*no_osc;
        CIJ = zeros(N);
        CIJ1 = ones(N);
        KK = 0;
        cnt = 0;
        seq = 1:N-1;
        while (KK<K)
            cnt = cnt + 1;
            dCIJ = triu(CIJ1,seq(cnt))-triu(CIJ1,seq(cnt)+1);
            dCIJ = dCIJ+dCIJ';
            CIJ = CIJ + dCIJ;
            KK = sum(sum(CIJ));
        end
        
        % compute recurrence plot
        recur_temp = corr(env(beg:end,:)');
        vert_lns = sum(abs(gradient(recur_temp)));
        vert_lns_sm = smooth(vert_lns,round(sample/meanfreq(1)*0.2));
        [~, loc] = findpeaks(vert_lns_sm,'MinPeakDistance',round((sample/meanfreq)));
        peak_loc = cat(1,peak_loc,beg + loc);
    end
    
end
fprintf('extracted borders for data driven windows based on recurrence plots \n')

end



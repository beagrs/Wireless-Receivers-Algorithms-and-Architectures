clear all
close all
clc

rng(12)

% Set parameters
os_factor = 4;
SNR  = 0;

load ber_pn_seq
load pn_sequence

data_length = length(ber_pn_seq)/2; % Number of QPSK data symbols

% convert SNR from dB to linear
SNRlin = 10^(SNR/10);

% add awgn channel
rx_signal = signal + sqrt(1/(2*SNRlin)) * (randn(size(signal)) + 1i*randn(size(signal)) );

% Matched filter
filtered_rx_signal = matched_filter(rx_signal, os_factor, 6);

% Frame synchronization
data_idx = frame_sync(filtered_rx_signal, os_factor); % Index of the first data symbol

data = zeros(1,data_length);
data2 = zeros(1,data_length);

cum_err = 0;
diff_err = zeros(1,data_length);
epsilon  = zeros(1,data_length);

frame_sync_length = 100;

for ii=1:data_length
    %indice va a passo oversample
     idx_start  = data_idx+(ii-1)*os_factor;
     %prendo il simbolo + sue copie oversampled per fare estimation of
     %epsilon
     idx_range  = idx_start:idx_start+os_factor-1;
     segment    = filtered_rx_signal(idx_range);
    
     % Estimate timing error epsilon
     [epsilon_tmp, cum_err] = epsilon_estimation_solution(segment, cum_err);
     epsilon(ii) = epsilon_tmp;
     
     % offset of sample
     sample_diff   = floor(epsilon(ii)*os_factor); % integer offset
     int_diff      = mod(epsilon(ii)*os_factor,1); % fractional offset (interval [0 1) )
    
     % linear
     y = filtered_rx_signal(idx_start+sample_diff:idx_start+sample_diff+1);
     y_hat = ((1-int_diff)*y(1)+int_diff*y(2));
     data(ii) = y_hat;
     
     % cubic
     y2 = filtered_rx_signal(idx_start+sample_diff-1:idx_start+sample_diff+2); %LOOK AT THE DEF
     c = (1/6) *[-1, 3,-3, 1; 3, -6, 3, 0; -2 ,-3 ,6, -1; 0, 6, 0 ,0] * y2;
     y_hat2 = c(1) *int_diff.^3 +c(2) *int_diff.^2 +c(3)*int_diff+c(4);
     data2(ii) = y_hat2;  

end

BER_lin = mean(demapper(data) ~= ber_pn_seq)
BER_cub = mean(demapper(data2) ~= ber_pn_seq)


% Plot epsilon
figure
plot(1:data_length, epsilon)
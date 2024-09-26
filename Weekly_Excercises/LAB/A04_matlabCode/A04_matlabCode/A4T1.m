clear all 
close all
clc

rng(32);
% Set parameters
SNR = 10;
os_factor = 4;

% Load image file
load task4.mat

data_length = prod(image_size) * 8 / 2; % Number of QPSK data symbols

% convert SNR from dB to linear
SNRlin = 10^(SNR/10);

% add awgn channel
rx_signal = signal + sqrt(1/(2*SNRlin)) * (randn(size(signal)) + 1i*randn(size(signal)) ); 

% Matched filter
filtered_rx_signal = matched_filter(rx_signal, os_factor, 6); % 6 is a good value for the one-sided RRC length (i.e. the filter has 13 taps in total)

% Frame synchronization
data_idx = frame_sync(filtered_rx_signal, os_factor); % Index of the first data symbol

data = zeros(1,data_length);
data2 = zeros(1,data_length);

cum_err = 0;
diff_err = zeros(1,data_length);
epsilon  = zeros(1,data_length);
vector_transform = [1, -1j, -1, 1j];
for i=1:data_length
    %sostanzialmente prendo i dati a passo L downsample
     idx_start  = data_idx+(i-1)*os_factor;
     %in idx range mi prendo il simbolo +copie oversamples
     idx_range  = idx_start:idx_start+os_factor-1;
     segment    = filtered_rx_signal(idx_range);
    
     % Estimate timing error epsilon
     % TODO
     %seguire esattamente la formula
     
     diff_err(i) = vector_transform*abs(segment).^2;
     cum_err     = cum_err + diff_err(i);
     epsilon(i)  = -1/(2*pi) *angle(cum_err);
     
     y_hat = linear_interpolation_solution(epsilon(i), filtered_rx_signal, os_factor, idx_start);
     data(i) = y_hat;
     
end

image = image_decoder(demapper(data), image_size);

figure;
plot(1:data_length, epsilon)
imshow(image)

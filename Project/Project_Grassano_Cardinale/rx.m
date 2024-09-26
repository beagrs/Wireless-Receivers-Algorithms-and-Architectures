function [rxbits conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

%% RF Downconversion

% Define the time frame
time = 0:(1/conf.f_s) : size(rxsignal,1)/conf.f_s - (1/conf.f_s);

% Perform the downconversion by removing the RF carrier
rxsignal_down = (rxsignal.*exp(-1i*2*pi*conf.f_c.*time.'));

%% Filter the Downconverted RX Signal

% Apply the low pass filter given
f_cutoff = conf.BW_BB + 0.05 * conf.BW_BB;      % Define the filter cutoff as 5% above the baseband BW
rx_signal_filtered = ofdmlowpass(rxsignal_down, conf, f_cutoff);

%% Frame Synchronization

OFDM_init = frame_sync(rx_signal_filtered, conf);    % First sample of the OFDM signal

% Extract only the OFDM part of the signal

% Each symbol is preceded by CP, then the total number of samples in the
% OFDM frame is the num of samples of an OFDM symbol + num of samples in
% the CP times the num of OFDM symbols in the frame
OFDM_signal_len = (conf.nb_tot_OFDM_symbs) * (conf.os_ifft_OFDM_symb_len + conf.N_samples_CP);  
OFDM_rx_signal = rx_signal_filtered(OFDM_init: OFDM_init + OFDM_signal_len - 1);

%% Remove The CP

% Define the length of the frame without CP
len_OFDM_no_CP = (conf.nb_tot_OFDM_symbs) * conf.os_ifft_OFDM_symb_len;

% Define a vector storing OFDM symbols without CP
OFDM_symbols_vec = zeros(len_OFDM_no_CP, 1);

for i = 1: (conf.nb_tot_OFDM_symbs)
    idx = (i-1) * conf.os_ifft_OFDM_symb_len + 1;       % index where an OFDM symbol starts in the resulting vector
    idx2 = (i-1) * (conf.os_ifft_OFDM_symb_len + conf.N_samples_CP) + conf.N_samples_CP +1;     % Index where an OFDM symbols starts in the initial vector (symbols + CP)
    OFDM_symbols_vec(idx: idx + conf.os_ifft_OFDM_symb_len - 1) = OFDM_rx_signal(idx2 : idx2 + conf.os_ifft_OFDM_symb_len -1);
end

%% Serial To Parallel Coversion

% We convert the OFDM vector signal into a matrix.
% Each column of the matrix will store an OFDM symbol, i.e. BPSK/QPSK symbols
% sent at the same time instant but at different carrier frequencies.
% Each row stores BPSK/QPSK symbols sent at different time instants but at the
% same carrier frequency.

OFDM_rx_matrix = reshape(OFDM_symbols_vec, [conf.os_ifft_OFDM_symb_len, (conf.nb_tot_OFDM_symbs)]);

%% FFT + Downsampling

% Make the OFDM symbols pass through the OSFFT block

% Defiine a new matrix with the same properties of the one before
OFDM_rx_matrix_down = zeros(conf.nb_subcarriers, (conf.nb_tot_OFDM_symbs));

for i = 1:(conf.nb_tot_OFDM_symbs)
    OFDM_rx_matrix_down(:,i) = osfft(OFDM_rx_matrix(:,i), conf.os_factor).';
end

conf.rx_frame_fft = OFDM_rx_matrix_down;

%% Channel Equalization

% (This task is performed after the P/S conversion as in the block diagram seen in class. 
% For coding purposes, here it is done before. It is simply because it seemed easier to work with a parallel
% signal where each OFDM symbol is defined in a column)

% Each column stores an OFDM symbol after the channel correction
OFDM_rx_matrix_down = channel_equalize_f(OFDM_rx_matrix_down, conf);


%% Parallel to Series

% Remove the training
OFDM_rx_matrix_down_no_training = remove_training(OFDM_rx_matrix_down, conf);

% Reshape into a vector
if strcmp(conf.training_type,'Comb') 
    
    % Find the number of data symbs in the frames before the last
    nb_symb_before = (size(OFDM_rx_matrix_down_no_training, 2) - 1) * size(OFDM_rx_matrix_down_no_training, 1);

    % Reshape into a vector the first columns minus the last
    OFDM_rx_vector_down = reshape(OFDM_rx_matrix_down_no_training(:, 1 : end-1), [nb_symb_before, 1]);

    % Find number of data symb in the last column
    nb_symb_last_col = (conf.OFDM_symbs_per_frame * conf.nb_subcarriers) - nb_symb_before;

    % Take the data in the last frame and append it to the vector of the
    % serial frame
    OFDM_rx_vector_down = [OFDM_rx_vector_down ; OFDM_rx_matrix_down_no_training(1 : nb_symb_last_col, end)];

else
    vector_len = (conf.OFDM_symbs_per_frame) * conf.nb_subcarriers;
    OFDM_rx_vector_down = reshape(OFDM_rx_matrix_down_no_training, [vector_len, 1]);
end


%Normalize the signal
avg_E = average_energy_f(OFDM_rx_vector_down);
OFDM_rx_vector_down = normalized_constellation_f(OFDM_rx_vector_down, avg_E);

%% Plot the constellation of the rx symbols
figure
plot(real(OFDM_rx_vector_down), imag(OFDM_rx_vector_down), '.');
hold on;
yline(0, 'r--', 'LineWidth', 2);
xline(0, 'r--', 'LineWidth', 2);
grid on;
axis square; box on;
%% Demapper

rxbits = demapper(OFDM_rx_vector_down);
end










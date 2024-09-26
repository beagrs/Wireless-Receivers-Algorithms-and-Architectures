function [txsignal, conf] = tx(txbits,conf, k)
% Digital Transmitter
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete transmitter
%   consisting of:
%       - modulator
%       - pulse shaping filter
%       - up converter
%   in digital domain.
%
%   txbits  : Information bits
%   conf    : Universal configuration structure
%   k       : Frame index
%

%% Generate Preamble

preamble = 1 - 2*lfsr_framesync(conf.npreamble);                    % Generate bits of preamble and map them in BPSK Gray coding
preamble_up = upsamplig_f(conf.preamble_os_factor, preamble).';     % Upsample the preamble
pulse = rrc(conf.preamble_os_factor, conf.alfa, conf.tx_filter_length * conf.preamble_os_factor);
preamble_pulse = conv(preamble_up, pulse, 'same');                  % Convolve the preamble with the pulse shaping filter
conf.preamble = preamble;
%% Training Symbol

% Define a known OFDM training symbol
% The OFDM symbol will have nb_subcarriers BPSK symbols
%Generate bits
nb_training_bits = conf.nb_subcarriers;
training_bits = randi([0, 1], nb_training_bits, 1);

% Map bits to BPSK Symbols with Gray Coding
training_BPSK_symbols = 1 - 2*training_bits;

% Assign the vector of BPSK symbols to the OFDM symbol saved in conf
conf.OFDM_training_symbol = training_BPSK_symbols;


%% Generate QPSK data symbols

% Convert txbits into a matrix where each row has 2 bits
txbits = reshape(txbits, [2, length(txbits)/2]).';

% Map them in QPSK with Gray coding
QPSK_data_symb = QPSK_grayMap(txbits);
 
%% Serial to Parallel Conversion

% Insert the training sequences
[OFDM_frame_vec, training_comb] = insert_training(QPSK_data_symb, conf);

% Save the comb training sequence in conf
conf.training_comb = training_comb;

%Generate the OFDM symbols of the frame (S to P conversion)

%Reshape the QPSK data symbols vector into an OFDM frame represented by a
%matrix.
%Each column of the matrix will have QPSK symbols sent in parallel in the same
%time frame T but in differet carrier frequencies.
%Each row will have QPSK symbols sent at the same carrier frequency but
%in different time frames.
nb_tot_OFDM_symbs = ceil(length(OFDM_frame_vec) / conf.nb_subcarriers);
conf.nb_tot_OFDM_symbs = nb_tot_OFDM_symbs;
padding_len = (nb_tot_OFDM_symbs * conf.nb_subcarriers) - length(OFDM_frame_vec);
OFDM_frame_vec = [OFDM_frame_vec ; zeros(padding_len, 1)];
OFDM_frame = reshape(OFDM_frame_vec, [conf.nb_subcarriers, nb_tot_OFDM_symbs]);
conf.OFDM_frame = OFDM_frame;

%% Inverse Fourier Transform

%We perform the IFFT + OS of each column. The result will be the pulse-shaped
%OFDM symbol
for i = 1: nb_tot_OFDM_symbs
    OS_IFFT_OFDM_frame(:,i) = osifft(OFDM_frame(:,i), conf.os_factor).';
end

%% CP Insertion

% Create a matrix storing in each column the samples of the CP.
% There is one CP per OFDM symbol, thus the choice of the number of columns.
CP_OFDM_samples = zeros(conf.N_samples_CP, nb_tot_OFDM_symbs);

% For every OFDM symbol take the last N_samples_CP samples and add them to
% the CP colummn.
for i = 1:nb_tot_OFDM_symbs
    CP_OFDM_samples(:,i) = OS_IFFT_OFDM_frame(end-conf.N_samples_CP+1:end,i);
end

%% Parallel - Series Conversion

% Vector that stores (serially) the OFDM symbols

os_ifft_OFDM_symb_len = length(OS_IFFT_OFDM_frame(:,1));    %length of the OFDM symbol after the os ifft block
os_ifft_CP_len = length(CP_OFDM_samples(:,1));              % length of the CP

conf.os_ifft_OFDM_symb_len = os_ifft_OFDM_symb_len;

% Define the length of the vector. It is the number of samples after the OS
% IFFT + the number of samples of the CP preceding every OFDM symbol.
serial_OFDM_signal_len = nb_tot_OFDM_symbs*os_ifft_OFDM_symb_len + nb_tot_OFDM_symbs*os_ifft_CP_len;
OFDM_signal = zeros(serial_OFDM_signal_len, 1);

for i = 1:nb_tot_OFDM_symbs
    idx = (i-1)*(os_ifft_OFDM_symb_len + os_ifft_CP_len) + 1;   % Starting idx of the CP of each OFDM symbol
    OFDM_signal(idx : idx + os_ifft_CP_len-1, 1) = CP_OFDM_samples(:,i);            % Add the CP preceding the i-th OFDM symbol
    idx = idx + os_ifft_CP_len;                                 %Starting idx of each OFDM symbol
    OFDM_signal(idx : idx + os_ifft_OFDM_symb_len-1, 1) = OS_IFFT_OFDM_frame(:,i);  % Add the OFDM symbol
end

%% Normalization

% Normalize OFDM signal
avg_E = average_energy_f(OFDM_signal);
OFDM_signal = normalized_constellation_f(OFDM_signal, avg_E);

% Normalize Preamble
avg_E_preamble = average_energy_f(preamble_pulse);
preamble_pulse = normalized_constellation_f(preamble_pulse, avg_E_preamble);

%% Building the final signal

% Initialize txsignal
txsignal = zeros(length(preamble_pulse) + serial_OFDM_signal_len,1);

%Assign the Preamble
txsignal(1:length(preamble_pulse)) = preamble_pulse;

%Assign the OFDM serial signal
txsignal(length(preamble_pulse) + 1 : end) = OFDM_signal;

%% RF Upconversion

time = 0:(1/conf.f_s) : size(txsignal,1)/conf.f_s - (1/conf.f_s);
txsignal = real(txsignal.*exp(1i*2*pi*conf.f_c.*time.'));

end

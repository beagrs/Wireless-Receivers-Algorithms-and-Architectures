clear all
close all
clc

% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal

% Configuration Values
conf.audiosystem = 'bypass'; % Values: 'matlab','native','bypass'

conf.f_s     = 48000;           % sampling rate  
conf.f_sym   = 100;             % symbol rate
conf.nframes = 1;               % number of frames to transmit
conf.nbits   = 2000;            % number of bits 
conf.modulation_order = 2;      % BPSK:1, QPSK:2
conf.f_c     = 8000;            % Carrier frequency [Hz]
conf.alfa    = 0.11;            % Roll-off factor of the rrc filter
conf.tx_filter_length = 20;     % Length of the rrc filter
conf.preamble_os_factor = 96;   % Preamble OS factor
conf.npreamble  = 100;          % Preamble length
conf.preamble = 0;
conf.bitsps     = 16;           % bits per audio sample
conf.offset     = 0;

% Variables to be defined by the user

conf.channel_type = 'Comb';              % Options: 'Block', 'Block_Viterbi', 'Comb'
conf.training_type = 'Comb';             % Options: 'Block', 'Comb'
conf.nb_subcarriers = 1024;               % Number of subcarriers
conf.OFDM_symbs_per_frame = 20;           % Number of OFDM data symbols
conf.nb_symbs_between_training = 2;      % Number of data symbols between each training symbol in BLOCK
conf.comb_insertion_rate = 4;             % number of QPSK data symbols between each training symbol in COMB

% Predefined variables that will be updated inside the functions

conf.nb_tot_OFDM_symbs = 0;               % Total num of OFDM symbs in a frame (training + data)
conf.OFDM_training_symbol = 0;            % OFDM training symbol. Initialize to 0
conf.OFDM_frame = 0;                      % OFDM TX frame (data + training)
conf.training_comb = 0;                   % Training symbols used in Comb

% Other useful values

conf.nb_QPSK_symbs_per_frame = conf.nb_subcarriers * conf.OFDM_symbs_per_frame;
conf.nb_bits_per_OFDM_symb = conf.nb_subcarriers * conf.modulation_order;
conf.nb_bits_per_frame = conf.nb_bits_per_OFDM_symb * conf.OFDM_symbs_per_frame; % Num of data bits in an OFDM frame
conf.os_ifft_OFDM_symb_len = 0;                                 % Length of the OFDM symbol after OS
conf.nb_training_symbs = ceil(conf.OFDM_symbs_per_frame / conf.nb_symbs_between_training);
conf.f_spacing = 5;                                             % 1/T [Hz]. T = duration of OFDM symbol
conf.BW_BB = ceil((conf.nb_subcarriers +1)/2)*conf.f_spacing;   % BaseBand bandwidth
conf.N_samples_CP = ceil(conf.f_s / (2*conf.f_spacing));        % Number of samples of the CP
conf.os_factor = ceil(conf.f_s / (conf.f_spacing * conf.nb_subcarriers));   % OS factor of our system. It will feed OSIFFT and OSFFT.

if mod(conf.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
end

conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

% Init Section
% all calculations that you only have to do once



% Initialize result structure with zero
res.biterrors   = zeros(conf.nframes,1);
res.rxnbits     = zeros(conf.nframes,1);

% TODO: To speed up your simulation pregenerate data you can reuse
% beforehand.


% Results


for k=1:conf.nframes
    %% Generate the Transmitted Signal
   
    % Generate the bits
    txbits = randi([0 1],(conf.nb_bits_per_frame),1);
    
    % TODO: Implement tx() Transmit Function
    [txsignal, conf] = tx(txbits,conf,k);

    %% PLOT Spectrum of Transmitted Signal
    TXSIGNAL = fftshift(fft(txsignal));
    f_range = - conf.f_s/2 : conf.f_s/length(txsignal) : conf.f_s/2 - conf.f_s/length(txsignal);
    figure(1)
    plot(f_range(1:end), (abs(TXSIGNAL)), 'b');
    ylabel("|X_{TX}|")
    xlabel("f, Hz")
    title("Spectrum of Transmitted Signal")
    %axis square; box on;
    grid on;

    %%
    
    % % % % % % % % % % % %
    % Begin
    % Audio Transmission
    %
    
    % normalize values
    peakvalue       = max(abs(txsignal));
    normtxsignal    = txsignal / (peakvalue + 0.3);
    
    % create vector for transmission
    rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
    rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
    txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
    
%     wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
    audiowrite('out.wav',rawtxsignal,conf.f_s)  
    
    % Platform native audio mode 
    if strcmp(conf.audiosystem,'native')
        
        % Windows WAV mode 
        if ispc()
            disp('Windows WAV');
            wavplay(rawtxsignal,conf.f_s,'async');
            disp('Recording in Progress');
            rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
            disp('Recording complete')
            rxsignal = rawrxsignal(1:end,1);

        % ALSA WAV mode 
        elseif isunix()
            disp('Linux ALSA');
            cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
            system(cmd); 
            disp('Recording in Progress');
            system('aplay  out.wav')
            pause(2);
            disp('Recording complete')
            rawrxsignal = audioread('in.wav');
            rxsignal    = rawrxsignal(1:end,1);
        end
        
    % MATLAB audio mode
    elseif strcmp(conf.audiosystem,'matlab')
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
        recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
        record(recobj);
        disp('Recording in Progress');
        playblocking(playobj)
        pause(0.5);
        stop(recobj);
        disp('Recording complete')
        rawrxsignal  = getaudiodata(recobj,'int16');
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
        
    elseif strcmp(conf.audiosystem,'bypass')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal; 
    end
    
    % Plot received signal for debugging
    figure;
    plot(rxsignal);
    title('Received Signal')
    
    %
    % End
    % Audio Transmission   
    % % % % % % % % % % % %
    
    % TODO: Implement rx() Receive Function
    [rxbits conf]       = rx(rxsignal,conf);
    
    res.rxnbits(k)      = length(rxbits);  
    res.biterrors(k)    = sum(rxbits ~= txbits);
    
    
end

per = sum(res.biterrors > 0)/conf.nframes
ber = sum(res.biterrors)/sum(res.rxnbits)

%% Plot Channel vs Time
if (strcmp(conf.channel_type,'Block') || strcmp(conf.channel_type,'Block_Viterbi'))
    H_hat = conf.rx_frame_fft ./ conf.OFDM_frame;
    Ts = 1/conf.f_spacing;  % Time to send a sample of an OFDM symbol; 
    time = 0 : Ts : (conf.OFDM_symbs_per_frame + conf.nb_training_symbs - 1) * Ts;
    figure
    subplot(2,1,1); 
    for i = 1 : 32 : conf.nb_subcarriers
    plot(time, 20*log10(abs(H_hat(i, :))./ max(abs(H_hat(i, :))) ));
    hold on;
    f = conf.f_spacing * i;
    end
    xlabel('Time, s');
    ylabel('Magnitude, dB')
    title('Channel Magnitude over Time');
    grid on;
    subplot(2,1,2);
    for i = 1 : 32 : conf.nb_subcarriers
    plot(time, unwrap(angle(H_hat(i, :))));
    hold on;
    f = conf.f_spacing * i;
    end
    xlabel('Time, s');
    ylabel('Phase, rad')
    title('Channel Phase over Time');
    grid on;
end


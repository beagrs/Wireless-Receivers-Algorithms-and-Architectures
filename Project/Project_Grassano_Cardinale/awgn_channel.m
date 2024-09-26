function [noisy_signal] = awgn_channel(signal, SNR)
    
    % Convert SNR from dB to linear
    SNRlin = 10^(SNR/10);
    var_w = 1/SNRlin;

    % Add AWGN
    noise = randn(length(signal),1)*sqrt(var_w/2) + 1i*randn(length(signal),1)*sqrt(var_w/2);
    noisy_signal = signal + noise;
    
end

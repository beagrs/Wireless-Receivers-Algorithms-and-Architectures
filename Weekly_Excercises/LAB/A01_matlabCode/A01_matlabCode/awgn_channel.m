function [bit_out, noisy_signal] = awgn_channel(signal, image_size, SNR)
    
    % Convert SNR from dB to linear
    SNRlin = 
    
    % Add AWGN
    noisy_signal = 
    
    % Demap
    bit_out = demapper(noisy_signal);
    
    % Decode and shown image
    image = image_decoder(bit_out, image_size);
    
end

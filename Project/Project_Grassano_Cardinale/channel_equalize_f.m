function [rx_matrix_equalized] = channel_equalize_f(rx_matrix, conf)
% CHANNEL_EQUALIZE_F: in this function we perform the estimation and
% correction of the channel.
% Channel is estimated using the training of the OFDM frame
% INPUT: rx_matrix = matrix storing in each column an OFDM symbol.
%        conf = structure storing all our global variables
% OUTPUT: rx_matrix_equalized = matrix storing in each column the OFDM data
%         symbols after equalization.
% In our first approximation the first symbol will always be the training

nb_tot_symbs = conf.OFDM_symbs_per_frame + conf.nb_training_symbs;  % Total number of symbols in frame (data + training)
training = (conf.OFDM_training_symbol);
h_hat = zeros(conf.nb_subcarriers, conf.nb_training_symbs);

% Block Channel Case

if strcmp(conf.channel_type,'Block')

    % At every loop we work on a training symbol +  all the data symbols between it
    % and the subsequent training (or the end of the frame). Each time we
    % use the training to estimate the channel and we correct the data
    % symbols by removing the channel
    
    for i = 1 : conf.nb_symbs_between_training + 1 : nb_tot_symbs
      
        % Extract the training
        training_rx = rx_matrix(:, i);
        
        % Divide the rx training by the original one and get the estimated channel
        H_hat = training_rx ./ training;
        h_hat(:,i) = fftshift(ifft(H_hat));
        %jj = jj+1;
        % Perform the channel removing from the symbols
        cols = i : min(i + conf.nb_symbs_between_training, ...
        nb_tot_symbs + (i + conf.nb_symbs_between_training == nb_tot_symbs) * realmax);
        % cols of the equalized matrix that will store the corrected data

        rx_matrix_equalized(:, cols) = rx_matrix(:, cols) .* exp(-1i * mod(angle(H_hat), 2*pi)) ./ abs(H_hat);

    end

elseif strcmp(conf.channel_type,'Block_Viterbi') 
    % Block Channel with Phase noise case: Viterbi - Viterbi
    
    % In this case we need to keep track of the rotation given to the symbols
    % by the phase noise

    % Each row of the matrix contains the vector with the shifting. Use it to
    % correct Viterbi - Viterbi
    shift = zeros(conf.nb_subcarriers, 6) + pi/2*(-1:4);
    
    % Initilize a matrix of phases: each column will store the values of the
    % estimated phase for each subcarrier of the corresponding OFDM symbol.
    % E.g. matrix(1,1) stores the phase shifting of the 1st subacarrier of the
    % 1st OFDM symbol
    theta_hat_matrix = zeros(conf.nb_subcarriers, nb_tot_symbs);
    
    % Loop over every set of symbols training + data before next training
    % (or end)
    for i = 1 : conf.nb_symbs_between_training + 1 : nb_tot_symbs

        % Extract the training
        training_rx = rx_matrix(:, i);
        
        % Divide the rx training by the original one and get the estimated channel
        H_hat = training_rx ./ training;
        h_hat(:, i) = fftshift(ifft(H_hat)); 

         % Remove channel on the training
        rx_matrix_equalized(:,i) = rx_matrix(:,i) .* exp(-1i * mod(angle(H_hat), 2*pi)) ./ abs(H_hat);
        
        % Add the phase shift of the training to the theta_hat_matrix
        theta_hat_matrix(:,i) = mod(angle(H_hat), 2*pi);
        for j = i+1 : min(i + conf.nb_symbs_between_training, nb_tot_symbs + (i + conf.nb_symbs_between_training == nb_tot_symbs)*realmax)
            
            % Take the previous phase shift
            theta_hat_prev = theta_hat_matrix(:,j-1); 
        
            % Generate a matrix storing in each column the vector theta_hat_prev.
            % This will be used a few lines below in the code for the comparison
            % between the estimated phase in all its shifted versions and the
            % estimated phase of each symbol
            theta_hat_prev_matrix = ones(conf.nb_subcarriers, 6) .* theta_hat_prev;
            
            % Estimate the new phase in the [-pi/4, pi/4] interval (QPSK modulation)
            theta_hat = (1/4) * angle(-(rx_matrix(:,j).^4));
        
            % Produce a matrix storing in each row the shifted versions of theta_hat
            % of a specific subcarrier
            theta_hat_shifted = shift + theta_hat;
            
            % Extract a vector containing the theta_hat that is closer to the one
            % of the previous symbol for every subcarrier
            [min_diff, idx] = min(abs(theta_hat_shifted - theta_hat_prev_matrix), [], 2);
            
            for k = 1:conf.nb_subcarriers
                theta_hat_correct(k,1) = theta_hat_shifted(k, idx(k));
            end

            % Filter the phase obtained to take into account noise and rescale to
            % the [0, 2*pi] interval. Add the result to the the culumn of
            % theta_hat_matrix
            theta_hat_matrix(:,j) = mod(0.01*theta_hat_correct + 0.99*theta_hat_prev, 2*pi);
          
            % Perform the channel removing from the rx OFDM symbol
            rx_matrix_equalized(:,j) = rx_matrix(:,j) .* exp(-1i * theta_hat_matrix(:,j)) ./ abs(H_hat);
        end
    end
elseif strcmp(conf.channel_type,'Comb') 
    
    % Loop over every column except the last one (we will deal with that later)
    for i = 1 : size(rx_matrix, 2)
        symb = rx_matrix(:, i);
        
        % Extract the training symbs
        Y = symb(1 : conf.comb_insertion_rate + 1 : end);
        
        % Compute the channel impulse response for the tones in which we
        % placed the training
        %h_hat = ((F'*T') * (T*F)) \ (F'*T' * Y);

        % Perform fft to get the channel transfer function
        H_hat = Y ./ conf.training_comb;

        % Perform the correction for the data symbols subsequent to a
        % training symbol and the training symbol itself
        idxs = 1 : conf.comb_insertion_rate + 1: size(rx_matrix, 1);    % Vector storing the indexes of the training symbs
        symb_equalized = [];
        for j = 1 : length(H_hat) - 1
            
            tmp = symb(idxs(j) : idxs(j+1)-1) .* exp(-1i * mod(angle(H_hat(j)), 2*pi)) ./ abs(H_hat(j));
            symb_equalized = [symb_equalized ; tmp];
        end

        tmp = symb(idxs(end) : end) .* exp(-1i * mod(angle(H_hat(end)), 2*pi)) ./ abs(H_hat(end)); 
        symb_equalized = [symb_equalized ; tmp];

        rx_matrix_equalized(:,i) = symb_equalized;

    end
    
    
%% Plot the impulse response

if (strcmp(conf.channel_type,'Block') || strcmp(conf.channel_type,'Block_Viterbi'))

    time = 0 : (1/conf.f_s) : conf.nb_subcarriers/conf.f_s - 1/conf.f_s;
    
    figure
    plot(time, abs(h_hat).^2);
    xlabel('Time, s');
    ylabel('Channel Magnitude')
    title('Channel Impulse Response')
    xlim([0 0.02])
    grid on;
    axis square; box on;
end
%%  Plot the channel magnitude and phase

if (strcmp(conf.channel_type,'Block') || strcmp(conf.channel_type,'Block_Viterbi'))
    frequency = 0 : conf.BW_BB / conf.nb_subcarriers : conf.BW_BB * (1 - 1/conf.nb_subcarriers);
    frequency = frequency + conf.f_c - conf.BW_BB / 2;
    
    figure
    subplot(2,1,1); plot(frequency, 20*log10(abs(H_hat)./ max(abs(H_hat)) ));
    xlabel('Frequency, Hz');
    ylabel('Magnitude, dB')
    title('Channel Magnitude over Frequency');
    grid on;
    subplot(2,1,2); plot(frequency, unwrap(angle(H_hat)));
    xlabel('Frequency, Hz');
    ylabel('Phase, rad')
    title('Channel Phase over Frequency');
    grid on;
end
end


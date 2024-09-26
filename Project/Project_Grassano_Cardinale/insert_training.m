function [frame_vec, training_comb] = insert_training(data_symbs,conf)
%INSERT_TRAINING: this function is used to insert the training sequences
%inside the OFDM sequences
% INPUT:  data_symbs = column vector storing the data symbols (not OFDM) of
%         my frame
%         conf = struct containing all the configuration variables
% OUTPUT: frame_vec = column vector storing serially all the non-OFDM symbols (data +
%         training) contained in the OfDM frame.
%         training_comb = vector storing the training symbols used for comb, taken from the initial 
%         BPSK training sequence storing nb_subcarriers symbols

% Block training type
if strcmp(conf.training_type,'Block')
    training_comb = 0;
    frame_len = conf.nb_training_symbs * length(conf.OFDM_training_symbol) + length(data_symbs);
    frame_vec = zeros(frame_len, 1);
    for i = 1 : conf.nb_training_symbs

        % Compute the starting index of each training sequence
        idx = (i-1) * conf.nb_symbs_between_training * conf.nb_subcarriers + 1 + (i-1) * length(conf.OFDM_training_symbol);

        % Find the data symbols that lie between the training we are
        % inserting and the (possible) following training
        idx2 = idx - (i-1) * length(conf.OFDM_training_symbol);
        data = data_symbs(idx2: min(idx2 + conf.nb_symbs_between_training * conf.nb_subcarriers -1, ...
            length(data_symbs) + (idx2 + conf.nb_symbs_between_training * conf.nb_subcarriers -1 == length(data_symbs) * realmax)));
        
        % Add the training +  the following data symbs into the output vector
        frame_vec(idx: min(idx + conf.nb_symbs_between_training * conf.nb_subcarriers + length(conf.OFDM_training_symbol) -1, frame_len), 1) = ...
        [conf.OFDM_training_symbol ; data];
    end

% Implement the comb training type. It is possible to insert N BPSK training
% symbols per OFDM frame. The training symbols are chosen starting from the
% ones for block training in the following way: if in comb the symbol is modulated
% into the i-th subcarrier than it will be taken as the i-th symbol of the
% BPSK training vector. This method allows to keep a random-ish aapproach
% in the choice of the training symbol, allowing for more diversity.

elseif strcmp(conf.training_type,'Comb')
    frame_vec = [];
    % Extract some BPSK training symbs to be used interleaved between data
    % symbs following it before the subsequent training (or the end)
    training_comb = conf.OFDM_training_symbol(1 : conf.comb_insertion_rate + 1 : end);
    idxs = 1 : conf.comb_insertion_rate + 1 : conf.nb_subcarriers;
    idx2 = 0;
   k = 1;
   tmp = [];
    while idx2 < length(data_symbs)
        t = training_comb(k);
        idx = idx2 +1;
        idx2 = min(idx + conf.comb_insertion_rate - 1, length(data_symbs) +...
            (idx + conf.comb_insertion_rate - 1 == length(data_symbs))*realmax);
        data = data_symbs(idx : idx2);
        if(length(tmp) + length(data) + 1 <= conf.nb_subcarriers)
            data_symbs(idx : idx2) = 0;
            tmp = [tmp; t ; data];
            
        else
            len_d = conf.nb_subcarriers - length(tmp) - 1;  % The number of data I can insert in that symbol
            tmp = [tmp ; t ; data(1 : len_d)];
            idx2 = idx + len_d - 1;
            data_symbs(idx : idx2) = 0;
        end
        
        if(length(tmp) == conf.nb_subcarriers || idx2 == length(data_symbs))
            frame_vec = [frame_vec ; tmp];
            tmp = [];
        end
        if(k == length(training_comb))
            k = 1;
        else
            k = k +1;
        end

    end
    
end
end


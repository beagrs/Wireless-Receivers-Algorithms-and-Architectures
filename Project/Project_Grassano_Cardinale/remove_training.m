function [data_only_matrix] = remove_training(frame_matrix, conf)
% REMOVE_TRAINING: function used to remove the training sequences from the
% the frame according to the type of training one wants to implement
% INPUT: frame_matrix = OFDM frame in form of a matrix: each column is an OFDM symbol
%        conf = configuration structure
% OUTPUT: data_only_matrix = matrix containing only data symbols (no training)

nb_tot_symbs = conf.OFDM_symbs_per_frame + conf.nb_training_symbs;
if strcmp(conf.training_type,'Block') 
    j = 1;  % Index to choose the columns of the output matrix that will store the data symbs after the training we are considering
    for i = 1 : conf.nb_symbs_between_training + 1 : nb_tot_symbs
        
        % Columns storing the data symbols after the considered training
        cols = j : min(j + conf.nb_symbs_between_training - 1, ...
            conf.OFDM_symbs_per_frame + (j + conf.nb_symbs_between_training - 1 == conf.OFDM_symbs_per_frame)*realmax);
        
        % Assign data to the new matrix
        data_only_matrix(: , cols) = frame_matrix(:, (i+1) : min(i + conf.nb_symbs_between_training, ...
            nb_tot_symbs + (i + conf.nb_symbs_between_training == nb_tot_symbs) * realmax));

        % Update j
        j = j + conf.nb_symbs_between_training;
    end

elseif strcmp(conf.training_type,'Block_Single')   
    
    data_only_matrix = frame_matrix(:, 2 : end);

elseif strcmp(conf.training_type,'Comb') 
    
    idxs = 1 : conf.comb_insertion_rate + 1 : conf.nb_subcarriers;
    data_only_matrix = [];
    for i = 1:  length(idxs) - 1
        tmp = frame_matrix(idxs(i) + 1 : idxs(i+1) - 1, :);
        data_only_matrix = [data_only_matrix ; tmp];
    end
    if(idxs(end) ~= size(frame_matrix, 1))
        tmp = frame_matrix(idxs(end) + 1 : end, :);
        data_only_matrix = [data_only_matrix ; tmp];
    end
end

end


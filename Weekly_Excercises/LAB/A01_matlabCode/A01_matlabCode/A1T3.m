SNR_range = -6:2:12;
L = 1e4;

rng(123)

% Initialize
BER_list_Gray = zeros(size(SNR_range));
BER_list_NonGray = zeros(size(SNR_range));
    




    
for ii = 1:numel(SNR_range) 
    % Convert SNR from dB to linear
    SNRlin = 
    
    % Generate source bitstream
    source = randi([0 1],L,2);
       
    % Map input bitstream using Gray mapping
    mappedGray = 
    if ii == 1
        mappedGray_record = mappedGray;
    end
      
    % Add AWGN
    mappedGrayNoisy = add_awgn_solution(mappedGray, SNRlin);
        
    % Demap
    
    demappedGray = 
    if ii == 1
        demappedGray_record = demappedGray;
    end
        
    % BER calculation for Gray mapping
    BER_list_Gray(ii) = mean(source(:) ~= demappedGray(:));
        
    % Map input bitstream using non-Gray mapping
    mappedNonGray = 
    if ii == 1
        mappedNonGray_record = mappedNonGray;
    end
          
          
    % Add AWGN
    mappedNonGrayNoisy = add_awgn_solution(mappedNonGray, SNRlin);
        
    % Demap
    
    demappedNonGray = 
    if ii == 1
        demappedNonGray_record = demappedNonGray;
    end
        
    % BER calculation for Gray mapping
    BER_list_NonGray(ii) = mean(source(:) ~= demappedNonGray(:));
end


%% uncomment this part for plot
%% graphical ouput
%figure;
% function with input (SNR_range, BER_list_Gray, 'bx-' ,'LineWidth',3)
%
%hold on
%function with input (SNR_range, BER_list_NonGray, 'r*--','LineWidth',3);

%xlabel('SNR (dB)')
%ylabel('BER')
%legend('Gray Mapping', 'Non-Gray Mapping')
%grid on


%% optional: you can write your function for Map/Demap here
% function output = mapGrayfunc(input)
% end
% function output = demapGrayfunc(input)
% end
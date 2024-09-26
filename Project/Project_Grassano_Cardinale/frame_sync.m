function beginning_of_data = frame_sync(rx_signal, conf)
    
   L = conf.preamble_os_factor;
   detection_threshold = 15;
   pulse = rrc(conf.preamble_os_factor, conf.alfa, conf.tx_filter_length * conf.preamble_os_factor);
%   phase_of_peak = 0;
%   magnitude_of_peak = -1;
%    rx_signal_mf_filterd = conv(pulse, rx_signal,'full');
%    rx_signal_mf_filterd = rx_signal_mf_filterd(1+2*conf.tx_filter_length : end - conf.tx_filter_length);
   
   % Filter it
   rx_signal_mf_filtered = conv(rx_signal, pulse, 'same');
   preamble_bpsk = 1 - 2*lfsr_framesync(conf.npreamble);
 
   current_peak_value = 0;
   samples_after_threshold = L;
   vecT = zeros(1,length(L * conf.npreamble + 1 : length(rx_signal_mf_filtered)));
   vecC = zeros(1,length(L * conf.npreamble + 1 : length(rx_signal_mf_filtered)));
   for i = L * conf.npreamble + 1 : length(rx_signal_mf_filtered)
       r = rx_signal_mf_filtered(i - L * conf.npreamble : L : i - L);
       %r = conv(r, pulse,'full');
       c = preamble_bpsk' * r; %r(1 + conf.tx_filter_length : end - conf.tx_filter_length);
       T = abs(c)^2 / abs(r' * r);
       vecT(i) = T;
       vecC(i) = c;
       if (T > detection_threshold || samples_after_threshold < L)
           samples_after_threshold = samples_after_threshold - 1;
           if (T > current_peak_value)
               beginning_of_data = i;
%               phase_of_peak = mod(angle(c),2*pi);
%               magnitude_of_peak = norm(c)/conf.npreamble;
               current_peak_value = T;
           end
           if (samples_after_threshold == 0)
               figure
               plot(abs(vecT))
%                hold on
%                plot(abs(vecC))
        	   return;
           end
       end
   end
   figure
   plot(abs(vecT))
%    hold on
%    plot(abs(vecC))
   if samples_after_threshold == L 
       fprintf("Error\n");
   end
end

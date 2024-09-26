
function filtered_signal = matched_filter(signal, os_factor, mf_length)

    rolloff_factor = 0.22; 
    h = rrc(os_factor, rolloff_factor, mf_length);
    plot(h);
    filtered_signal = conv(h, signal);

end

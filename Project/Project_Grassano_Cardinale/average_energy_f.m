function [average_energy] = average_energy_f(constellation)
%Compute the average energy of a constellation or signal

average_energy = sum(abs(constellation).^2)/length(constellation);

end


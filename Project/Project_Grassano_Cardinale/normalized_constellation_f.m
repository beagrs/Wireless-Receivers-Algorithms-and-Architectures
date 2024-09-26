function [constellation_norm] = normalized_constellation_f(constellation, average_energy_const)
%Output a normalized version of a constellation: s.t. its average energy
%becomes 1
constellation_norm = (1/sqrt(average_energy_const))*constellation;
end


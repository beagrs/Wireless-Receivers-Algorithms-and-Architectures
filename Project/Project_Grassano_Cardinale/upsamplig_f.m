function [symbol_up] = upsamplig_f(os_factor,symbol)
%upsampling a sequence of symbols by putting (os_factor - 1) symbols in i
%for a QPSK
%symbol = sequence of symbols
len_symbol = length(symbol);
symbol_up = zeros(os_factor*len_symbol,1).';
symbol_up(1:os_factor:length(symbol_up)) = symbol;
end


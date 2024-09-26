function [preamble] = preamble_generate(length)
% preamble_generate() 
% input : length: a scaler value, desired length of preamble.
% output: preamble: preamble bits
preamble = zeros(length, 1);
% TODO:
%create buffer
buffer_gen = ones(8,1);

for k = 1:100
    preamble(k) = buffer_gen(8);
    pre_zero = xor(xor(xor(buffer_gen(8), buffer_gen(6)),buffer_gen(5)), buffer_gen(4));
    buffer_gen = circshift(buffer_gen,1)
    buffer_gen(1) = pre_zero;
end

end



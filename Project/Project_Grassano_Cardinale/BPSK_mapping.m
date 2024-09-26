function [BPSK_symbs] = BPSK_mapping(in_bits)
%BPSK_MAPPING: convert bits to BPSK symbols
%in_bits: vector of bits to be mapped
%BPSK_symbs: column vector of mapped bits

BPSK_symbs = 2*in_bits - 1;
end


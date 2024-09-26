clc
clear all
close all

load('task2.mat', 'image_size', 'signal'); % load variable 'image_size' and 'signal'
num_sample_img = image_size(1)* image_size(2)*4;
detector_threshold=30;

preamble = preamble_gen_solution(100);
preamble_bpsk = -2*(preamble) + 1;

start = detector_solution(preamble_bpsk, signal, detector_threshold);
% read the payload here
%main issue: every pixel = 8bit
payload = signal(start:start+prod(image_size)*4-1);

image = image_decoder(demapper(payload), image_size);
imshow(image)

%threshold leading to false positive?
%
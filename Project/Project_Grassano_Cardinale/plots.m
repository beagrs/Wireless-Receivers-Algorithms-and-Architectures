clear all
close all
clc

%% BER vs Training Rate

nb_training = [1, 2, 5, 10, 20, 50];
ber = [0.0065, 0.0042, 0.0044, 0.0032, 0.0027, 0.0024];
figure(1)
semilogy(nb_training, ber, 'b.-', 'LineWidth', 2, 'MarkerSize', 20);
ylabel('BER')
xlabel('Number of training symbols');
title('BER vs Training symbols for 100 OFDM data symbols');
axis square; box on;
grid on

%% BER vs N Symbs. Comparison Viterbi and non Viterbi

nb_symbs = [100 200 300];
ber_nonViterbi = [0.0055, 0.0082, 0.0089];
ber_Viterbi    = [0.0038, 0.0044, 0.0057];

figure(2)
semilogy(nb_symbs, ber_nonViterbi, 'b.-', 'LineWidth', 2, 'MarkerSize', 20);
hold on
semilogy(nb_symbs, ber_Viterbi, 'r.-', 'LineWidth', 2, 'MarkerSize', 20);
xlabel('Number of data symbols')
ylabel('BER')
title('BER vs Number of symbols sent with and without Viterbi')
legend('No Viterbi', 'Viterbi');
grid on




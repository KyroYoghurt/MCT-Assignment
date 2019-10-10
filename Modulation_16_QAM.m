clc;
clear all;
close all;
N=10^6; % Number of Bits to be processed
M=16; %Modulation Order
k=log2(M); %No of bits per symbol
nSymbols=N/k; %No of symbols
generatedata= randi([0 1],N,1); % Generate bits as 0 & 1
Modulator = comm.GeneralQAMModulator; % Define modulator for 16 QAM
Demodulator = comm.GeneralQAMDemodulator;   % Define demodulator for 16 QAM
modulatedata = step(Modulator,generatedata);  % Modulate bits using 16 QAM
E_bN_0=-5:10;   % Define range of Eb/N0 for calculations
SNR = E_bN_0 + 10*log10(k); %Conversion into SNR
for j=1:length(SNR)
receivedSig(:,j)= awgn(modulatedata,SNR(j),'measured');    % Add AWGN noise to modulated Data     
Demodulatedata = step(Demodulator,receivedSig(:,j)); % Demodulate noisy data 
[numerr(j),errrate(j)] = biterr(generatedata,Demodulatedata); %Find number of Error Bits and Error Rate
end
tber = berawgn(E_bN_0,'qam',16,'nondiff');   % Theoretical BER of 16 QAM in AWGN Channel 
semilogy(E_bN_0,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
hold on
semilogy(E_bN_0,errrate,'b-','linewidth',2)  %Plot BER of Simulated Data
axis([-6 11 10^-5 1])
legend ('Theoretical BER','Simulation BER');
xlabel('E_b/N_0 (dB)')
ylabel('Bit Error Rate')
title('16 QAM Modulation in AWGN')
grid on
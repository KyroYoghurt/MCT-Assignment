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
figure(1)
semilogy(E_bN_0,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
hold on
semilogy(E_bN_0,errrate,'b-','linewidth',2)  %Plot BER of Simulated Data
axis([-6 11 10^-5 1])
legend ('Theoretical BER','Simulation BER');
xlabel('E_b/N_0 (dB)')
ylabel('Bit Error Rate')
title('16 QAM Modulation in AWGN')
grid on


M = 16;                         % Modulation order
x = (0:15)';                    % Integer input
y1 = qammod(x,16,'bin');        % 16-QAM output
%Use the scatterplot function to plot the constellation diagram and annotate it with binary representations of the constellation points.
figure(2)
scatterplot(y1)
text(real(y1)+0.1, imag(y1), dec2bin(x))
title('16-QAM, Binary Symbol Mapping')
axis([-4 4 -4 4])

%PSD
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t)+randn(size(t));
figure(3)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
%% Set the simulation parameters.
M = 16;         % Modulation order
k = log2(M);   % Bits per symbol
 EbNo = 5;      % Eb/No (dB)
Fs = 340;       % Sample rate (Hz)
nsamp = 128;     % Number of samples per symbol
freqsep = 20;  % Frequency separation (Hz)
ebnoVec = 1:16;
%% Generate random data symbols.

data = randi([0 M-1],5000,1);

%% Apply FSK modulation.
txsig = fskmod(data,M,freqsep,nsamp,Fs);

%% Pass the signal through an AWGN channel


rxSig  = awgn(txsig,EbNo+10*log10(k)-10*log10(nsamp),'measured',[],'dB');
%% Demodulate the received signal.
dataOut = fskdemod(rxSig,M,freqsep,nsamp,Fs);
%% Calculate the bit error rate.
[num,BER] = biterr(data,dataOut);


%% Determine the theoretical BER and compare it to the estimated BER.
BER_theory = berawgn(EbNo,'fsk',M,'noncoherent');
[BER BER_theory]


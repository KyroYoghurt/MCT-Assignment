custMap = [0 2 4 6 8 10 12 14 15 13 11 9 7 5 3 1];

pskModulator = comm.PSKModulator(16,'BitInput',true,'SymbolMapping','Custom', 'CustomSymbolMapping',custMap);
pskDemodulator = comm.PSKDemodulator(16,'BitOutput',true,'SymbolMapping','Custom','CustomSymbolMapping',custMap);


constellation(pskModulator)
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t)+randn(size(t));
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));

errorRate = comm.ErrorRate;

ebnoVec = 6:18;
ber = zeros(size(ebnoVec));

for k = 1:length(ebnoVec)
    
    % Reset the error counter for each Eb/No value
    reset(errorRate)
    % Reset the array used to collect the error statistics
    errVec = [0 0 0];
    % Set the channel Eb/No
    awgnChannel.EbNo = ebnoVec(k);
    
    while errVec(2) < 200 && errVec(3) < 1e7
        % Generate a 1000-symbol frame
        data = randi([0 1],4000,1);
        % Modulate the binary data
        modData = pskModulator(data);
        % Pass the modulated data through the AWGN channel
        rxSig = awgnChannel(modData);
        % Demodulate the received signal
        rxData = pskDemodulator(rxSig);
        % Collect the error statistics
        errVec = errorRate(data,rxData);
    end
    
    % Save the BER data
    ber(k) = errVec(1);
end


berTheory = berawgn(ebnoVec,'psk',16,'nondiff');

figure
semilogy(ebnoVec,[ber; berTheory])
xlabel('Eb/No (dB)')
ylabel('BER')
grid
legend('Simulation','Theory','location','ne')
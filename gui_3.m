function varargout = gui_3(varargin)
% GUI_3 MATLAB code for gui_3.fig
%      GUI_3, by itself, creates a new GUI_3 or raises the existing
%      singleton*.
%
%      H = GUI_3 returns the handle to a new GUI_3 or the handle to
%      the existing singleton*.
%
%      GUI_3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_3.M with the given input arguments.
%
%      GUI_3('Property','Value',...) creates a new GUI_3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_3

% Last Modified by GUIDE v2.5 17-Oct-2019 10:35:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_3_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_3 is made visible.
function gui_3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_3 (see VARARGIN)

% Choose default command line output for gui_3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in random_bits_gen.
function random_bits_gen_Callback(hObject, eventdata, handles)
% hObject    handle to random_bits_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of random_bits_gen
%function x = randBits(N)
N = 256 ;
x = zeros(1,N);
a = rand(1,N);
for i = 1:N
    if a(i)>0.5;
        x(i) = 1;
    end
end
disp(x);
handles.source = x;
guidata(hObject, handles);
%end

% --- Executes on button press in text2bin_gen.
function text2bin_gen_Callback(hObject, eventdata, handles)
% hObject    handle to text2bin_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of text2bin_gen
%function bits = textinput(fileName)
str = fileread('test1.txt');
ascii = unicode2native(str);
binary = dec2bin(ascii);
[i,l] = size(binary);
bits = zeros(1, i*l);
k = 1;
while k<=i
    a = (k-1)*7 + 1;
    bits(a:a+6) = binary(i,:);
    k = k+1;
end
bits = bits-48;
disp(bits);
handles.source = bits;
guidata(hObject, handles);
%end

% --- Executes on button press in img2bits_gen.
function img2bits_gen_Callback(hObject, eventdata, handles)
% hObject    handle to img2bits_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of img2bits_gen
img = imread("tiger_grayscale.jpg");
a = dec2bin(img);

[i l] = size(a);
a = a';
a = reshape(a,[1,i*l]);
disp(a);
handles.size = size(a);
handles.source = a;
guidata(hObject, handles);
% --- Executes on button press in aud2bits_gen.
function aud2bits_gen_Callback(hObject, eventdata, handles)
% hObject    handle to aud2bits_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aud2bits_gen

%function y = audioBits(audiofilename)
[y,Fs] = audioread('Audio-1.mp3');
nbits = 8;
data = y.*(2^(nbits-1));
data = data(:);
data_uint = typecast(data,'uint8');
binary = dec2bin(data_uint);
binary_trans = binary';
y = binary_trans-'0';
y = reshape(y,1,[]);
disp(y);
handles.source = y;
guidata(hObject, handles);
%end


% --- Executes on button press in hamming_encode.
function hamming_encode_Callback(hObject, eventdata, handles)
% hObject    handle to hamming_encode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hamming_encode
%function [y] = Hamming_code(x) 
x = handles.source;
n = 7 %codeword bits
k = 4 %message bits
A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
G = [ eye(k) A ] %Generator matrix
H = [ A' eye(n-k) ] %Parity-check matrix
%x = [1 0 1 1 0 1 0 1 1 1 1 1]
len = length(x)-3;
j=1;
%ENCODING%
for i=1:4:len
    inp = x(1,i:i+3)
    msg = inp; %Message block vector-change to any 4 bit sequence
    code = mod(msg*G,2)%Encode message
    y(1,j:j+6) = code;
    j = j+7;
end
disp(y); %coded bit stream
handles.code = y;
guidata(hObject, handles);
% --- Executes on button press in bch_encode.
function bch_encode_Callback(hObject, eventdata, handles)
% hObject    handle to bch_encode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bch_encode
% m = 7;
% n = (2^m)- 1;
% k=64;
% m = handles.source; %msg=input('Enter 64 bit message:');
% msg = gf(m);
% disp('Message:');
% disp (msg);
% % t the error correction capability
% [genpoly,t] = bchgenpoly(n,k);
% disp('Error Correcting Capability:');
% disp(t);
% %the coding
% code = bchenc(msg,n,k);
% disp('encoded message:');
% c=gf(code);
% disp(c);

% --- Executes on button press in convolution_encode.
function convolution_encode_Callback(hObject, eventdata, handles)
% hObject    handle to convolution_encode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convolution_encode

%function [y tblen] = ConvolutionEncode(x)
x = handles.source;
K=3;
G1=5;
G2=7;
tblen = length(x);
%msg=[1 1 0 0 1 0]
trel=poly2trellis(K,[G1 G2]);
y=convenc(x,trel);
disp(y);
handles.code = y;
guidata(hObject, handles);



function ebno_user_input_Callback(hObject, eventdata, handles)
% hObject    handle to ebno_user_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebno_user_input as text
%        str2double(get(hObject,'String')) returns contents of ebno_user_input as a double
ebno = get(handles.ebno_user_input, 'String');
ebno = str2num(ebno);
disp(ebno);
handles.ebno = ebno;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ebno_user_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebno_user_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fsk_mod.
function fsk_mod_Callback(hObject, eventdata, handles)
% hObject    handle to fsk_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fsk_mod
%% Set the simulation parameters.
M = 16;         % Modulation order
k = log2(M);   % Bits per symbol
 EbNo = 5;      % Eb/No (dB)
Fs = 340;       % Sample rate (Hz)
nsamp = 128;     % Number of samples per symbol
freqsep = 20;  % Frequency separation (Hz)
ebnoVec = 1:16;
%% Generate random data symbols.

data1 = handles.code;
data2 = handles.source;
%% Apply FSK modulation.
txsig1 = fskmod(data1,M,freqsep,nsamp,Fs);
txsig2 = fskmod(data2,M,freqsep,nsamp,Fs);

Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = txsig1;
figure(1)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
    
    
t1=0:length(txsig1)-1;
figure(2)
plot(t1,txsig1);
title('Time domain plot');
%% Pass the signal through an AWGN channel


rxSig1  = awgn(txsig1,EbNo+10*log10(k)-10*log10(nsamp),'measured',[],'dB');
rxSig2  = awgn(txsig2,EbNo+10*log10(k)-10*log10(nsamp),'measured',[],'dB');



%% Demodulate the received signal.
dataOut1 = fskdemod(rxSig1,M,freqsep,nsamp,Fs);
dataOut2 = fskdemod(rxSig2,M,freqsep,nsamp,Fs);


%% Calculate the bit error rate.
[num1,BER1] = biterr(data1,dataOut1);
[num2,BER2] = biterr(data2,dataOut2);


%% Determine the theoretical BER and compare it to the estimated BER.
BER_theory = berawgn(EbNo,'fsk',M,'noncoherent');


% figure
 semilogy(ebnoVec,[BER1;  BER_theory])
 xlabel('Eb/No (dB)')
 ylabel('BER')
 grid
 legend('Simulation-Coded message','Theory','location','ne')


% --- Executes on button press in psk_mod.
function psk_mod_Callback(hObject, eventdata, handles)
% hObject    handle to psk_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of psk_mod
custMap = [0 2 4 6 8 10 12 14 15 13 11 9 7 5 3 1];

pskModulator = comm.PSKModulator(16,'BitInput',true,'SymbolMapping','Custom', 'CustomSymbolMapping',custMap);
pskDemodulator = comm.PSKDemodulator(16,'BitOutput',true,'SymbolMapping','Custom','CustomSymbolMapping',custMap);


%constellation(pskModulator)

awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
errorRate = comm.ErrorRate;

ebnoVec = 1:30;
ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

for z = 1:length(ebnoVec)
    
    % Reset the error counter for each Eb/No value
    reset(errorRate)
    % Reset the array used to collect the error statistics
    errVec1 = [0 0 0];
    % Set the channel Eb/No
    awgnChannel.EbNo = ebnoVec(z);
    
    while errVec1(2) < 200 && errVec1(3) < 1e7
        % Generate a 1000-symbol frame
        data1 = handles.code;
        if mod(length(data1),4) ~= 0
            for a = 3:-1:mod(length(data1),4)
            data1 = [data1,0] 
            end
        end
        data1 = data1';
        data2 = handles.source;
        if mod(length(data2),4) ~= 0
            for a = 3:-1:mod(length(data2),4)
            data2 = [data2,0] 
            end
        end
        data2 = data2';
        % Modulate the binary data
        modData1 = pskModulator(data1);
        modData2 = pskModulator(data2);
        
        % Pass the modulated data through the AWGN channel
        rxSig1 = awgnChannel(modData1);
        rxSig2 = awgnChannel(modData2);
        
        % Demodulate the received signal
        rxData1 = pskDemodulator(rxSig1);
        rxData2 = pskDemodulator(rxSig2);
        
        % Collect the error statistics
        errVec1 = errorRate(data1,rxData1);
        errVec2 = errorRate(data2,rxData2);
        
    end
    
    % Save the BER data
    ber1(z) = errVec1(1);
    ber2(z) = errVec2(1);
    
end
%Plot PSD of modulated signal
Fs = 10000;
t = 0:1/Fs:1-1/Fs;
x = modData1;
figure(3)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));

berTheory = berawgn(ebnoVec,'psk',16,'nondiff');

figure
semilogy(ebnoVec,[ber1; ber2; berTheory])
xlabel('Eb/No (dB)')
ylabel('BER')
grid
legend('Simulation-Coded message','Simulation-Uncoded message','Theory','location','ne')

M = 16;             % Modulation alphabet size
phOffset = 0;       % Phase offset
symMap = 'binary';  % Symbol mapping (either 'binary' or 'gray')
%Construct the modulator object.

pskModulator = comm.PSKModulator(M,phOffset,'SymbolMapping',symMap);
%Plot the constellation.

constellation(pskModulator);
handles.dem1 = rxData1;
handles.dem2 = rxData2;
guidata(hObject, handles);


data1 = handles.code;
data1 =data1';
% Modulate the binary data
modData1 = pskModulator(data1);

t1=0:length(modData1)-1;
figure(4)
plot(t1,modData1)

% --- Executes on button press in qam_mod.
function qam_mod_Callback(hObject, eventdata, handles)
% hObject    handle to qam_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Constellation diagram
% Hint: get(hObject,'Value') returns toggle state of qam_mod
M = 16;                         % Modulation order
x = (0:15)';                    % Integer input
y1 = qammod(x,16,'bin');        % 16-QAM output
%Use the scatterplot function to plot the constellation diagram and annotate it with binary representations of the constellation points.
%figure(3)
scatterplot(y1)
text(real(y1)+0.1, imag(y1), dec2bin(x))
title('16-QAM, Binary Symbol Mapping')
axis([-4 4 -4 4]);

x1 = handles.code;
x2 = handles.source;
x1 =x1';
x2 = x2';
N1=length(x1); % Number of coded Bits to be processed
N2=length(x2); % Number of uncoded Bits to be processed


M=16; %Modulation Order
k=log2(M); %No of bits per symbol
nSymbols1=N1/k; %No of coded symbols
nSymbols2=N2/k; %No of uncoded symbols

%generatedata= randi([0 1],N,1); % Generate bits as 0 & 1
Modulator = comm.GeneralQAMModulator; % Define modulator for 16 QAM
Demodulator = comm.GeneralQAMDemodulator;   % Define demodulator for 16 QAM
modulatedata1 = step(Modulator,x1);  % Modulate coded bits using 16 QAM
modulatedata2 = step(Modulator,x2);  % Modulate uncoded bits using 16 QAM

E_bN_0=-5:10;   % Define range of Eb/N0 for calculations
SNR = E_bN_0 + 10*log10(k); %Conversion into SNR
for j=1:length(SNR)
receivedSig1(:,j)= awgn(modulatedata1,SNR(j),'measured');    % Add AWGN noise to modulated Data     
Demodulatedata1 = step(Demodulator,receivedSig1(:,j)); % Demodulate noisy coded data 
[numerr1(j),errrate1(j)] = biterr(x1,Demodulatedata1); %Find number of Error Bits and Error Rate
end
tber = berawgn(E_bN_0,'qam',16,'nondiff');   % Theoretical BER of 16 QAM in AWGN Channel 
figure(2)
semilogy(E_bN_0,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
hold on
semilogy(E_bN_0,errrate1,'b-','linewidth',2)  %Plot BER of Simulated Data
axis([-6 11 10^-5 1])
legend ('Theoretical BER','Simulation BER');
xlabel('E_b/N_0 (dB)')
ylabel('Bit Error Rate')
title('16 QAM Modulation in AWGN for coded data')
grid on

for j=1:length(SNR)
receivedSig2(:,j)= awgn(modulatedata2,SNR(j),'measured');    % Add AWGN noise to modulated Data     
Demodulatedata2 = step(Demodulator,receivedSig2(:,j)); % Demodulate noisy coded data 
[numerr2(j),errrate2(j)] = biterr(x2,Demodulatedata2); %Find number of Error Bits and Error Rate
end
tber = berawgn(E_bN_0,'qam',16,'nondiff');   % Theoretical BER of 16 QAM in AWGN Channel 
figure(3)
semilogy(E_bN_0,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
hold on
semilogy(E_bN_0,errrate2,'b-','linewidth',2)  %Plot BER of Simulated Data
axis([-6 11 10^-5 1])
legend ('Theoretical BER','Simulation BER');
xlabel('E_b/N_0 (dB)')
ylabel('Bit Error Rate')
title('16 QAM Modulation in AWGN for uncoded data')
grid on


Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = receivedSig1;
figure(4)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))




handles.dem1 = Demodulatedata1;
handles.dem2 = Demodulatedata2;
guidata(hObject, handles);

% --- Executes on button press in hamming_decode.
function hamming_decode_Callback(hObject, eventdata, handles)
% hObject    handle to hamming_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hamming_decode
%function [x] = Hamming_Decode(y)
%y = [1 0 1 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 1 1 1];
y = handles.dem1;
y =y';
 n = 7 %codeword bits
 k = 4 %message bits
 A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
 G = [ eye(k) A ] %Generator matrix
 H = [ A' eye(n-k) ] %Parity-check matrix
 len = length(y)-6;
 j=1;
% %DECODING%
 for i=1:7:len
     codewrd = y(1,i:i+6)
%     %codewrd(2)=~codewrd(2);  %Introducing some error
     syndrome = mod(codewrd * H',2);
%     %Find position of the error in codeword (index)
%      find = 0;
%      for ii = 1:n
%          if ~find
%              errvect = zeros(1,n);
%              errvect(ii) = 1;
%              search = mod(errvect * H',2);
%              if search == syndrome
%                  find = 1;
%                  index = ii;
%              end
%          end
%      end
%      disp(['Position of error in codeword=',num2str(index)]);
      correctedcode = codewrd;
%      correctedcode(index) = mod(codewrd(index)+1,2)%Corrected codeword
%     %Strip off parity bits
     msg_decoded=correctedcode;
     msg_decoded=msg_decoded(1:4)
% 
%      msg = inp; %Message block vector-change to any 4 bit sequence
%     code = mod(msg*G,2)%Encode message
     x(1,j:j+3) = msg_decoded;
     j = j+4;
 end
 x
 handles.decode = x;
 guidata(hObject, handles);

% --- Executes on button press in bch_decode.
function bch_decode_Callback(hObject, eventdata, handles)
% hObject    handle to bch_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bch_decode


% --- Executes on button press in convolutional_decode.
function convolutional_decode_Callback(hObject, eventdata, handles)
% hObject    handle to convolutional_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of convolutional_decode
%function [x] = ConvolutionDecode(y, tblen)
trellis = poly2trellis(3,[5 7]);

% Generate random binary data, convolutionally encode the data, and decode the data using the Viterbi algorithm.


codedData = handles.dem1';
decodedData = vitdec(codedData,trellis,34,'trunc','hard')

% --- Executes on button press in bits2text.
 function bits2text_Callback(hObject, eventdata, handles)
% hObject    handle to bits2text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bits2text
%function x = convert2Text(bitstream)
 bitstream = handles.decode;
 l = size(bitstream,2)
 nchar = 7
 x = char(zeros(nchar,1))
 k =1;
 for i = 1:7:l
     c = native2unicode([64 32 16 8 4 2 1].*((bitstream(i:i+6)-48)'));
     x(k) = c;
     k=k+1;
 end
 disp(x);
% end

% --- Executes on button press in bit2image.
% function bit2image_Callback(hObject, eventdata, handles)
% hObject    handle to bit2image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bit2image
%function image = bits2img(bits,x,y,z)
% [x y z] = handles.size;
% bits  = handles.dem2;
% a = zeros(1,size(bits,2)/8);
% k=1;
% for i = 1:8:size(bits,2)
%     c = bin2dec(bits(i:i+7));
%     a(k) = c;
%     k=k+1;
% end
% image = reshape(a,[x,y,z]);
% image = image/255;
% imshow(image);
% %image = uint8(image);
% %end

% --- Executes on button press in bits2audio.
% function bits2audio_Callback(hObject, eventdata, handles)
% hObject    handle to bits2audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bits2audio

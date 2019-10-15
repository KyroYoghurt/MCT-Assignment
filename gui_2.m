function varargout = gui_2(varargin)
% GUI_2 MATLAB code for gui_2.fig
%      GUI_2, by itself, creates a new GUI_2 or raises the existing
%      singleton*.
%
%      H = GUI_2 returns the handle to a new GUI_2 or the handle to
%      the existing singleton*.
%
%      GUI_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_2.M with the given input arguments.
%
%      GUI_2('Property','Value',...) creates a new GUI_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_2

% Last Modified by GUIDE v2.5 16-Oct-2019 00:31:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_2_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_2_OutputFcn, ...
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


% --- Executes just before gui_2 is made visible.
function gui_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_2 (see VARARGIN)

% Choose default command line output for gui_2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in bin2text_reconstruction.
function bin2text_reconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to bin2text_reconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bin2text_reconstruction
%function x = convert2Text(handles.noise)
l = size(handles.noise,2)
nchar = l/7
x = char(zeros(1,nchar))
k =1;
for i = 1:7:l
    c = native2unicode([64 32 16 8 4 2 1]*((handles.noise(i:i+6)-48)'));
    x(k) = c;
    k=k+1;
end
%end
disp(x)
% --- Executes on button press in bits2img_reconstruction.
function bits2img_reconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to bits2img_reconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bits2img_reconstruction


% --- Executes on button press in bin2audio_reconstruction.
function bin2audio_reconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to bin2audio_reconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bin2audio_reconstruction


% --- Executes on button press in ber_with_coding.
function ber_with_coding_Callback(hObject, eventdata, handles)
% hObject    handle to ber_with_coding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ber_with_coding
%% Theoretical Coded BERS

SNR=10.^(EbNo./10);
p=qfunc(sqrt(SNR));  % Given for BPSK. Change the value of p for the three schemes(16-PSK, 16-FSK, 16-QAM)
BER_HARD=zeros(1,length(SNR));
for j=2:n
    BER_HARD=BER_HARD+nchoosek(n,j).*(p.^j).*((1-p).^(n-j));
end
semilogy(SNRdB,BER_HARD,'k-','linewidth',2.0);         % Theroritical BER in decoding

%% Simulated Coded BERS

EbN0=-5:10;   % Define range of Eb/N0 for calculations
trans_msg=handles.source;%[1 0 0 1 1 0 1]; % from function
decoded_msg_hamm=[1 0 0 0 1 0 0]; % from hamming decoding function
decoded_msg_BCH=[1 0 0 1 0 1 0 ]; % from BCH decoding function
decoded_msg_conv=[0 1 0 1 1 1 1]; % from convolutional decoding function

for j=1:length(EbN0)
    
    [numerrhamm(j),errratehamm(j)] = biterr(handles.source,decoded_msg_hamm); %Find number of Error Bits and Error Rate
    [numerrBCH(j),errrateBCH(j)] = biterr(handles.source,decoded_msg_BCH);
    [numerrconv(j),errrateconv(j)] = biterr(handles.source,decoded_msg_conv);

end

figure
semilogy(EbN0,errratehamm,'b-','linewidth',2)  %Plot BER of Simulated Data
hold on;
semilogy(EbN0,errrateBCH,'b-','linewidth',2)  %Plot BER of Simulated Data
semilogy(EbN0,errrateconv,'b-','linewidth',2)  %Plot BER of Simulated Data
hold off;
% --- Executes on button press in ber_without_coding.
function ber_without_coding_Callback(hObject, eventdata, handles)
% hObject    handle to ber_without_coding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ber_without_coding

%% Uncoded Theoretical BERs

EbNo = (0:10)';
figure;
M = 16;
berP = berawgn(EbNo,'psk',M,'nondiff');
berQ = berawgn(EbNo,'qam',M);
berF = berawgn(EbNo,'fsk',M,'coherent');
semilogy(EbNo,[berP berQ berF])
xlabel('Eb/No (dB)')
ylabel('BER')
legend('16-PSK','16-QAM','16-FSK')
grid
% --- Executes on button press in Hamming_decode.
function Hamming_decode_Callback(hObject, eventdata, handles)
% hObject    handle to Hamming_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hamming_decode
%function [x] = Hamming_Decode(y)
%y = [1 0 1 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 1 1 1];

n = 7 %codeword bits
k = 4 %message bits
A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
G = [ eye(k) A ] %Generator matrix
H = [ A' eye(n-k) ] %Parity-check matrix
len = length(y)-6;
j=1;
%DECODING%
for i=1:7:len
    codewrd = y(1,i:i+6)
    %codewrd(2)=~codewrd(2);  %Introducing some error
    syndrome = mod(codewrd * H',2);
    %Find position of the error in codeword (index)
    find = 0;
    for ii = 1:n
        if ~find
            errvect = zeros(1,n);
            errvect(ii) = 1;
            search = mod(errvect * H',2);
            if search == syndrome
                find = 1;
                index = ii;
            end
        end
    end
    disp(['Position of error in codeword=',num2str(index)]);
    correctedcode = codewrd;
    correctedcode(index) = mod(codewrd(index)+1,2)%Corrected codeword
    %Strip off parity bits
    msg_decoded=correctedcode;
    msg_decoded=msg_decoded(1:4)

    %msg = inp; %Message block vector-change to any 4 bit sequence
    %code = mod(msg*G,2)%Encode message
    x(1,j:j+3) = msg_decoded;
    j = j+4;
end
x

% --- Executes on button press in BCH_decode.
function BCH_decode_Callback(hObject, eventdata, handles)
% hObject    handle to BCH_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BCH_decode


% --- Executes on button press in Convolution_decode.
function Convolution_decode_Callback(hObject, eventdata, handles)
% hObject    handle to Convolution_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convolution_decode


% --- Executes on button press in Constellation.
function Constellation_Callback(hObject, eventdata, handles)
% hObject    handle to Constellation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Constellation


% --- Executes on button press in Time_domain.
function Time_domain_Callback(hObject, eventdata, handles)
% hObject    handle to Time_domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Time_domain


% --- Executes on button press in PSD.
function PSD_Callback(hObject, eventdata, handles)
% hObject    handle to PSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PSD


% --- Executes on button press in FSK_encoding.
function FSK_encoding_Callback(hObject, eventdata, handles)
% hObject    handle to FSK_encoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FSK_encoding

%% Set the simulation parameters.
M = 16;         % Modulation order
k = log2(M);   % Bits per symbol
 EbNo = 5;      % Eb/No (dB)
Fs = 340;       % Sample rate (Hz)
nsamp = 128;     % Number of samples per symbol
freqsep = 20;  % Frequency separation (Hz)
ebnoVec = 1:16;
%% Input the channel coded data

data = handles.code;

%% Apply FSK modulation.
txsig = fskmod(data,M,freqsep,nsamp,Fs);

t1=0:length(txsig)-1;
figure(1)
plot(t1,txsig);
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t)+randn(size(t));
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))

handles.mod_sig = abs(txsig);
guidata(hObject, handles);


% --- Executes on button press in QAM_encoding.
function QAM_encoding_Callback(hObject, eventdata, handles)
% hObject    handle to QAM_encoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of QAM_encoding
% clc;
% clear all;
% close all;
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
% figure(1)
% semilogy(E_bN_0,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
% hold on
% semilogy(E_bN_0,errrate,'b-','linewidth',2)  %Plot BER of Simulated Data
% axis([-6 11 10^-5 1])
% legend ('Theoretical BER','Simulation BER');
% xlabel('E_b/N_0 (dB)')
% ylabel('Bit Error Rate')
% title('16 QAM Modulation in AWGN')
% grid on


figure(2)
t1=0:length(modulatedata)-1;
plot(t1,modulatedata)
xlabel('Time')
ylabel('Modulated Signal')
title('Time domain plot for 16 QAM Modulation')

M = 16;                         % Modulation order
x = (0:15)';                    % Integer input
y1 = qammod(x,16,'bin');        % 16-QAM output
%Use the scatterplot function to plot the constellation diagram and annotate it with binary representations of the constellation points.
%figure(3)
scatterplot(y1)
text(real(y1)+0.1, imag(y1), dec2bin(x))
title('16-QAM, Binary Symbol Mapping')
axis([-4 4 -4 4])

% --- Executes on button press in PSK_encoding.
function PSK_encoding_Callback(hObject, eventdata, handles)
% hObject    handle to PSK_encoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PSK_encoding
custMap = handles.code;

pskModulator = comm.PSKModulator(16,'BitInput',true,'SymbolMapping','Custom', 'CustomSymbolMapping',custMap);
pskDemodulator = comm.PSKDemodulator(16,'BitOutput',true,'SymbolMapping','Custom','CustomSymbolMapping',custMap);


%constellation(pskModulator)

Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t)+randn(size(t));
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));

errorRate = comm.ErrorRate;

ebnoVec = 6:18;
ber = zeros(size(ebnoVec));

% for k = 1:length(ebnoVec)
%     
%     % Reset the error counter for each Eb/No value
%     reset(errorRate)
%     % Reset the array used to collect the error statistics
%     errVec = [0 0 0];
%     % Set the channel Eb/No
%     awgnChannel.EbNo = ebnoVec(k);
%     
%     while errVec(2) < 200 && errVec(3) < 1e7
%         % Generate a 1000-symbol frame
%         data = randi([0 1],4000,1);
%         % Modulate the binary data
%         modData = pskModulator(data);
%         % Pass the modulated data through the AWGN channel
%         rxSig = awgnChannel(modData);
%         % Demodulate the received signal
%         rxData = pskDemodulator(rxSig);
%         % Collect the error statistics
%         errVec = errorRate(data,rxData);
%     end
%     
%     % Save the BER data
%     ber(k) = errVec(1);
% end


% berTheory = berawgn(ebnoVec,'psk',16,'nondiff');

% figure
% semilogy(ebnoVec,[ber; berTheory])
% xlabel('Eb/No (dB)')
% ylabel('BER')
% grid
% legend('Simulation','Theory','location','ne')

M = 16;             % Modulation alphabet size
phOffset = 0;       % Phase offset
symMap = 'binary';  % Symbol mapping (either 'binary' or 'gray')
%Construct the modulator object.

pskModulator = comm.PSKModulator(M,phOffset,'SymbolMapping',symMap);
%Plot the constellation.

constellation(pskModulator)

data = randi([0 1],100,1);
% Modulate the binary data
modData = pskModulator(data);
t1=0:length(modData)-1;
figure(4)
plot(t1,modData)

% --- Executes on button press in Hamming_code.
function Hamming_code_Callback(hObject, eventdata, handles)
% hObject    handle to Hamming_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hamming_code
%function [y] = Hamming_code(x) 
x = handles.source ;
n = 7 ;%codeword bits
k = 4 ;%message bits
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
handles.code = y

guidata(hObject, handles);

% --- Executes on button press in BCH_code.
function BCH_code_Callback(hObject, eventdata, handles)
% hObject    handle to BCH_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BCH_code
%the message
% clear;
% clc;
m = 7;
n = (2^m)- 1;
k=64;
m = handles.source; %msg=input('Enter 64 bit message:');
x = gf(m);
%k = length(msg);
disp('Message:');
disp (x);
% t the error correction capability
[genpoly,t] = bchgenpoly(n,k);
disp('Error Correcting Capability:');
disp(t);
%the coding


% disp('encoded message:');
% c=gf(code);
% disp(c);


for i=1:length(x):1
    inp = x(1,i:i+63)
    msg = inp; %Message block vector-change to any 4 bit sequence
    code = bchenc(msg,n,k);%Encode message
    disp('Random')
    c=gf(code);
%disp(c);

%     y(1,j:j+6) = code;
%     j = j+7;
end
c
% --- Executes on button press in Convolution_code.
function Convolution_code_Callback(hObject, eventdata, handles)
% hObject    handle to Convolution_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convolution_code
%function [y tblen] = ConvolutionEncode(x)
x = handles.source
K=3;
G1=5;
G2=7;
tblen = length(x);
%msg=[1 1 0 0 1 0]
trel=poly2trellis(K,[G1 G2]);
y=convenc(x,trel);
handles.code = y;
guidata(hObject, handles);

% --- Executes on button press in Random_bit_gen.
function Random_bit_gen_Callback(hObject, eventdata, handles)
% hObject    handle to Random_bit_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Random_bit_gen
%function x = randBits(N)
N = 256;
x = zeros(1,N);
a = rand(1,N);
for i = 1:N
    if a(i)>0.5;
        x(i) = 1;
    end
end
x
handles.source = x;
guidata(hObject, handles);
%end

% --- Executes on button press in Text2binary.
function Text2binary_Callback(hObject, eventdata, handles)
% hObject    handle to Text2binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Text2binary
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
bits = bits-48
handles.source = bits
guidata(hObject, handles);
%end


% --- Executes on button press in img2binary.
function img2binary_Callback(hObject, eventdata, handles)
% hObject    handle to img2binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of img2binary
img = imread("tiger_grayscale.jpg");
source = dec2bin(img);

[i l] = size(source);
source = source';
source = reshape(source,[1,i*l])
handles.source = source;
guidata(hObject, handles);
% --- Executes on button press in aud2bin.
function aud2bin_Callback(hObject, eventdata, handles)
% hObject    handle to aud2bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
% Hint: get(hObject,'Value') returns toggle state of aud2bin



% --- Executes on button press in snr_10db.
function snr_10db_Callback(hObject, eventdata, handles)
% hObject    handle to snr_10db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of snr_10db





function snr_usr_inp_Callback(hObject, eventdata, handles)
% hObject    handle to snr_usr_inp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr_usr_inp as text
%        str2double(get(hObject,'String')) returns contents of snr_usr_inp as a double
%function x = AddAWGN(vec, snr);
vec = handles.mod_sig;
snr = get(handles.snr_usr_inp, 'String');
snr = str2num(snr)
x = awgn(vec,snr);
handles.noise = x;
handles.snr = snr;
guidata(hObject, handles);
%end

% --- Executes during object creation, after setting all properties.
function snr_usr_inp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr_usr_inp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in covolution_decode.
function covolution_decode_Callback(hObject, eventdata, handles)
% hObject    handle to covolution_decode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of covolution_decode


% --- Executes on button press in FSK_demod.
function FSK_demod_Callback(hObject, eventdata, handles)
% hObject    handle to FSK_demod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FSK_demod
M = 16;         % Modulation order
k = log2(M);
Fs = 340;       % Sample rate (Hz)
nsamp = 128;     % Number of samples per symbol
freqsep = 20;  % Frequency separation (Hz)

txsig = handles.mod_sig;
rxSig  = awgn(txsig,handles.snr+10*log10(k)-10*log10(nsamp),'measured',[],'dB');
dataOut = fskdemod(rxSig,M,freqsep,nsamp,Fs);
disp(dataOut);

% --- Executes on button press in QAM_demod.
function QAM_demod_Callback(hObject, eventdata, handles)
% hObject    handle to QAM_demod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of QAM_demod


% --- Executes on button press in PSK_demod.
function PSK_demod_Callback(hObject, eventdata, handles)
% hObject    handle to PSK_demod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PSK_demod

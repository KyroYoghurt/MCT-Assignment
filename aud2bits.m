[y1, Fs] = audioread('Audio-1.mp3');

y2 = uint8(y1*Fs);
[x y z] = size(y1);
a = dec2bin(y2);
[i l] = size(a);
a = reshape(a,[1,i*l]);



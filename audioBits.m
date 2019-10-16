function y = audioBits(audiofilename)
[y,Fs] = audioread(audiofilename);
nbits = 8;
data = y.*(2^(nbits-1));
data = data(:);
data_uint = typecast(data,'uint8');
binary = dec2bin(data_uint);
binary_trans = binary';
y = binary_trans-'0';
y = reshape(y,1,[]);
end

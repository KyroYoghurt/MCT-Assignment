function x = audioBits(audiofilename)
fileId= fopen(audiofilename,'r');
data = fread(fileId,'uint8');
bits = uint8(rem(floor((1./[128 64 32 16 8 4 2 1]).*data),2));
bits_trans = bits';
[i l] = size(bits);
x = reshape(bits_trans,[1,i*l]);
end

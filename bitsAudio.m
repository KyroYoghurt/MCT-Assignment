function x = bitsAudio(bits)
Fs = 48000;
bin_new = bits';
bin_char = char('0' + bin_new);
rebin_char = reshape(bin_char,8,[]);
bin_char_trans = rebin_char';
dec = uint8(bin2dec(bin_char_trans));
dec_new = typecast(dec,'double');
x = dec_new./(2^7);
x = reshape(x,[],2);
audiowrite('Reconstructed.wav',x,Fs);
sound(x, Fs);
end

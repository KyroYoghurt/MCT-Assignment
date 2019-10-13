function audio = bits2aud(bits,x,y,z,Fs)
a = zeros(1,size(bits,2)/8);
k=1;
for i = 1:8:size(bits,2)
    c = bin2dec(bits(i:i+7));
    a(k) = c;
    k=k+1;
end
audio = reshape(a,[x,y,z]);
audio = audio/255;
sound(audio)

end


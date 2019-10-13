function image = bits2img(bits,x,y,z)
a = zeros(1,size(bits,2)/8);
k=1;
for i = 1:8:size(bits,2)
    c = bin2dec(bits(i:i+7));
    a(k) = c;
    k=k+1;
end
image = reshape(a,[x,y,z]);
image = image/255;
imshow(image);
%image = uint8(image);
end


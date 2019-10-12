function x = convert2Text(bitstream)
l = size(bitstream,2)
nchar = l/7
x = char(zeros(1,nchar))
k =1;
for i = 1:7:l
    c = native2unicode([64 32 16 8 4 2 1]*((bitstream(i:i+6)-48)'));
    x(k) = c;
    k=k+1;
end
end

   

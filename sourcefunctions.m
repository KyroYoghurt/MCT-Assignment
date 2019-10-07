%% call functions in this section
%%randomBITS = randBits(200);
%%textBits = textinput('text.txt');


function x = randBits(N)
x = zeros(1,N);
a = rand(1,N);
for i = 1:N
    if a(i)>0.5;
        x(i) = 1;
    end
end
end

function bits = textinput(fileName)
str = fileread(fileName);
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
end

    
    

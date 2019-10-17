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
bin_trans = binary';
bits = reshape(bin_trans,1,[]);
end

    
    

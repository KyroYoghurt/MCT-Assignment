clear
n = 7 %codeword bits
k = 4 %message bits
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
y %coded bit stream

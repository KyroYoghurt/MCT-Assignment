clear
n = 7 %codeword bits
k = 4 %message bits
A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
G = [ A eye(k)  ] %Generator matrix
H = [ eye(n-k) A'] %Parity-check matrix
%x = [1 0 1 1 0 1 0 1 1 1 1 1]
len = length(x)-3;
%ENCODING%
for i=1:4:len
    inp = x(1,i:i+3)
    msg = inp; %Message block vector-change to any 4 bit sequence
    code = mod(msg*G,2)%Encode message
end

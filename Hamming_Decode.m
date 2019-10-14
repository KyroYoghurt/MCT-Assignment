function [x] = Hamming_Decode(y)
%y = [1 0 1 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 1 1 1];
n = 7 %codeword bits
k = 4 %message bits
A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
G = [ eye(k) A ] %Generator matrix
H = [ A' eye(n-k) ] %Parity-check matrix
len = length(y)-6;
j=1;
%DECODING%
for i=1:7:len
    codewrd = y(1,i:i+6)
    codewrd(2)=~codewrd(2);  %Introducing some error
    syndrome = mod(codewrd * H',2);
    %Find position of the error in codeword (index)
    find = 0;
    for ii = 1:n
        if ~find
            errvect = zeros(1,n);
            errvect(ii) = 1;
            search = mod(errvect * H',2);
            if search == syndrome
                find = 1;
                index = ii;
            end
        end
    end
    disp(['Position of error in codeword=',num2str(index)]);
    correctedcode = codewrd;
    correctedcode(index) = mod(codewrd(index)+1,2)%Corrected codeword
    %Strip off parity bits
    msg_decoded=correctedcode;
    msg_decoded=msg_decoded(1:4)

    %msg = inp; %Message block vector-change to any 4 bit sequence
    %code = mod(msg*G,2)%Encode message
    x(1,j:j+3) = msg_decoded;
    j = j+4;
end
x

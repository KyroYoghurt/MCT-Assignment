%% BCH Encoding

%the message
clear;
clc;
m = 7;
n = (2^m)- 1;
k=64;
m=[1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 0 1 0 0]; %msg=input('Enter 64 bit message:');
msg = gf(m);
disp('Message:');
disp (msg);
% t the error correction capability
[genpoly,t] = bchgenpoly(n,k);
disp('Error Correcting Capability:');
disp(t);
%the coding
code = bchenc(msg,n,k);
disp('encoded message:');
c=gf(code);
disp(c);
noisycode = code + randerr(1,n,1:t);
disp('received codeword:');
disp(noisycode);
%decoding
[newmsg,err,ccode] = bchdec(noisycode,n,k);
disp('decoded codeword:');
disp(newmsg);
if msg==newmsg
    disp('the message was recovered perfectly');
else
    disp('error in recovery of message');
end

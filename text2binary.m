%% MATLAB code for converting text to binary bits
clc 
clear all 
close all
[f_name,path] = uigetfile({'*.txt';'*.*'},'Select the Text (.txt) file'); 
%above line will let you to browse for the file using a GUI window
path=[path f_name]; %building the path of the text file
A = importdata(path); %import a text data as a structure
x=A(1);
for i=2:size(A,1)
x=strcat(x,'\n'); %preserving a line brake like this, it will be helpful while reconstruction
x=strcat(x, A(i));
end %this loop is to concatenate all the characters in a continuous sequence
x=x{1}; %converting structure to a character array
c=dec2bin(x,8); %create a cell of x character & 8bit
y=c(1,:); %for a continoues bit stream
for i=2:size(c,1)
y=strcat(y,c(i,:));
end
for i=1:size(y,2)
z(i)=double(y(i));
end
z(find(z(:)==48))=0;
z(find(z(:)==49))=1;%z is a bitstream now
disp('The Resultant Bit Stream is:');
disp(z);

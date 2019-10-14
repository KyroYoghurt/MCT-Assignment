function [x] = ConvolutionDecode(y, tblen)
K=3;
G1=5;
G2=7;
trel=poly2trellis(K,[G1 G2]);
x=vitdec(y,trel,tblen,'trunc','hard')
%tblen = length(x);
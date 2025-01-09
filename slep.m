function [matice_vystup] = slep(matice_vstup1,matice_vstup2)
%SLEP Summary of this function goes here
%   Detailed explanation goes here
[n,m]=size(matice_vstup1);
[k,l]=size(matice_vstup2);
if n<k
    matice_vstup1=[matice_vstup1;zeros(k-n,m)];
elseif k<n
    matice_vstup2=[matice_vstup2;zeros(n-k,l)];
end
matice_vystup=[matice_vstup1,matice_vstup2];
end


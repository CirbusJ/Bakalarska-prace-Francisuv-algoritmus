function [vystupni_matice] = vytvor_testovaci_matici(rozmer)
%VYTVOR_TESTOVACI_MATICI Summary of this function goes here
%   Detailed explanation goes here
diagonala=zeros(1,rozmer);
for j=1:rozmer
    diagonala(j)=j+j*1i;
end
S=rand(rozmer);
vystupni_matice=inv(S)*diag(diagonala)*S;
end


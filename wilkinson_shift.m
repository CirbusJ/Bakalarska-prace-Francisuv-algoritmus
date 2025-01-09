function [spektrum] = wilkinson_shift(vstupni_matice,tolerance_deflace,max_pocet_iteraci,spektrum)
%WILKINSON_SHIFT Najde vlastni cisla matice maleho rozmeru Francisovym algoritmem s jednonasobnym shiftem
if nargin<=1
    spektrum=[];
    tolerance_deflace=10^-6;
    max_pocet_iteraci=10000;
end

n=size(vstupni_matice,1);
if n==1
    spektrum=[spektrum,vstupni_matice];
    return
end
if n==2
    a=vstupni_matice(1,1);
    b=vstupni_matice(1,2);
    c=vstupni_matice(2,1);
    d=vstupni_matice(2,2);
    x1=(a+d+sqrt((a+d)^2-4*(a*d-(c*b))))/2;
    x2=(a+d-sqrt((a+d)^2-4*(a*d-(c*b))))/2;
    spektrum=[spektrum,x2,x1];
    return
end

for i=1:max_pocet_iteraci
    for j=1:n-1
        if abs(vstupni_matice(j+1,j))<=tolerance_deflace
            %disp(i)
            [spektrum]=wilkinson_shift(vstupni_matice(1:j,1:j),tolerance_deflace,max_pocet_iteraci,spektrum);
            [spektrum]=wilkinson_shift(vstupni_matice(j+1:n,j+1:n),tolerance_deflace,max_pocet_iteraci,spektrum);
            return
        end
    end
    diagonala=diag(vstupni_matice);
    shift=diagonala(n);
    vstupni_matice=create_bulge(vstupni_matice,shift,1);
    vstupni_matice=bulge_chasing(vstupni_matice,1);
    if i==max_pocet_iteraci
        disp("dosli iterace ve Wilkinsonove shiftu")
    end
end
end


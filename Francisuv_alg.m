function [chyba_vl_cisel,spektrum,celkovy_pocet_iteraci] = Francisuv_alg(matice_vstup,velikost_shiftu,max_pocet_iteraci,tolerance_deflace,chyba_vl_cisel,spektrum)
%FRANCISUV_ALG Najde vsechna vlastni cisla obecne matice

%VSTUP: matice_vstup-matice, ktere chceme zjistit vl. cisla
%       velikost_shiftu-nasobnost shiftu
%       uzivatel zadava pouze dva parametry na vstupu a to matice_vstup a
%       velikost_shiftu
%       dalsi vstupy jsou pro rekurzivni volani funkce

%VYSTUP:chyba_vl_cisel-vraci prazdny vektor, artefakt z druhe implementace
%       spektrum-vsechna vlastni cisla vstupni matice
%       celkovy_pocet_iteraci-pocet provedenych vytvoreny bulge a
%       nasledneho bulge chasingu potrebnych k nalezeni celeho spektra

%prevede obecnou ctvercovou matici na horni Hessenberguv tvar
%a pote pomoci funkci create_bulge a bulge_chasing po nekolika
%iteracich aproximuje jedno z vlastnich cisel vstupni matice
%pokud nejaky z poddiagonalnich prvku je mensi nez pozadovana tolerance,
%pak algoritmus provede deflaci implementovanou rekurzi a tim redukuje
%problem na 2 mensi, timto postupem algoritmus najde vsechna vlastni cisla

H=matice_vstup;
global celkovy_pocet_iteraci;
if nargin<=2                    %tento if blok se provede pri prvnim zavolani funkce
    %velikost_shiftu=2;
    H=hess(H);
    max_pocet_iteraci=10000;      %zde je mozne zmenit parametry algoritmu
    tolerance_deflace=10^-6;
    celkovy_pocet_iteraci=0;
    chyba_vl_cisel=[];          %inicializace vystupnich argumentu
    spektrum=[];
end
n=size(H,1);
if n==1                         %pokud je matice 1x1 pridej toto cislo do spektra
    spektrum=[spektrum,H];
    
    return
elseif n==2                     %pokud je matice 2x2 najdi vlastni cisla teto matice explicitnim vzorcem
    a=H(1,1);
    b=H(1,2);
    c=H(2,1);
    d=H(2,2);
    x1=(a+d+sqrt((a+d)^2-4*(a*d-(c*b))))/2;
    x2=(a+d-sqrt((a+d)^2-4*(a*d-(c*b))))/2;
    spektrum=[spektrum,x1,x2];
    return
else                             
diagonala=diag(H);

if n<=6                                     
    velikost_shiftu=1;
end

shift=diagonala(n-velikost_shiftu+1:n);     %vybere prvni shift

for i=1:max_pocet_iteraci       %proces nulovani jednoho z poddiagonalnich prvku
    for j=1:n-1 
        if abs(H(j+1,j))<=tolerance_deflace %kontroluje zda uz je nektery z poddiagonalnich prvku dostanecne maly
            celkovy_pocet_iteraci=celkovy_pocet_iteraci+i-1;
            [chyba_vl_cisel,spektrum]=Francisuv_alg(H(1:j,1:j),velikost_shiftu,max_pocet_iteraci,tolerance_deflace,chyba_vl_cisel,spektrum);      %rekurzivni zavolani algoritmu na matice nizsich rozmeru
            [chyba_vl_cisel,spektrum]=Francisuv_alg(H(j+1:n,j+1:n),velikost_shiftu,max_pocet_iteraci,tolerance_deflace,chyba_vl_cisel,spektrum);
            return
        end
    end
    H=create_bulge(H,shift,velikost_shiftu);           %vytvori hrbol na pozadovanem miste
    H=bulge_chasing(H,velikost_shiftu);                %vyzene hrbol z matice a tim ji vrati do horniho Hessenbergova tvaru

    diagonala=diag(H);
    shift=diagonala(n-velikost_shiftu+1:n);

    if i==max_pocet_iteraci
        disp("dosli iterace :(")
    end
end
end
end


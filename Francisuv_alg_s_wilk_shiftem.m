function [konvergence_vl_cisel,chyba_vl_cisel,spektrum,celkovy_pocet_iteraci] = Francisuv_alg_s_wilk_shiftem(matice_vstup,velikost_shiftu,max_pocet_iteraci,tolerance_deflace,konvergence_vl_cisel,chyba_vl_cisel,spektrum)
%FRANCISUV_ALG Najde vsechna vlastni cisla obecne matice

%VSTUP: matice_vstup-matice, ktere chceme zjistit vl. cisla
%       uzivatel zadava pouze jeden vstup a to matice_vstup
%       dalsi vstupy jsou pro rekurzivni volani funkce
%VYSTUP:konvergence_vl_cisel-ve sloupcich jsou shifty, ktere se postupne
%priblizuji nejakemu vlastnimu cislu, radky odpovidaji indexu iterace
%       chyba_vl_cisel-ve sloupcich jsou poddiagonalni prvky a radky
%odpovidaji indexu iterece
%       spektrum-vsechna vlastni cisla vstupni matice


%prevede obecnou ctvercovou matici na horni Hessenberguv tvar
%a pote pomoci funkci create_bulge a bulge_chasing po pozadovem poctu
%iteraci aproximuje jedno z vlastnich cisel vstupni matice
%pokud nejaky z poddiagonalnich prvku je mensi nez pozadovana tolerance,
%pak algoritmus provede deflaci implementovanou rekurzi a tim redukuje
%problem na 2 mensi, timto postupem algoritmus najde vsechna vlastni cisla

H=matice_vstup;
global celkovy_pocet_iteraci;
if nargin<=2                    %tento if blok se provede pri prvnim zavolani funkce
    %velikost_shiftu=2;
    H=hess(H);
    max_pocet_iteraci=100000;      %zde je mozne zmenit parametry algoritmu
    tolerance_deflace=10^-6;    
    celkovy_pocet_iteraci=0;
    %chyba_vl_cisel=[];          %inicializace vystupnich argumentu
    spektrum=[];
    konvergence_vl_cisel=[];
end
n=size(H,1);
if n==1                         %pokud je matice 1x1 pridej toto cislo do spektra
    spektrum=[spektrum,H];
    aproximace_vl_cisla=H;
    %konvergence_vl_cisel=slep(konvergence_vl_cisel,vlastni_cislo);
    return
elseif n==2                     %pokud je matice 2x2 najdi vlastni cisla teto matice explicitnim vzorcem
    a=H(1,1);
    b=H(1,2);
    c=H(2,1);
    d=H(2,2);
    x1=(a+d+sqrt((a+d)^2-4*(a*d-(c*b))))/2;
    x2=(a+d-sqrt((a+d)^2-4*(a*d-(c*b))))/2;
    spektrum=[spektrum,x1,x2];
    aproximace_vl_cisla=[x1,x2];
    konvergence_vl_cisel=slep(konvergence_vl_cisel,aproximace_vl_cisla);
    %chyba_vl_cisel=slep(chyba_vl_cisel,[0,0]);
    return
else                             
diagonala=transpose(diag(H));
poddiagonala=abs(transpose(diag(H,-1)));

if n<=6
    velikost_shiftu=1;
end

if velikost_shiftu==1
    %aproximace_vl_cisla=[kandidati_na_shift(2)]; %poznamka pro Jana, zde muzes zlepsit volbu shiftu
    aproximace_vl_cisla=serazeni_kandidatu(wilkinson_shift(H(n-1:n,n-1:n)),diagonala(n));
    if celkovy_pocet_iteraci==0
    chyba=[poddiagonala(n-5:n-1)];
    end  
else
    %aproximace_vl_cisla=[aproximace_vl_cisla;diagonala(n-velikost_shiftu+1:n)];
    aproximace_vl_cisla=serazeni_kandidatu(wilkinson_shift(H(n-velikost_shiftu+1:n,n-velikost_shiftu+1:n)),diagonala(n-velikost_shiftu+1:n));
    if celkovy_pocet_iteraci==0
    chyba=[poddiagonala(n-5:n-1)];
    end
end
priklad_otravne_matice=H;
for i=1:max_pocet_iteraci       %proces nulovani jednoho z poddiagonalnich prvku
    for j=1:n-1 
        if abs(H(j+1,j))<=tolerance_deflace %kontroluje zda uz je nektery z poddiagonalnich prvku dostanecne maly
            %disp(i-1)
            if celkovy_pocet_iteraci==0
            chyba_vl_cisel=chyba;
            end
            celkovy_pocet_iteraci=celkovy_pocet_iteraci+i-1;
            konvergence_vl_cisel=slep(konvergence_vl_cisel,aproximace_vl_cisla);
            %chyba_vl_cisel=slep(chyba_vl_cisel,chyba);
            [konvergence_vl_cisel,chyba_vl_cisel,spektrum]=Francisuv_alg_s_wilk_shiftem(H(1:j,1:j),velikost_shiftu,max_pocet_iteraci,tolerance_deflace,konvergence_vl_cisel,chyba_vl_cisel,spektrum);      %rekurzivni zavolani algoritmu na matice nizsich rozmeru
            [konvergence_vl_cisel,chyba_vl_cisel,spektrum]=Francisuv_alg_s_wilk_shiftem(H(j+1:n,j+1:n),velikost_shiftu,max_pocet_iteraci,tolerance_deflace,konvergence_vl_cisel,chyba_vl_cisel,spektrum);
            return
        end
    end
    H=create_bulge(H,aproximace_vl_cisla(i,:),velikost_shiftu);           %vytvori hrbol na pozadovanem miste
    H=bulge_chasing(H,velikost_shiftu);                                 %vyzene hrbol z matice a tim ji vrati do horniho Hessenbergova tvaru
    
    diagonala=transpose(diag(H));
    poddiagonala=abs(transpose(diag(H,-1)));

    if velikost_shiftu==1
    %kandidati_na_shift=wilkinson_shift(H(n-1:n,n-1:n));
    aproximace_vl_cisla=[aproximace_vl_cisla;serazeni_kandidatu(wilkinson_shift(H(n-1:n,n-1:n)),diagonala(n))]; %poznamka pro Jana, zde muzes zlepsit volbu shiftu
    if celkovy_pocet_iteraci==0
    chyba=[chyba;poddiagonala(n-5:n-1)];
    end
    else
    %aproximace_vl_cisla=[aproximace_vl_cisla;diagonala(n-velikost_shiftu+1:n)];
    aproximace_vl_cisla=[aproximace_vl_cisla;serazeni_kandidatu(wilkinson_shift(H(n-velikost_shiftu+1:n,n-velikost_shiftu+1:n)),diagonala(n-velikost_shiftu+1:n))];
    if celkovy_pocet_iteraci==0
    chyba=[chyba;poddiagonala(n-5:n-1)];
    end
    end

    if i==max_pocet_iteraci
        disp("dosli iterace :(")
        disp(aproximace_vl_cisla(i,:))
        %disp(poddiagonala)
        %disp(H)
        %size(H)
        %disp(diagonala)
        %disp(priklad_otravne_matice)
    end
end
end
end


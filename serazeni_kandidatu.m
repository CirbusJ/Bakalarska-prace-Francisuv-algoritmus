function [spravne_serazeny_vektor] = serazeni_kandidatu(vektor_na_serazeni,vzor)
%SERAZENI_KANDIDATU Seradi kandidaty na shifty do spravneho poradi
n=size(vektor_na_serazeni,2);
m=size(vzor,2);
spravne_serazeny_vektor=zeros(1,m);
if m==1
    if abs(vzor-vektor_na_serazeni(1))<abs(vzor-vektor_na_serazeni(2))
        spravne_serazeny_vektor(1)=vektor_na_serazeni(1);
    else
        spravne_serazeny_vektor(1)=vektor_na_serazeni(2);
    end
else
    for i=1:n
        minimum=inf;
        index=0;
        for j=1:n
            if abs(vzor(i)-vektor_na_serazeni(j))<minimum
                index=j;
                minimum=abs(vzor(i)-vektor_na_serazeni(j));
            end
        end
        spravne_serazeny_vektor(i)=vektor_na_serazeni(index);
    end
end
end


function [vystupni_matice] = podobnosti_transformace_core_transformaci(vstupni_matice,core_tranformation,index_core_tranformation)
%PODOBNOSTI_TRANSFORMACE_CORE_TRANSFORMACI provede efektivne operaci B=Q^*AQ
velikost_matice=size(vstupni_matice,1);
i=index_core_tranformation;
vystupni_matice=vstupni_matice;

if i==size(core_tranformation,1)
    return
end
c=core_tranformation(i,i);
s=core_tranformation(i+1,i);

%nejprve operace AQ
for k=1:velikost_matice
    prvek_ki=vystupni_matice(k,i);
    prvek_ki1=vystupni_matice(k,i+1);
    %if prvek_ki==0 && prvek_ki1==0
    %    vystupni_matice(k,i)=0;
    %    vystupni_matice(k,i+1)=0;
    %else
        vystupni_matice(k,i)=prvek_ki*c+prvek_ki1*s;
        vystupni_matice(k,i+1)=prvek_ki*(-conj(s))+prvek_ki1*conj(c);
    %end
end
%nyni operace Q^*(AQ)
for l=1:velikost_matice
    prvek_il=vystupni_matice(i,l);
    prvek_i1l=vystupni_matice(i+1,l);
    %if prvek_il==0 && prvek_i1l==0
    %    vystupni_matice(i,l)=0;
    %    vystupni_matice(i+1,l)=0;
    %else
        vystupni_matice(i,l)=prvek_il*conj(c)+prvek_i1l*conj(s);
        vystupni_matice(i+1,l)=prvek_il*(-s)+prvek_i1l*c;
    %end
end
end
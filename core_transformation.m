function [vystupni_matice] = core_transformation(radek_kde_operuje,sloupec_kde_operuje,vstupni_matice)
%CORE_TRANSFORMATION vytvori matici, ktera reprezentuje core
%transformation
i=radek_kde_operuje;
j=sloupec_kde_operuje;
if(i+1>size(vstupni_matice,1)||j>size(vstupni_matice,1))
    vystupni_matice=eye(size(vstupni_matice,1));
    return
end
r=norm(vstupni_matice(i:i+1,j));
if r~=0
c=vstupni_matice(i,j)/r;
s=vstupni_matice(i+1,j)/r;
else
c=1;
s=0;
end
Q=eye(size(vstupni_matice,1));
Q(i,i)=c;
Q(i+1,i+1)=conj(c);
Q(i,i+1)=-conj(s);
Q(i+1,i)=s;
vystupni_matice=Q;
end


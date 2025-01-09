function [vystupni_matice] = bulge_chasing(vstupni_matice,velikost_shiftu)
%BULGE_CHASING vstupni matici, ktera obsahuje bulge, prevede pomoci core
%transformations zpet na horni Hessenberguv tvar
vystupni_matice=vstupni_matice;
n=size(vstupni_matice,1);
for i=2:n-1
    if i+velikost_shiftu-1>n
    for j=n-1:-1:i
    Q=core_transformation(j,i-1,vystupni_matice);
    %vystupni_matice=Q'*vystupni_matice*Q;
    vystupni_matice=podobnosti_transformace_core_transformaci(vystupni_matice,Q,j);
    end
    else
    for j=i+velikost_shiftu-1:-1:i
    Q=core_transformation(j,i-1,vystupni_matice);
    %vystupni_matice=Q'*vystupni_matice*Q;
    vystupni_matice=podobnosti_transformace_core_transformaci(vystupni_matice,Q,j);
    end
    end
    %spy(vystupni_matice)
    %pause
    %disp(vystupni_matice)
end
end


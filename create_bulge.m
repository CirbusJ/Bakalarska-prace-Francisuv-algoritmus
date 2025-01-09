function [vystupni_matice] = create_bulge(vstupni_matice,shift,velikost_shiftu)
%CREATE_BULGE do vstupni matice vytvori bulge pomoci core transformation s
%pozadovanym shiftem
for i=1:velikost_shiftu
vstupni_matice_shift(:,:,i)=vstupni_matice-shift(i)*eye(size(vstupni_matice,1));
end
%vytvor vektor x=f(A)*e_1
x=eye(size(vstupni_matice,1),1);
for i=1:velikost_shiftu
x=vstupni_matice_shift(:,:,i)*x;
end
vystupni_matice=vstupni_matice;
for i=1:velikost_shiftu
    Q=core_transformation(velikost_shiftu-i+1,1,x);
    x=Q'*x;
    %vystupni_matice=Q'*vystupni_matice*Q;
    vystupni_matice=podobnosti_transformace_core_transformaci(vystupni_matice,Q,velikost_shiftu-i+1);
end
%spy(vystupni_matice)
%pause
%disp(vystupni_matice)
%disp("matice nade mnou je vytvoreny bulge")
end
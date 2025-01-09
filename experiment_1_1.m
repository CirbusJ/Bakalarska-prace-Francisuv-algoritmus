%Q=rand(100);
%testovaci_matice=inv(Q)*diag([2:2:200])*Q;

testovaci_matice=vytvor_testovaci_matici(100);

[~,~,spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg(testovaci_matice,1);

vlastni_cisla=1:100;
for j=1:100
    vlastni_cisla(j)=(1+1i)*j;
end


spravne_serazene_aproximace=serazeni_kandidatu(spektrum_ja,vlastni_cisla);
disp("Maximová norma rozdílu")
disp(norm(vlastni_cisla-spravne_serazene_aproximace,"inf"))


f1=figure;

%set(gcf,'Position',[0 0 600 600])
plot(real(vlastni_cisla),imag(vlastni_cisla),"o");
hold on
plot(real(spektrum_ja),imag(spektrum_ja),"*");
FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('přesná vlastní čísla','aproximace vlastních čísel');
leg.Location='southeast';
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,FS,fs);
set(leg,FS,fs);
print('-dpng','test.png')
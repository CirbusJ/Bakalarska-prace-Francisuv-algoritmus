%testovaci_matice=vytvor_testovaci_matici(100);

Q=rand(100);
testovaci_matice=inv(Q)*diag([1:1:100])*Q;
%priklad = matfile('nahodna_matice_exp_2_1.mat');
%nahodna_matice_exp_2_1=priklad.Q;
%testovaci_matice=inv(nahodna_matice_exp_2_1)*diag([1:1:100])*nahodna_matice_exp_2_1;
tic
[~,~,spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg(testovaci_matice,2);
toc
tic
[~,~,spektrum_ja_wilk,celkovy_pocet_iteraci_wilk]=Francisuv_alg_s_wilk_shiftem(testovaci_matice,2);
toc

vlastni_cisla=1:100;
for j=1:100
    vlastni_cisla(j)=j;
end


spravne_serazene_aproximace=serazeni_kandidatu(spektrum_ja,vlastni_cisla);
disp("Maximová norma rozdílu Rayleigh")
disp(norm(vlastni_cisla-spravne_serazene_aproximace,"inf"))

spravne_serazene_aproximace_wilk=serazeni_kandidatu(spektrum_ja_wilk,vlastni_cisla);
disp("Maximová norma rozdílu Wilkinson")
disp(norm(vlastni_cisla-spravne_serazene_aproximace_wilk,"inf"))

disp("Celkový počet iterací Rayleigh")
disp(celkovy_pocet_iteraci)

disp("Celkový počet iterací Wilkinson")
disp(celkovy_pocet_iteraci_wilk)

f1=figure;


plot(real(vlastni_cisla),imag(vlastni_cisla),"o");
hold on
plot(real(spektrum_ja),imag(spektrum_ja),"*");
p=plot(real(spektrum_ja_wilk),imag(spektrum_ja_wilk),"+");
p.Color="#77AC30";
FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('přesná vlastní čísla','aproximace vlastních čísel-Rayleigh','aproximace vlastních čísel-Wilkinson');
leg.Location='southeast';
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,FS,fs);
set(leg,FS,fs);
print('-dpng','test.png')
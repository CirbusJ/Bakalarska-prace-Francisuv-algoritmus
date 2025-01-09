testovaci_matice=randn(100);


tic
[~,~,spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg(testovaci_matice,2);
toc
tic
[~,~,spektrum_ja_wilk,celkovy_pocet_iteraci_wilk]=Francisuv_alg_s_wilk_shiftem(testovaci_matice,2);
toc

spektrum_eig=eig(testovaci_matice);

spravne_serazene_aproximace=serazeni_kandidatu(spektrum_ja,transpose(spektrum_eig));
disp("Maximová norma rozdílu Rayleigh")
disp(norm(transpose(spektrum_eig)-spravne_serazene_aproximace,"inf"))

spravne_serazene_aproximace_wilk=serazeni_kandidatu(spektrum_ja_wilk,transpose(spektrum_eig));
disp("Maximová norma rozdílu Wilkinson")
disp(norm(transpose(spektrum_eig)-spravne_serazene_aproximace_wilk,"inf"))

disp("Celkový počet iterací Rayleigh")
disp(celkovy_pocet_iteraci)

disp("Celkový počet iterací Wilkinson")
disp(celkovy_pocet_iteraci_wilk)

disp("Počet nazelených vlastních čísel Rayleigh")
disp(size(spektrum_ja,2))

f1=figure;


plot(real(spektrum_eig),imag(spektrum_eig),"o");
hold on
plot(real(spektrum_ja),imag(spektrum_ja),"*");
p=plot(real(spektrum_ja_wilk),imag(spektrum_ja_wilk),"+");
p.Color="#77AC30";
FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('aproximace vlastních čísel-eig','aproximace vlastních čísel-Rayleigh','aproximace vlastních čísel-Wilkinson');
leg.Location='northeast';
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,FS,fs);
set(leg,FS,fs);
print('-dpng','test.png')
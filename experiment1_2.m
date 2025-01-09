%Experiment 1 příklad 2
%Francisuv algoritmus s jedním shiftem hledá vlastní čísla porušené a neporušené matice s
%citlivými vlastními čísly a porovnává je s aproximacemi z vestavěné funkce eig


%Výběr testovací matice
%testovaci_matice=gallery('lesp',50);
testovaci_matice=gallery('grcar',50);

E=rand(50)*10^(-6);
porusena_testovaci_matice=testovaci_matice+E;

[~,spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg(testovaci_matice,1);
[~,porusene_spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg(porusena_testovaci_matice,1);

spektrum_eig=eig(testovaci_matice);
porusene_spektrum_eig=eig(porusena_testovaci_matice);

spravne_serazene_aproximace=serazeni_kandidatu(spektrum_ja,transpose(spektrum_eig));

disp("Maximová norma rozdílu")
disp(norm(transpose(spektrum_eig)-spravne_serazene_aproximace,"inf"))

disp("Celkový počet iterací")
disp(celkovy_pocet_iteraci)


f1=figure;
plot(real(spektrum_eig),imag(spektrum_eig),"o");
hold on
plot(real(spektrum_ja),imag(spektrum_ja),"*");

FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('neporušená matice-eig','neporušená matice-Francis');
leg.Location='northeast';
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,FS,fs);
set(leg,FS,fs);

f2=figure;
plot(real(porusene_spektrum_eig),imag(porusene_spektrum_eig),"o");
hold on
plot(real(porusene_spektrum_ja),imag(porusene_spektrum_ja),"*");
FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('porušená matice-eig','porušená matice-Francis');
leg.Location='northeast';
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,FS,fs);
set(leg,FS,fs);

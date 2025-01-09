
testovaci_matice=gallery('lesp',50);

E=rand(50)*10^(-6);
porusena_testovaci_matice=testovaci_matice+E;

[~,~,spektrum_ja,~]=Francisuv_alg(testovaci_matice,1);
[~,~,porusene_spektrum_ja,~]=Francisuv_alg(porusena_testovaci_matice,1);

spektrum_eig=eig(testovaci_matice);
porusene_spektrum_eig=eig(porusena_testovaci_matice);

spravne_serazene_aproximace=serazeni_kandidatu(spektrum_ja,transpose(spektrum_eig));
disp("Maximová norma rozdílu")
disp(norm(transpose(spektrum_eig)-spravne_serazene_aproximace,"inf"))


f1=figure;

plot(real(spektrum_eig),imag(spektrum_eig),"o");
hold on
plot(real(spektrum_ja),imag(spektrum_ja),"*");
plot(real(porusene_spektrum_eig),imag(porusene_spektrum_eig),"diamond");
plot(real(porusene_spektrum_ja),imag(porusene_spektrum_ja),"+");


FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('neporušená matice-eig','neporušená matice-Francis','porušená matice-eig','porušená matice-Francis');
leg.Location='northeast';
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,FS,fs);
set(leg,FS,fs);
print('-dpng','test2.png')
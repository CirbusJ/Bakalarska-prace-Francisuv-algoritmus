%Experiment 3 příklad 1
%Zde testujeme konvergenci posledních pěti poddiagonálních prvků před první
%deflací ve vztahu k velikosti shiftu, kterou musí uživatel nastovovat ručně
%Testujeme implementaci s Wilkinsonovým shiftem na matici s vlastními čísly 
%rozloženými rozvnoměrně na intervalu [1,100]
%Měříme zároveň výpočetní čas

%!!!Pozor, aby uživatel sledoval stejnou matici pro různé velikostí shiftů
%musí po prvním spuštění zakomentovat náhodnou matici S !!!
S=rand(100);
testovaci_matice=inv(S)*diag([1:1:100])*S;

disp('Výpočetní čas algoritmu:')
tic
[konvergence_poddiagonich_prvku,spektrum_ja_wilk,celkovy_pocet_iteraci_wilk]=Francisuv_alg_s_wilk_shiftem(testovaci_matice,1); %<- Zde přenastavit velikost shiftu
toc

spektrum_eig=eig(testovaci_matice);

pocet_iteraci_v_prvnim_kroku=size(konvergence_poddiagonich_prvku,1);

spravne_serazene_aproximace_wilk=serazeni_kandidatu(spektrum_ja_wilk,transpose(spektrum_eig));
disp("Maximová norma rozdílu Wilkinson")
disp(norm(transpose(spektrum_eig)-spravne_serazene_aproximace_wilk,"inf"))


disp("Celkový počet iterací Wilkinson")
disp(celkovy_pocet_iteraci_wilk)

f1=figure;
plot(real(spektrum_eig),imag(spektrum_eig),"o");
hold on
p=plot(real(spektrum_ja_wilk),imag(spektrum_ja_wilk),"*");

f2=figure;

semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,1),"red",LineWidth=1.4)
hold on
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,2),Color="#77AC30",LineWidth=1.4)
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,3),"blue",LineWidth=1.4)
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,4),Color="#4DBEEE",LineWidth=1.4)
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,5),Color="#7E2F8E",LineWidth=1.4)
p.Color="#77AC30";
FS='FontSize';fs=12;
title=('\fontsize{15}Vlastní čísla');
leg=legend('$a_{n-4,n-5}$','$a_{n-3,n-4}$','$a_{n-2,n-3}$','$a_{n-1,n-2}$','$a_{n,n-1}$');
leg.Location='southwest';
xlabel('Index kroku');
ylabel('Velikost poddiagonálních prvků');
a = axis;
set(gca, 'XTick', 1:pocet_iteraci_v_prvnim_kroku)
set(leg,'Interpreter','latex')
set(gca,FS,fs);
set(leg,FS,fs);

S=rand(100);
testovaci_matice=inv(S)*diag([1:1:100])*S;
%testovaci_matice=gallery('wilk',21);
%porusena_testovaci_matice=testovaci_matice+10^(-6)*rand(50);
%testovaci_matice=randi([-1 1],200);
%testovaci_matice=vytvor_testovaci_matici(200);
%testovaci_matice=[0.8483   -1.0720    0.1006   -0.3596   -0.4892;
%    0.6264    0.9592    0.2961    0.5414   -0.0388;
%   -0.0000    0.0294    1.3111    0.4267   -0.8964;
%   -0.0000    0.0000    0.0444    0.3442   -1.1004;
%   -0.0000   -0.0000         0    0.9334    0.3200];
%testovaci_matice=magic(200);
%testovaci_matice=[  227.9546    5.3751   -0.4616   13.2422;
%                    1.2017   79.9471   -9.5897    0.2507;
%                    0.0000    1.9581  227.5181   -5.5492;
%                    0.0000         0   10.2580   80.5803];
%testovaci_matice=[-2.0471   -2.0884   -0.2163    0.0410;
%                2.4483   -2.1988    0.1005    0.0316;
%                0.0000    0.0029   -0.2541   -2.7750;
%                0.0000         0    3.1012   -0.2514];
cond(testovaci_matice)
tic
[~,konvergence_poddiagonich_prvku,spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg_s_wilk_shiftem(testovaci_matice,2);
%[~,~,spektrum_ja_2,celkovy_pocet_iteraci_2]=Francisuv_alg_s_wilk_shiftem(porusena_testovaci_matice,2);
%[~,~,spektrum_ja,celkovy_pocet_iteraci]=Francisuv_alg(testovaci_matice,1);
toc
tic
spektrum_eig=eig(testovaci_matice);
toc

pocet_iteraci_v_prvnim_kroku=size(konvergence_poddiagonich_prvku,1);

f1=figure;
%set(gcf,'Position',[0 0 600 600])
plot(real(spektrum_eig),imag(spektrum_eig),"o");
hold on
plot(real(spektrum_ja),imag(spektrum_ja),"*");
%plot(real(spektrum_ja_2),imag(spektrum_ja_2),"+");

f2=figure;
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,1),"red")
hold on
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,2),"green")
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,3),"blue")
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,4),"cyan")
semilogy(1:pocet_iteraci_v_prvnim_kroku,konvergence_poddiagonich_prvku(1:pocet_iteraci_v_prvnim_kroku,5),"magenta")

size(spektrum_eig)
size(spektrum_ja)
disp(celkovy_pocet_iteraci)

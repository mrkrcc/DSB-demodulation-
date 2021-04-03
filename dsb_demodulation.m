%% A Zaman Domeni Analizleri
clear all, close all; clc 
[x , Fs] = audioread('davul_gitar.wav'); 
Ts = 1/Fs; 
t = [0:numel(x)-1]*Ts; 
figure,plot(t,x),grid on 
xlabel('Zaman [sn]'), ylabel('Genlik'); 
title(['DSB Modüleli Ses Isareti Fs:', num2str(Fs)]) 
sound(x,Fs) 
%% B Frekans Domeni Analizleri
X = fftshift(abs(fft(x))); 
F = linspace(-Fs/2 , Fs/2 , numel(X)); 
figure,plot(F,X),grid on,xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]'), ylabel('Genlik [V]')
title(['Frekans Analizi Ýle Taþýyýcý Frekanslarýn Kesitirimi'])
% Figuru inceledigimizde 2 tane tasiyici frekans oldugunu gormekteyiz.
% 5e3 ve 15e3 
%% C Demodülasyon
Fc = 5e3; % davul icin secilmis tasýyýcý frekans
c = cos(2*pi*Fc*t)'; % tasiyici isaret  
s = x.*c; % demodulasyon icin module isareti lokal osilator ile carptik
%% C Demodülasyon
X = fftshift(abs(fft(x))); 
S = fftshift(abs(fft(s))); 
F = linspace(-Fs/2 , Fs/2 , numel(X)); 
figure, 
subplot(211),plot(F,X,'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaret'])
subplot(212),plot(F,S,'g-'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaretin Lokal OsilatOr Ile Carpilmasi'])
%% B frekans domeni analizleri
F_davul = 2500; % Fc >= 2 * F_davul olmasý gerektiði için F_davul = 2500(5e3 taþýyýcý frekans için) 
H = zeros(numel(F),1); 
for i = 1:numel(F) 
    if abs(F(i))<F_davul
        H(i) = 1; 
    end
end
figure,plot(F,X/max(X)),hold on,plot(F,H,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'),xlim([-Fs/2 , Fs/2])
title(['Filtre Ýçin F max:', num2str(F_davul)])
%% B frekans domeni analizleri
X = abs(fft(x)); 
Hs = fftshift(H); %filtreyi uyumlu hale getirmek için kaydýrýyoruz
figure,plot(X/max(X)),hold on, plot(Hs,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]')
%% C Demodülasyon
figure,plot(F,S/max(S),'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'), 
hold on,plot(F,H,'g-','linewidth',2) 
s_recover = ifft(Hs.*fft(s)); 
s_recover = real(s_recover); s_recover = s_recover/max(abs(s_recover));
davul = s_recover ;
sound(s_recover,Fs)
%% B frekans domeni analizleri
F_gitar = 7500; 
H = zeros(numel(F),1); 
for i = 1:numel(F) 
    if abs(F(i))<F_gitar
        H(i) = 1; 
    end
end
figure,plot(F,X/max(X)),hold on,plot(F,H,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'),xlim([-Fs/2 , Fs/2])
title(['Filtre Ýçin F max:', num2str(F_gitar)])
%% B frekans domeni analizleri
X = abs(fft(x)); 
Hs = fftshift(H);
figure,plot(X/max(X)),hold on, plot(Hs,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]')
%% C Demodülasyon
Fc = 15e3; % gitar icin secilmis tasýyýcý frekans
c = cos(2*pi*Fc*t)'; % tasiyici isaret  
s = x.*c; % demodulasyon icin module isareti lokal osilator ile carptik
%% C Demodülasyon
X = fftshift(abs(fft(x))); 
S = fftshift(abs(fft(s))); 
F = linspace(-Fs/2 , Fs/2 , numel(X)); 
figure,  
subplot(211),plot(F,X,'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaret'])
subplot(212),plot(F,S,'g-'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaretin Lokal Osilator Ile Carpilmasi'])
%% C Demodülasyon
figure,plot(F,S/max(S),'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'), 
hold on,plot(F,H,'g-','linewidth',2) 
s_recover = ifft(Hs.*fft(s)); 
s_recover = real(s_recover); s_recover = s_recover/max(abs(s_recover));
gitar = s_recover ;
davul_gitar = davul + gitar ;
%sound(davul,Fs)
sound(gitar,Fs)
%sound(davul_gitar,Fs) % clear sound
%% ...

...     0dB Gurultu 


%% A Zaman Domeni Analizleri
clear all, close all; clc 
[x , Fs] = audioread('davul_gitar.wav'); 
Ts = 1/Fs; 
t = [0:numel(x)-1]*Ts;
x = awgn(x,0,'measured'); %0dB gürültü eklendi

figure,plot(t,x),grid on 
xlabel('Zaman [sn]'), ylabel('Genlik'); 
title(['DSB Modüleli 0dB Gurultulu Ses Isareti Fs:', num2str(Fs)]) 
sound(x,Fs)
%% B Frekans Domeni Analizleri
X = fftshift(abs(fft(x))); 
F = linspace(-Fs/2 , Fs/2 , numel(X)); 
figure,plot(F,X),grid on,xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]'), ylabel('Genlik [V]')
title(['Frekans Analizi Ýle Taþýyýcý Frekanslarýn Kesitirimi'])
% Figuru inceledigimizde 2 tane tasiyici frekans oldugunu gormekteyiz.
% 5e3 ve 15e3 
%% C Demodülasyon
Fc = 5e3; % davul icin secilmis tasýyýcý frekans
c = cos(2*pi*Fc*t)'; % tasiyici isaret  
s = x.*c; % demodulasyon icin module isareti lokal osilator ile carptik
%% C Demodülasyon
X = fftshift(abs(fft(x))); 
S = fftshift(abs(fft(s))); 
F = linspace(-Fs/2 , Fs/2 , numel(X)); 
figure, 
subplot(211),plot(F,X,'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaret'])
subplot(212),plot(F,S,'g-'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaretin Lokal OsilatOr Ile Carpilmasi'])
%% B frekans domeni analizleri
F_davul = 2500; % Fc >= 2 * F_davul olmasý gerektiði için F_davul = 2500(5e3 taþýyýcý frekans için) 
H = zeros(numel(F),1); 
for i = 1:numel(F) 
    if abs(F(i))<F_davul
        H(i) = 1; 
    end
end
figure,plot(F,X/max(X)),hold on,plot(F,H,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'),xlim([-Fs/2 , Fs/2])
title(['Filtre Ýçin F max:', num2str(F_davul)])
%% B frekans domeni analizleri
X = abs(fft(x)); 
Hs = fftshift(H); %filtreyi uyumlu hale getirmek için kaydýrýyoruz
figure,plot(X/max(X)),hold on, plot(Hs,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]')
%%  C Demodülasyon
figure,plot(F,S/max(S),'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'), 
hold on,plot(F,H,'g-','linewidth',2) 
s_recover = ifft(Hs.*fft(s)); 
s_recover = real(s_recover); s_recover = s_recover/max(abs(s_recover));
davul = s_recover ;
sound(s_recover,Fs)
%% B frekans domeni analizleri
F_gitar = 7500; 
H = zeros(numel(F),1); 
for i = 1:numel(F) 
    if abs(F(i))<F_gitar
        H(i) = 1; 
    end
end
figure,plot(F,X/max(X)),hold on,plot(F,H,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'),xlim([-Fs/2 , Fs/2])
title(['Filtre Ýçin F max:', num2str(F_gitar)])
%% B frekans domeni analizleri
X = abs(fft(x)); 
Hs = fftshift(H);
figure,plot(X/max(X)),hold on, plot(Hs,'m','linewidth',2),grid on, 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]')
%% C Demodülasyon
Fc = 15e3; % gitar icin secilmis tasýyýcý frekans
c = cos(2*pi*Fc*t)'; % tasiyici isaret  
s = x.*c; % demodulasyon icin module isareti lokal osilator ile carptik
%% C Demodülasyon
X = fftshift(abs(fft(x))); 
S = fftshift(abs(fft(s))); 
F = linspace(-Fs/2 , Fs/2 , numel(X)); 
figure,  
subplot(211),plot(F,X,'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaret'])
subplot(212),plot(F,S,'g-'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]') 
title(['Modüleli Isaretin Lokal Osilator Ile Carpilmasi'])
%% C Demodülasyon
figure,plot(F,S/max(S),'r'),grid on,xlim([-Fs/2 , Fs/2]) 
xlabel('Frekans [Hz]'), ylabel('Genlik [V]'), 
hold on,plot(F,H,'g-','linewidth',2) 
s_recover = ifft(Hs.*fft(s)); 
s_recover = real(s_recover); s_recover = s_recover/max(abs(s_recover));
gitar = s_recover ;
davul_gitar = davul + gitar ;
%sound(davul,Fs)
sound(gitar,Fs)
%sound(davul_gitar,Fs) % clear sound
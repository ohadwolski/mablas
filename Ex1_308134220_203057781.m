% Or Avner      203057781
% Ohad Wolski   308134220 55555

clear all; close all; clc;

F_c = 16000;


% Question no. 1
w0 = 2*pi*250;
w = linspace(-300*2*pi,300*2*pi,F_c);
X_F_w0 = (1/(2j))*(1./(2+1j*(w-w0)))-(1/(2j))*(1./(2+1j*(w+w0)));
figure(1);
subplot(3,1,1);
plot(w, abs(X_F_w0));
hold on;

% Question no. 2
w1 = 2*pi*250;
w2 = 2*pi*315;
w3 = 2*pi*375;
w4 = 2*pi*500;
t = linspace(0,4,4*F_c);
t_sec = linspace(0,1,F_c);
x_wn = @(t,wn) sin(wn.*t).*heaviside(t).*exp(-2.*t);
x_t = [arrayfun(@(t) x_wn(t,w1), t_sec) arrayfun(@(t) x_wn(t,w2), t_sec) arrayfun(@(t) x_wn(t,w3), t_sec) arrayfun(@(t) x_wn(t,w4), t_sec)];

soundsc(x_t,F_c);

% Question no. 3
index1 = 2.99*F_c;
index2 = 3.01*F_c;
t_cut = t(index1:index2);
x_t_cut = x_t(index1:index2);

figure(2);
subplot(2,1,1);
plot(t_cut, x_t_cut, ':');
xlim([2.99 3.01]);
ylim([-1 1]);
hold on;

subplot(2,1,2);
plot(t_cut, x_t_cut, ':');
xlim([2.99 3.01]);
ylim([-1 1]);
hold on;

% Question no. 4
Fs1 = 2000;
t_n1 = downsample(t, F_c/Fs1);
x_n1 = downsample(x_t, F_c/Fs1);
subplot(2,1,1);
stem(t_n1(2.99*Fs1:3.01*Fs1),x_n1(2.99*Fs1:3.01*Fs1));

Fs2 = 800;
t_n2 = downsample(t, F_c/Fs2);
x_n2 = downsample(x_t, F_c/Fs2);
subplot(2,1,2);
stem(t_n2(2.99*Fs2:3.01*Fs2),x_n2(2.99*Fs2:3.01*Fs2));


% Question no. 5
figure(1);
subplot(3,1,2);
X_k1 = fftshift(fft(x_n1));
w_k1 = linspace(-Fs1/2,Fs1/2,8000);
plot(w_k1, abs(X_k1));
xlim([-800 800]);

subplot(3,1,3);
X_k2 = fftshift(fft(x_n2));
w_k2 = linspace(-Fs2/2,Fs2/2,3200);
plot(w_k2, abs(X_k2));
xlim([-800 800]);


% Question no. 7
figure(2);
subplot(2, 1, 1);
x_comb1 = upsample(x_n1, F_c/Fs1);
t_sinc = linspace(-255/F_c,256/F_c,512);
sinc_rec1 = sinc((Fs1 * t_sinc));
x_ideal1 = conv(x_comb1, sinc_rec1, 'same');
plot(t_cut,x_ideal1(index1:index2));

subplot(2, 1, 2);
x_comb2 = upsample(x_n2, F_c/Fs2);
sinc_rec2 = sinc((Fs2 * t_sinc));
x_ideal2 = conv(x_comb2, sinc_rec2, 'same');
plot(t_cut,x_ideal2(index1:index2));


% Question no. 8
subplot(2, 1, 1);
zoh1 = [zeros(1,(F_c/Fs1)/2) ones(1,F_c/Fs1) zeros(1,(F_c/Fs1)/2)]; 
x_zoh1 = conv(x_comb1, zoh1, 'same');
plot(t_cut,x_zoh1(index1:index2));

subplot(2, 1, 2);
zoh2 = [zeros(1,(F_c/Fs2)/2) ones(1,F_c/Fs2) zeros(1,(F_c/Fs2)/2)]; 
x_zoh2 = conv(x_comb2, zoh2, 'same');
plot(t_cut,x_zoh2(index1:index2));


% Question no. 9
subplot(2, 1, 1);
foh1 = triang(2*F_c/Fs1); 
x_foh1 = conv(x_comb1, foh1, 'same');
plot(t_cut,x_foh1(index1:index2));

subplot(2, 1, 2);
foh2 = triang(2*F_c/Fs2); 
x_foh2 = conv(x_comb2, foh2, 'same');
plot(t_cut,x_foh2(index1:index2));


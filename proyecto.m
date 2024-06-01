%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEPARADOR POR FILTRADO INVERSO
% FUNCIONA EN BLOQUES
% Parametros:
% Nombre del fichero de entrada: n_fichero
% Polo del filtro de preenfasis: polo
% Longitud de la ventana:            l_v
% Numero de coeficientes del modelo LPC: n_coef
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 close all;
 clear;

 n_fichero = 'aeiou2.wav';

[x,fs]=audioread(n_fichero);

l_v = 512;                                  % longitud de la ventana
desplaza = l_v/2;                           % desplazamiento
n_trozos = floor(length(x)/desplaza) -1;    % trozos sobre los que evaluar la LPC
n_coef=14;                                  % Numero coef LPC
n_formantes = min(5,n_coef/2-1);
trama = floor(n_trozos/2);              % Trama individual a pintar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_f_pulso=[n_fichero,'_pulso.wav'];
n_f_tracto=[n_fichero,'_tracto.wav'];
n_f_onda=[n_fichero,'_onda.wav'];

a = 1;
polo=0.95;
b = [1,-polo];

x_fil = filter(b,a,x);
x=x_fil;

matr_lpc=zeros(l_v/2+1,n_trozos);

frecs_f=zeros(n_coef/2,n_trozos);

Zi=zeros(1,n_coef); % Condiciones iniciales de los filtros
tracto = []; % señal de tracto vocal (respuesta al impulso del tracto)
pulso  = []; % señal de pulso glotico
recompuesto = [];
recompuesto2 = [];
recompuesto3 = [];
Zi_rec=zeros(1,n_coef); % Condiciones iniciales de los filtros para la señal recompuesta
matriz_a_lpc = zeros(n_coef, n_trozos);
for i=1:n_trozos
    trozo=x((i-1)*desplaza+1 : (i-1)*desplaza+l_v);
    umbral=0.01;
    %Coeficientes LPC
    a_lpc=real(lpc(trozo,n_coef));
    %Se guardan los coeficientes en una matriz para tener toda la
    %informacion del tracto
    matr_a_lpc(:,i)=a_lpc;
    %[H,f]=freqz(1,a_lpc,l_v,fs);
    [H,~]=freqz(1,a_lpc,l_v/2+1,fs);
    matr_lpc(:,i)=20*log10(abs(H)+eps);
    %Extracción de polos
    r=roots(a_lpc);
    % Selección de polos a efectos de representación de formantes
    r=r(imag(r)>umbral);
    frecs=sort(atan2(imag(r),real(r))*fs/(2*pi));
    frecs_f(1:length(frecs),i)=frecs;
    % Filtro inverso del traco vocal
    b_inv=a_lpc;
    a_inv=1;
    % Proceso de filtrado con condiciones iniciales
    [trozo_pulso, Zi]=filter(b_inv,a_inv,trozo(1:desplaza),Zi);
    pulso=[pulso, (trozo_pulso(1:desplaza))'];
    
    [H_tracto,~]=freqz(1,a_lpc,l_v,fs,'whole');
    trozo_tracto=real(ifft(H_tracto));
    tracto=[tracto, (trozo_tracto(1:desplaza))'];
end

% Recomposición de voz con las señales de pulso y tracto
pulso_r= [pulso(2000:length(pulso)-2000),pulso(2000:length(pulso)-2000),pulso(2000:length(pulso)-2000)];
pulso_r1= resample (pulso_r,6,8);
pulso_r2= resample (pulso_r,4,2);

for i=1:n_trozos-1
    trozo = pulso((i-1)*desplaza+1 : (i-1)*desplaza+l_v);
    a_lpc = matr_a_lpc(:,i);
    [trozo_recompuesto, Zi_rec]=filter(1,a_lpc,trozo(1:desplaza),Zi_rec);
    recompuesto = [recompuesto, (trozo_recompuesto(1:desplaza))];
end

for i=1:n_trozos-1
    trozo = pulso_r1((i-1)*desplaza+1 : (i-1)*desplaza+l_v);
    a_lpc = matr_a_lpc(:,i);
    [trozo_recompuesto, Zi_rec]=filter(1,a_lpc,trozo(1:desplaza),Zi_rec);
    recompuesto2 = [recompuesto2, (trozo_recompuesto(1:desplaza))];
end

for i=1:n_trozos-1
    trozo = pulso_r2((i-1)*desplaza+1 : (i-1)*desplaza+l_v);
    a_lpc = matr_a_lpc(:,i);
    [trozo_recompuesto, Zi_rec]=filter(1,a_lpc,trozo(1:desplaza),Zi_rec);
    recompuesto3 = [recompuesto3, (trozo_recompuesto(1:desplaza))];
end


figure(1); 
specgram(tracto,l_v/2,fs,hamming(256),0) 
title('Espectro LPC');
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
zlabel('Magnitud (dB)');

figure(2);
eje_t=(0:n_trozos-1)/(fs/desplaza);
eje_lpc=(0:l_v/2)*fs/l_v;
surf(eje_t,eje_lpc,matr_lpc);
shading interp;
title('Espectro LPC');
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
zlabel('Magnitud (dB)');

[B,f,t]=specgram(x,l_v,fs);
figure(3);
eje_f=1:l_v/2+1;
eje_f=eje_f-1;
eje_f=eje_f*fs/l_v;
trama_lpc=matr_lpc(:,trama);
trama_fourier=20*log10(abs(B(:,trama)));
trama_lpc=trama_lpc-max(trama_lpc);
trama_fourier=trama_fourier-max(trama_fourier);
a_v=1;
b_v=[1 -1];
D_envolv=filter(b_v,a_v,trama_lpc);
D_D_envolv=filter(b_v,a_v,D_envolv);

plot(eje_lpc,trama_lpc,'r',eje_f,trama_fourier,'b', eje_f, D_envolv, 'k', eje_f, D_D_envolv,'g');
title('fft y envolvente lpc')
xlabel('Frecuencia (Hz)');
ylabel('Amplitud (dB)');
legend('lpc', 'fft', 'Dlpc', 'DDlpc');


figure(4);
surf(matr_lpc-20*log10(abs(B)));
title('Diferencia entre espectro LPC y FFT')
xlabel('Tiempo (s)');
ylabel('Frecuencia (Hz)');
zlabel('Amplitud (dB)');
shading interp;



figure(5);
eje_t2=0:1/fs:(length(pulso)-1)/fs;
plot(eje_t2,pulso);
title('Pulso glotico');
xlabel('Tiempo (s)');
ylabel('Amplitud');

figure(6);
specgram(pulso,512,fs);
title('Espectro del pulso glotico');

% Integramos el pulso para obtener la onda
a=[1,-0.99];
b=0.01;
onda_glot=filter(b,a,pulso);
figure(7);
plot(eje_t2,onda_glot);
title('Onda gl�tica');
xlabel('Tiempo (s)');
ylabel('Amplitud');
figure(8);
specgram(onda_glot,512,fs);
title('Espectro de la onda glotica');


figure(9);
eje = 1:length(onda_glot);
plot(eje,pulso,'b', eje, onda_glot*100,'r');
title('Pulso y onda');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('pulso','onda');

pulso=pulso/max(abs(pulso));
onda_glot=onda_glot/max(abs(onda_glot));

pulso=pulso/max(abs(pulso));
audiowrite(n_f_pulso,pulso,fs);
audiowrite(n_f_onda,onda_glot,fs);
tracto=tracto/max(abs(tracto));
audiowrite(n_f_tracto,tracto,fs);



sound(x,fs);
pause;
figure(10);
spectrogram(x,512,256,512,fs,'yaxis');
sound(recompuesto,fs);
pause;
figure(11);
spectrogram(recompuesto,512,256,512,fs,'yaxis');
sound(recompuesto2,fs);
pause;
figure(12);
spectrogram(recompuesto2,512,256,512,fs,'yaxis');

sound(recompuesto3,fs);
pause;
figure(13);
spectrogram(recompuesto3,1024,256,1024,fs,'yaxis');



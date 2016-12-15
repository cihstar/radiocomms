function BER = DigiTrans_with_IF(SNR, plot_on, diff_on)
% Matlab script file:	DigiTrans.m
%	
% This file simulates a transmission of 16-QAM data. 
%
% All functions which are called from this script file are available
% from http://www.ecs.soton.ac.uk/~sw1/ez622/ez622.html
%
% S. Weiss, 10/11/2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bit stream generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LB = 1000;			% number of bits
B = BitStream(LB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion to 16-QAM symbol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = QPSK_mod(B);			% to implement 2^(2*Nbits)-QAM

if(diff_on == 1)
    X = Differential_Modulation(X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upsampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1024; 	% oversampling factor for transmitted data
Xup = zeros(1,length(X)*N);
for i = 1:N
    Xup(i:N:end) = X;
end

if(plot_on==1)
figure(8); clf;
plot(Xup(20:end),'.');
title('constellation'); xlabel('I'); ylabel('Q');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulation onto carrier and transmit filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fc = 70e3; %70kHz carrier
t= 0:1/(length(Xup)):(1-1/(length(Xup)));

s =  (real(Xup) .* cos(2*pi*Fc*t) + imag(Xup) .* sin(2*pi*Fc*t));
s_star = s;

if(plot_on==1)
figure(3);
periodogram(s_star,[],[],length(s_star),'power','centered');
title('s^*');
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering with channel impulse response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = [1 0 0];			        % CIR (at N x symbol rate)
s_hat = filter(c,1,s_star);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additive white Gaussian noise (AWGN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_x = std(s_hat);
Ls = length(s_hat);
noise = (randn(1,Ls) + sqrt(-1)*randn(1,Ls))*sqrt(N)/sqrt(2);
s_hat = sqrt(2) * s_hat + sigma_x*10^(-SNR/20)*noise;
% line above WAS: (incorrectly) s_hat = s_hat + sigma_x*10^(-SNR/20)*sqrt(N)*noise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% receive filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bandpass filter signal
if(plot_on==1)
figure(4);  clf;
hold on
periodogram(s_hat);
title('IF filter');
<<<<<<< HEAD
order    = 5;
fcutlow  = 69000;
fcuthigh = 71000;
Fs1 = N*length(X);
[b,a]    = butter(order,[fcutlow,fcuthigh]/(Fs1/2), 'bandpass');
s2_hat = filter(b,a,s_hat);
[h,w] = freqz(b,a,100);
if(plot_on==1)
plot(w/pi,20*log10(abs(h)),'r');
title('BPF for received signal');
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF SUBSAMPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dec_fac = 64;

sample_rate = N/dec_fac;

Fs = N*length(X);
Fs_altDDC = sample_rate*length(X);   % Sampling frequency
adc = s2_hat(1:sample_rate:end);        % "sample" the signal 
if(plot_on == 1)
figure(5);
[Hys,Fys] = periodogram(adc,[],[],Fs,'power','centered');
plot(Hys, Fys);
title('ADC periodogram');

Nc = length(Fys);
figure(7); clf;
periodogram(s2_hat,[],[],Fs,'power','centered');
title('Received signal');
clear y;
hold on;
plot((-(ceil(Nc/2*90)-1):floor(Nc/2*90))/Nc*Fs_altDDC/1000, ...
        repmat(10*log10(Hys),90,1),'r:');
axis([0 40 -300 0])
legend('Input of A/D Converter','Aliased Output of A/D Converter', ...
    'Location','NorthEast');
end

t2 = t(1:sample_rate:end);

s2_star = 2*(adc .* cos(2*pi*6e3*t2) + j*adc.* sin(2*pi*6e3*t2));


s2_star_i = real(s2_star);
s2_star_q = imag(s2_star);

N3   = 120;
Fs  = length(s2_star);
Fp  = 1e3;
Ap  = 0.01;
Ast = 80;

Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);

NUM = firceqrip(N3,Fp/(Fs/2),[Rp Rst],'passedge');
s2_star = filter(NUM,1,s2_star_i) + j*filter(NUM,1,s2_star_q);
if(plot_on==1)
figure(8); clf;
periodogram(s2_star,[],[],Fs_altDDC,'power','centered');
title('demod & filtered signal');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(diff_on == 1)
    s2_star = Differential_Demodulation(s2_star, sample_rate*4);
end

if(plot_on == 1)
figure(24); clf;
plot(s2_star(20:end),'.');
title('constellation recieved after differential detection'); xlabel('I'); ylabel('Q');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EYE = zeros(32,200); 
EYE(:) = s2_star(sample_rate*20+1:sample_rate*20+32*200)';
if(plot_on==1)
figure(1); clf;
plot(real(EYE));	% I-component only
title('eye diagram of received data');
xlabel('wrapped time'); ylabel('I-component amplitude');
end

EYE = zeros(32,200); 
EYE(:) = s2_star(sample_rate*20+1:sample_rate*20+32*200)';
if(plot_on==1)
figure(9); clf;
plot(imag(EYE));	% Q-component only
title('eye diagram of received data');
xlabel('wrapped time'); ylabel('Q-component amplitude');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampling at symbol rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ninit = 2;%EstimateNinit(s2_star, sample_rate*4, 400, 0.03, 0.9, 0.03); 			% determine sampling point (0<Ninit<=N)
X_hat = 6*s2_star(Ninit:sample_rate*4:end);        % "sample" the signal 

% plot received constellation
if(plot_on==1)
figure(2); clf;
plot(X_hat(20:end),'.');
title('constellation'); xlabel('I'); ylabel('Q');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion from QPSK to bits stream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B2 = QPSK_demod(X_hat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate bit errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[corr, lag] = xcorr(B2, B);%perform cross correlation
[~,I] = max(abs(corr)); %The index with the highest value of correlation
delayB = lag(I); % to compensate dalays in channel & TX/RX
diff = B(1:end-delayB) - B2(delayB+1:end);
BER = sum(abs(diff))/(length(B)-delayB);
disp(sprintf('bit error probability = %f',BER));
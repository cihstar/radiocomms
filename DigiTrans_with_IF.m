function BER = DigiTrans_with_IF(SNR, plot_on, diff_on, if_on)
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

timespan = 1;
bps = LB/timespan;

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
N = 1024*timespan; 	% oversampling factor for transmitted data
Xup = zeros(1,length(X)*N);
Xup(1:N:end) = X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulation onto carrier and transmit filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = sqrt(N)*firrcos(N/timespan,1/(N/timespan),.5,2,'rolloff','sqrt');
s = filter(h,1,Xup);

if(if_on == 1)
    Fc = 70e3; %70kHz carrier
    t = 0:timespan/(length(s)):timespan - 1/length(s);   
    s = (real(s) .* cos(2*pi*Fc*t) + imag(s) .* sin(2*pi*Fc*t));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtering with channel impulse response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = [1 0 0];			        % CIR (at N x symbol rate)
s_hat = filter(c,1,s);
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
%s2_star = filter(h,1,s_hat);   
if(if_on == 1)
    %bandpass filter signal
    order    = 4;
    fwidth = bps*2;
    Fs1 = length(s_hat)/ timespan;
    
    
    
    [b,a] = butter(order,[Fc-fwidth/2,Fc+fwidth/2]/(Fs1/2), 'bandpass');
    s2_hat = filter(b,a,s_hat);  
    
    if(plot_on==1)
        figure(1);  clf;
        [pxx, f] = periodogram(s_hat,[],[],(length(s_hat))/timespan,'centered');
        title('IF filter');
        [h,w] = freqz(b,a,100);
        plot(f,10*log10(pxx),'b',w*(Fs1/(2*pi)),20*log10(abs(h)),'r');
        axis([-inf inf -200 0])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IF SUBSAMPLING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dec_fac = 32; %*timespan;
    sample_rate = N/(dec_fac*timespan);

    Fs = length(s2_hat)/timespan
    Fs_altDDC = Fs/sample_rate   % Sampling frequency
    adc = s2_hat(1:sample_rate:end);        % "sample" the signal 
    
    if(plot_on == 1)
        [Hys,Fys] = periodogram(adc,[],[],Fs,'power','centered');

        Nc = length(Fys);
        figure(2); clf;
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

    %RECEIVE FILTERING
    
    t2 = t(1:sample_rate:end);
    
    s2_star = 2*(adc .* cos(2*pi*6e3*t2) + j*adc.* sin(2*pi*6e3*t2));
    
    h = sqrt(N)*firrcos((N/dec_fac),1/(N/dec_fac),.5,2,'rolloff','sqrt');
else
    s2_star = s_hat;
    h = sqrt(N)*firrcos(N,1/N,.5,2,'rolloff','sqrt');
end

s2_star = filter(h,1,s2_star);
if(if_on==1)
    s2_star_i = real(s2_star);
    s2_star_q = imag(s2_star);

    %low filter signal
    if(plot_on==1)
        figure(3);  clf;
        hold on
        [pxx, f] = periodogram(s2_star_i,[],[],Fs_altDDC,'centered');
        title('LPF filter');
    end
    
    order    = 11;
    fcutlow  = bps*2;
    [b,a]    = butter(order,[fcutlow]/Fs_altDDC);
    s2_star = filter(b,a,s2_star_i) + 1i* filter(b,a,s2_star_q);
    
    if(plot_on==1)
        [h,w] = freqz(b,a,100);
        plot(f,10*log10(pxx),'b',w*(Fs_altDDC/(2*pi)),20*log10(abs(h)),'r');
        axis([-inf inf -200 0])
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(if_on == 1)
    samples_per_symbol = (N/dec_fac)
else
    samples_per_symbol = N;
end

if(diff_on == 1)
    s2_star = Differential_Demodulation(s2_star, samples_per_symbol);
    if(plot_on == 1)
        figure(4); clf;
        plot(s2_star(65:end),'.');
        title('constellation recieved after differential detection'); xlabel('I'); ylabel('Q');
    end
else
    s2_star = s_hat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_on==1)
    eyediagram(s2_star(samples_per_symbol*10:end),samples_per_symbol*2,samples_per_symbol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampling at symbol rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ninit = 24;
X_hat = s2_star(Ninit:samples_per_symbol:end);        % "sample" the signal 
length(X_hat)

% plot received constellation
if(plot_on==1)
    figure(6); clf;
    plot(X_hat(40:end),'.');
    title('Received constellation'); xlabel('I'); ylabel('Q');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion from QPSK to bits stream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B2 = QPSK_demod(X_hat(20:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate bit errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_t = B(39:end);

[corr, lag] = xcorr(B2, B_t);%perform cross correlation
[~,I] = max(abs(corr)); %The index with the highest value of correlation
delayB = lag(I); % to compensate dalays in channel & TX/RX
diff = B_t(1:end-delayB) - B2(delayB+1:end);
BER = sum(abs(diff))/(length(B)-delayB);
disp(sprintf('bit error probability = %f\tNinit = %d',BER, Ninit));
end
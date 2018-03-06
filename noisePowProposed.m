%%%% propose SPP algorithm to estimate the spectral noise power
%%%% papers: "Unbiased MMSE-Based Noise Power Estimation with Low Complexity and Low Tracking Delay", IEEE TASL, 2012 
%%%% "Noise Power Estimation Based on the Probability of Speech Presence", Timo Gerkmann and Richard Hendriks, WASPAA 2011
%%%% Output:  noisePowMat:  matrix with estimated noise power for each frame
clear all; close all;
% [s,fs]=wavread('Dev_2ch_3src_Sq_Ce_A_Src1.wav');
% s=s(1:80000,1);
% 
% s=s./max(s);
% t=1/fs.*(1:length(s));
% 
% figure; plot(t,s); 
% title('Clean Speech'); xlabel('Time (sec)');

%% white noise
% n= rand(length(s),1); title('White Noise');

%% babble noise

% [n,fs2]=wavread('Dev_2ch_3src_Sq_Ce_A_Noise.wav');
% n=n(1:80000,1);
% n=n./max(n);
% % figure;plot(t,n);
% % title('Babble Noise'); xlabel('Time (sec)');
% %% plot noisy signal
% noisy= s + 0.9*n;
% wavwrite(noisy,fs,'noisyJapanese.wav');
% 
% figure; plot(t,noisy);
% xlabel('Time (sec)'); title('Noisy signal');

[noisy,fs]=wavread('noisyJapanese.wav');
%% some constants
frLen   = 32e-3*fs;  % frame size
fShift  = frLen/2;   % fShift size
nFrames = floor(length(noisy)/fShift)-1; % number of frames

anWin  = hanning(frLen,'periodic'); %analysis window

%% allocate some memory
noisePowMat = zeros(frLen/2+1,nFrames);
sig=zeros(frLen/2+1,nFrames);
%% initialize
noisePow = init_noise_tracker_ideal_vad(noisy,frLen,frLen,fShift, anWin); % This function computes the initial noise PSD estimate. It is assumed that the first 5 time-frames are noise-only.
noisePowMat(:,1)=noisePow;
%
PH1mean  = 0.5;
alphaPH1mean = 0.9;
alphaPSD = 0.8;

%% constants for a posteriori SPP
    q          = 0.3; % a priori probability of speech presence:
    priorFact  = q./(1-q);
    % optimal fixed a priori SNR for SPP estimation
    ksiOptDb    = 15;    
    ksiOpt      = 10.^(ksiOptDb./10);
    logGLRFact = log(1./(1+ksiOpt));
%   logGLRFact = log((1+xiOpt));
    GLRexp     = ksiOpt./(1+ksiOpt);    
    c=sqrt(pi)/2;
%%

shat=zeros(size(noisy));

for indFr = 1:nFrames
    indices       = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;
    noisy_frame   = anWin.*noisy(indices);
    noisyDftFrame1 = fft(noisy_frame,frLen);
    noisyDftFrame = noisyDftFrame1(1:frLen/2+1);
	noisyPer= abs(noisyDftFrame).^2;
   % noisyPer = noisyDftFrame.*conj(noisyDftFrame);
    snrPost1 =  noisyPer./(noisePow);% a posteriori SNR based on old noise power estimate

    %% noise power estimation
       
    % GLR from 7.191 Loizou page 258
	GLR     = priorFact .* exp(min(logGLRFact + GLRexp.*snrPost1,200));	
    % PH1 from 7.187 Loizou page 258
    PH1     = GLR./(1+GLR); % a posteriori speech presence probability    
    %% Recursive smoothing of P(H1 | y) over time
    %  equation (23)
	PH1mean  = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
    
    %% If the smoothed quantity is larger than 0.99 
    %  the update is stagnated
    %  equation (24)
	stuckInd = PH1mean > 0.99;
	PH1(stuckInd) = min(PH1(stuckInd),0.99);
    
    %% soft weighting between the noisy observation |y|^2
    %  and the previous noise power estimate (?_?)^2
    %  equation (22)
	estimate = (1-PH1) .* noisyPer + PH1 .* noisePow   ;
    
    %% recursive smoothing equation 
    %  equation (8)
	noisePow = alphaPSD *noisePow+(1-alphaPSD)*estimate;
	noisePowMat(:,indFr) = noisePow;
    noiseAng(:,indFr)=angle(noisyDftFrame1);
    
    
    
    vk=ksiOpt.*snrPost1./(1+ksiOpt);
    j0=besseli(0,vk/2);
    j1=besseli(1,vk/2);

 hw=( c * vk.^0.5 .* exp(-0.5*vk) ./snrPost1 ).*((1+vk).*j0+vk.*j1);
    sig(:,indFr) =PH1.* hw.*abs(noisyDftFrame);
    sig2=[sig;flipud(conj(sig(2:end-1,:)))];
        xi_w=real( ifft( sig2(:,indFr) .* exp(1i*angle(noisyDftFrame1)),512));  
   shat(indices)= shat(indices) +  xi_w;

end

plot(shat)

% %%
% % consider the 128 symmetric components 
%  sig2=[sig;flipud(conj(sig(2:end-1,:)))];
% % %noiseAng2=[noiseAng;flipud(conj(noiseAng(2:end-1,:)))];
% % noiseMatrix=noisePowMat2.*exp(1i*noiseAng);
%  sTime= real(ifft(sig2.* exp(1i*noiseAng)+eps));
%  shat=zeros(size(noisy));
% % 
% for indFr = 1:nFrames
%  indices     = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;
%  shat(indices)= shat(indices) +   sTime(:,indFr);
% end
% % 
% % figure; plot(t,nhat);
% % xlabel('Time (sec)'); title('Estimated noise signal');
% 
% shat=noisy-nhat;
% figure; plot(t,shat);
% xlabel('Time (sec)'); title('Estimated clean signal');

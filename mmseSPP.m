%function xfinal= mmseSPP(filename,outfile)

clear all;close all;
filename='noisyJapanese.wav';
outfile='out_mmse.wav';
% e.g.     mmse('noisyJapanese.wav','out_mmse.wav');

%  References:
%   [1] Ephraim, Y. and Malah, D. (1985). Speech enhancement using a minimum 
%       mean-square error log-spectral amplitude estimator. IEEE Trans. Acoust., 
%       Speech, Signal Process., ASSP-23(2), 443-445.

[noisy, fs]= wavread( filename);	
 t=1/fs.*(1:length(noisy));
plot(t,noisy); title('noisy');
% =============== Initialize variables ===============

frLen   = 32e-3*fs;  % frame size
fShift=frLen/2;

win=hanning(frLen);  % define window
win = win*fShift/sum(win);  % normalize window for equal level output 

%--- allocate memory and initialize various variables
nFrames=floor(length(noisy)/fShift)-1;
xfinal=zeros(size(noisy));

% --------------- Initialize parameters ------------
aa=0.98;
c=sqrt(pi)/2;
qk=0.3;
priorFact=(1-qk)/qk;
ksi_min=10^(-25/10); % note that in Chap. 7, ref. [17], ksi_min (dB)=-15 dB is recommended
noise_muProp=ProposedOrig(noisy,fs);
noisePowMat=[noise_muProp;flipud(conj(noise_muProp(2:end-1,:)))];

%===============================  Start Processing =======================================================
%

for indFr=1:nFrames
    indices       = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;
    noisy_frame=win.*noisy(indices);
    noisyDftFrame1=fft(noisy_frame,frLen);
    sig=abs(noisyDftFrame1); % compute the magnitude
    noisyPer=sig.^2;
    noisePow=noisePowMat(:,indFr);
    snrPost1=min(noisyPer./noisePow,40);  % posteriori SNR
    
    % Loizou 7.55
    if indFr==1
        ksi=aa+(1-aa)*max(snrPost1-1,0);
    else
        ksi=aa*noisyPerPrevious./noisePow + (1-aa)*max(snrPost1-1,0);     
        % decision-direct estimate of a priori SNR
        ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
    end
    
    vk=ksi.*snrPost1./(1+ksi);
    j0=besseli(0,vk/2);
    j1=besseli(1,vk/2);
    
        B=(1+vk).*j0+vk.*j1;
        hw=( c*(vk.^0.5).* exp(-0.5*vk) ./snrPost1 ).*B;    

    % --- estimate speech presence probability
        GLR=priorFact*exp(vk)./(1+ksi);
        pSAP=GLR./(1+GLR);
        
      
        sig=sig.*hw.*pSAP;
 
    noisyPerPrevious=sig.^2;  % save for estimation of a priori SNR in next frame
    xi_w=real( ifft( sig .* exp(1i*angle(noisyDftFrame1)),frLen));  
   xfinal(indices)= xfinal(indices) +  xi_w;
end
%========================================================================================
figure;
 plot(t,xfinal);title('estimate');
wavwrite(xfinal,fs,16,outfile);


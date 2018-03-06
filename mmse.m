function [xfinal,Srate]=mmse(filename,outfile)

%
%  Implements the MMSE algorithm [1].
% 
%  Usage:  mmse(noisyFile, outputFile, SPU)
%           
%         infile - noisy speech file in .wav format
%         outputFile - enhanced output file in .wav format
%         SPU  - if 1, includes speech-presence uncertainty
%                if 0, doesnt include speech-presence uncertainty
%  
%
%  Example call:  mmse('sp04_babble_sn10.wav','out_mmse.wav',1);
%
%  References:
%   [1] Ephraim, Y. and Malah, D. (1985). Speech enhancement using a minimum 
%       mean-square error log-spectral amplitude estimator. IEEE Trans. Acoust., 
%       Speech, Signal Process., ASSP-23(2), 443-445.
%   
% Authors: Philipos C. Loizou
%
% Copyright (c) 2006 by Philipos C. Loizou
% $Revision: 0.0 $  $Date: 10/09/2006 $
%-------------------------------------------------------------------------




[x, Srate, bits]= wavread( filename);	


% =============== Initialize variables ===============

len=floor(20*Srate/1000); % Frame size in samples
if rem(len,2)==1, len=len+1; end;
PERC=50; % window overlap in percent of frame size
len1=floor(len*PERC/100);
len2=len-len1;

win=hanning(len);  % define window
win = win*len2/sum(win);  % normalize window for equal level output 

% Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%
nFFT=2*len;
j=1;
noise_mean=zeros(nFFT,1);
for k=1:6
    noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
    j=j+len;
end
noise_mu=noise_mean/6;
noise_mu2=noise_mu.^2;

%--- allocate memory and initialize various variables

k=1;
img=sqrt(-1);
x_old=zeros(len1,1);
Nframes=floor(length(x)/len2)-1;
xfinal=zeros(Nframes*len2,1);

% --------------- Initialize parameters ------------
%
k=1;
aa=0.98;
frLen   = 32e-3*Srate;  % frame size
fShift=frLen/2;
nFFT=frLen;
win=hanning(frLen);  % define window
win = win*fShift/sum(win);  % normalize window for equal level output 

%--- allocate memory and initialize various variables
nFrames=floor(length(x)/fShift)-1;
xfinal=zeros(size(x));
c=sqrt(pi)/2;
qk=0.3;
qkr=(1-qk)/qk;
ksi_min=10^(-25/10); % note that in Chap. 7, ref. [17], ksi_min (dB)=-15 dB is recommended
 noise_muProp=ProposedOrig(x,Srate);
 noisePowMat=[noise_muProp;flipud(conj(noise_muProp(2:end-1,:)))];

%===============================  Start Processing =======================================================
%
for indFr=1:nFrames
    indices       = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;

  insign=win.*x(indices);
    %--- Take fourier transform of  frame
    %
    spec=fft(insign,nFFT);
    sig=abs(spec); % compute the magnitude
    sig2=sig.^2;
    noise_mu2=noisePowMat(:,indFr);

    gammak=min(sig2./noise_mu2,40);  % posteriori SNR
    if indFr==1
        ksi=aa+(1-aa)*max(gammak-1,0);
    else
        ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
        % decision-direct estimate of a priori SNR
        ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
    end


    % ===end of vad===

    vk=ksi.*gammak./(1+ksi);
   j0=besseli(0,vk/2);
   j1=besseli(1,vk/2);
  
        C=exp(-0.5*vk);
        A=((c*(vk.^0.5)).*C)./gammak;
        B=(1+vk).*j0+vk.*j1;
        hw=A.*B;
    


    % --- estimate speech presence probability
    %
    
        evk=exp(vk);
        Lambda=qkr*evk./(1+ksi);
        pSAP=Lambda./(1+Lambda);
       sig=sig.*hw.*pSAP;
    
        sig=sig.*hw;
    
    
    Xk_prev=sig.^2;  % save for estimation of a priori SNR in next frame


    xi_w=real( ifft( sig .* exp(1i*angle(spec)),nFFT));  
   xfinal(indices)= xfinal(indices) +  xi_w;
    
end
%========================================================================================


wavwrite(xfinal,Srate,16,outfile);


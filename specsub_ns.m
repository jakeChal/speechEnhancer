function xfinal=specsub_ns(filename,outfile)

% e.g. 
% xfinal=specsub_ns('noisyJapanese.wav','estJap.wav')
[x,Srate,nbits]=wavread(filename);


% =============== Initialize variables ===============
%

len=floor(32*Srate/1000); % Frame size in samples
fShift=len/2;
if rem(len,2)==1, len=len+1; end;
PERC=50; % window overlap in percent of frame size
len1=floor(len*PERC/100);
len2=len-len1; 


alpha=2.0; % power exponent
FLOOR=0.002;
win=hamming(len); %tukey(len,PERC);  % define window


%--- allocate memory and initialize various variables
   

nFFT=len;
%nFFT=2*len;
img=sqrt(-1);
Nframes=floor(length(x)/len2)-1;
%xfinal=zeros(Nframes*len2,1);

%===============================  Start Processing =======================================================
%
xfinal=zeros(size(x));

noise_ps = ProposedOrig(x,Srate);
for n=1:Nframes 
       indices       = (n-1)*fShift+1:(n-1)*fShift+len;
 insign=win.*x(indices);
  % insign=win.*x(k:k+len-1);     %Windowing  
   spec=fft(insign,nFFT);     %compute fourier transform of a frame
   sig=abs(spec); % compute the magnitude

   noise_psd = [noise_ps(:,n);flipud(noise_ps(2:end-1,n))];
   noise_mu=sqrt(noise_psd);  % magnitude spectrum
   % ---------------------------------------------------------
   
    %save the phase information for each frame.
    theta=angle(spec);  
   
   SNRseg=10*log10(norm(sig,2)^2/norm(noise_mu,2)^2);
   
   if alpha==1.0
      beta=berouti1(SNRseg);
   else
     beta=berouti(SNRseg);
  end
   
   
   %&&&&&&&&&
   sub_speech=sig.^alpha - beta*noise_mu.^alpha;
   diffw = sub_speech-FLOOR*noise_mu.^alpha;
   
   % Floor negative components
   z=find(diffw <0);  
   if~isempty(z)
      sub_speech(z)=FLOOR*noise_mu(z).^alpha;
   end
   
    
   sub_speech(nFFT/2+2:nFFT)=flipud(sub_speech(2:nFFT/2));  % to ensure conjugate symmetry for real reconstruction
   %multiply the whole frame fft with the phase information
   x_phase=(sub_speech.^(1/alpha)).*(cos(theta)+img*(sin(theta)));
  
   
   % take the IFFT 
   xi=real(ifft(x_phase));         
 
  % --- Overlap and add ---------------
  % 
  
   xfinal(indices)= xfinal(indices) +   xi;
  %xfinal(k:k+len2-1)=x_old+xi(1:len1);
%  x_old=xi(1+len1:len);
  
% k=k+len2;
end
%========================================================================================

wavwrite(xfinal,Srate,16,outfile);

%-------------------------------- E N D --------------------------------------
function a=berouti1(SNR)

if SNR>=-5.0 && SNR<=20
   a=3-SNR*2/20;
else
   
  if SNR<-5.0
   a=4;
  end

  if SNR>20
    a=1;
  end
  
end

function a=berouti(SNR)

if SNR>=-5.0 && SNR<=20
   a=4-SNR*3/20; 
else
   
  if SNR<-5.0
   a=5;
  end

  if SNR>20
    a=1;
  end
  
end

snrPostdB=10*log10(snrPost1);

snrPostMin=min(snrPostdB);

xaxis=linspace(snrPostMin,-snrPostMin,100);

q2=0.1;
    priorFact2  = q2./(1-q2);

 ksiOptDb2    =15;    
    ksiOpt2      = 10.^(ksiOptDb2./10);
    logGLRFact2 = log(1./(1+ksiOpt2));
for i=1:100
    
    Lamda    = priorFact2 .* exp(min(logGLRFact2 + GLRexp.*xaxis(i),200));	
     PH1peir (i)    = Lamda./(1+Lamda);
end

plot(xaxis,PH1peir)
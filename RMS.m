function rms=synergi(Window,tempsimu,y);
%% calculo valor RMS 1
 j=0; k=0; RR=0;
 for i=1:tempsimu
    j=j+1;
    RR=RR+y(i)^2;
    if(j>Window-1)
    k=k+1;    
    RMS(k)=sqrt((RR/j));
    j=0; RR=0;
    end
 end
 rms=RMS;
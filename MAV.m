function mav=MAV(Window,tempsimu,y);
 j=0; k=0; RR=0;
  %% Calculo valor MAV
 for i=1:tempsimu
    j=j+1;
    RR=RR+abs(y(i));
    if(j>=Window)
        k=k+1;    
        MAVs(k)=(RR/j);
        j=0; RR=0;
    end
  
 end
 mav=MAVs;    
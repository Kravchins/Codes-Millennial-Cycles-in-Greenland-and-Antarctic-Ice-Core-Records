function [coi]=coi_calc(fn_resamp,dt)
ff=(4*pi)/(6+sqrt(2+6^2));
coi=ff/sqrt(2);
n=length(fn_resamp);
coi=coi*dt*[1E-5,1:((n+1)/2-1),fliplr((1:(n/2-1))),1E-5];
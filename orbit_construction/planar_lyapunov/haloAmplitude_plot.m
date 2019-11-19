function haloAmplitude_plot

% Shane Ross (revised 7.13.04)

gam = 1.497531218885853e+06 ; % scale factor for Sun-Earth L1

l1 =-15.9650314   ;
l2 =  1.740900800 ;
Delta = 0.29221444425 ;

i = 0 ;
for Az = 0:0.01e5/gam:5e5/gam,
    i = i + 1 ;

    x(i) = gam*Az ;
    y(i) = gam*sqrt( (-Delta-l2*Az^2)/l1 ) ;

end


h=plot(x/1e5,y/1e5,'k') ;
axis([0 5e5 1.9e5 2.7e5]/1e5) ;
xlabel('A_z');
ylabel('A_x');

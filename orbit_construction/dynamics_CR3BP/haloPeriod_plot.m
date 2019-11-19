function haloAmplitude_plot

% Shane Ross (revised 7.13.04)

gam = 1.497531218885853e+06 ; % scale factor for Sun-Earth L1

l1 =-15.9650314   ;
l2 =  1.740900800 ;
Delta = 0.29221444425 ;

s1 =-8.246608317e-1 ;
s2 = 1.210985938e-1 ;
 
wp = 2.086453455 ;

i = 0 ;
for Az = 0:0.01e5/gam:5e5/gam,
    i = i + 1 ;

    Ax = sqrt( (-Delta-l2*Az^2)/l1 ) ;
    nu = 1 + s1*Ax^2 + s2*Az^2 ;

    x(i) = gam*Az ;
    y(i) = 2*pi/(wp*nu) ;

end


h=plot(x/1e5,y*365.25/(2*pi),'k') ;
axis([0 5 177 177.8]) ;
xlabel('A_z');
ylabel('T');

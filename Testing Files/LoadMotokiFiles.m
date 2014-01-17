motmodel = ncdfread('/Users/JacobWenegrat/Documents/LWM/Integration/output/SSH_out.nc');

motin = ncdfread('/Users/JacobWenegrat/Documents/LWM/Wforc/ECMWF_TAUX_IO.nc');
%%

for i=1:1682
   
    pcolor(motmodel.lon, motmodel.lat, double(squeeze(motmodel.SSH(:,:,i))'));
    shading interp;
    ylim([-15 15]);
    drawnow;
end
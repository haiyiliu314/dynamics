Ebind = 4.18;
p_fort_input = importdata('/home/liuhai/dynamics/with_coulomb_Fortran/p_freq.dat', ',');
E_fort_input = importdata('/home/liuhai/dynamics/with_coulomb_Fortran/E_freq.dat', ',');
p_fort = p_fort_input(:,1) + 1i*p_fort_input(:,2);
E_fort = E_fort_input(:,1) + 1i*E_fort_input(:,2);
plot(E001'*Ebind+4*Ebind,imag(p_fort./E_fort))
plot(E001'*Ebind+4*Ebind, abs(imag(p_fort./E_fort)/max(imag(p_fort./E_fort)) - E1/max(E1)))























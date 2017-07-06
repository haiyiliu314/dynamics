Ebind = 4.18;
t_end1 = 7;  %ps
Nt = 7000;
dt = t_end1/Nt;
p_fort_input = importdata('/home/liuhai/dynamics/Formal_V1_0703/run3/p_freq.dat', ',');
E_fort_input = importdata('/home/liuhai/dynamics/Formal_V1_0703/run3/E_freq.dat', ',');
p_fort = p_fort_input(:,1) + 1i*p_fort_input(:,2);
E_fort = E_fort_input(:,1) + 1i*E_fort_input(:,2);
test = importdata('/home/liuhai/Downloads/simpleRK4(12).dat', ',');
num = (1:280)*2500;

figure
plot(E001'*Ebind+4*Ebind, imag(p_fort./E_fort)/max(imag(p_fort./E_fort)))
hold on
plot(E001'*Ebind+4*Ebind, E1/max(E1))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('susceptibility')
legend('Haiyi', 'Ben')
head = sprintf('%s normalized susceptibility', date);
title(head)
print('susceptibility', '-dpdf')

% p_fort_input = importdata('/home/liuhai/dynamics/dynamics-20170619-seems-right/p_freq.dat', ',');
% E_fort_input = importdata('/home/liuhai/dynamics/dynamics-20170619-seems-right/E_freq.dat', ',');
% p_fort1 = p_fort_input(:,1) + 1i*p_fort_input(:,2);
% E_fort1 = E_fort_input(:,1) + 1i*E_fort_input(:,2);
% plot(E001'*Ebind+4*Ebind, (imag(p_fort./E_fort)./max(imag(p_fort./E_fort))))
% hold on
% plot(E001'*Ebind+4*Ebind, ( imag(p_fort1./E_fort1)./max(imag(p_fort1./E_fort1))))
% plot(E001'*Ebind+4*Ebind, (imag(p_fort./E_fort)./max(imag(p_fort./E_fort)) - imag(p_fort1./E_fort1)./max(imag(p_fort1./E_fort1))))
figure
semilogy(E001'*Ebind+4*Ebind, abs(imag(p_fort./E_fort)/max(imag(p_fort./E_fort)) - E1/max(E1))./(abs(E1/max(E1))))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('relative error')
head = sprintf('%s Comparison on susceptibility', date);
title(head)
print('susceptibility_error', '-dpdf')

pt1 = importdata('/home/liuhai/dynamics/Formal_V1_0703/run3/pt.dat', ',');
figure
subplot(2,1,1)
semilogy((num)*dt - 5, abs((real(pt1(num,1)) - test(:,3))./test(:,3)))
xlabel('time(ps)')
ylabel('relative error')
head = sprintf('%s Comparison on Re[P(t)]', date);
title(head)
subplot(2,1,2)
semilogy((num)*dt - 5, abs(((pt1(num,2)) - test(:,4))./test(:,4)))
xlabel('time(ps)')
ylabel('relative error')
head = sprintf('%s Comparison on Im[P(t)]', date);
title(head)
print('pt_error', '-dpdf')

Et1 = importdata('/home/liuhai/dynamics/Formal_V1_0703/run3/Et.dat', ',');
figure
plot((num-1)*dt - 5, (real(Et1(num-1,1)') - test(2:end,5)') ./ test(2:end,5)')
xlabel('time(ps)')
ylabel('relative error')
head = sprintf('%s Comparison on E(t)', date);
title(head)
print('Et_error', '-dpdf')

figure
subplot(2,1,1)
pk_Re_t = importdata('/home/liuhai/Downloads/z_P_Re(1).dat',',');
pk_Re = pk_Re_t(:,2:end);
pk1 = importdata('/home/liuhai/dynamics/Formal_V1_0703/run3/pk.dat',',');
% pk2 = zeros(800, 280);
% for i = 0:279
%     pk2(:,i+1) = pk1((800*i+1):800*(i+1), 1)+ 1i*pk1((800*i+1): 800*(i+1), 2);
% end
plot(y, pk_Re(21,:))
hold on
plot(y, real(pk(:, 45)))
legend('Ben', 'Haiyi')
xlabel('y')
ylabel('Re[P_y]')
head = sprintf('%s Comparison on Re[P_y] at t = 0.625 ps', date);
title(head)
subplot(2,1,2)
semilogy(y, abs(pk_Re(21,:) - real(pk(:, 45)'))./(real(pk(:, 45)')))
xlabel('y')
ylabel('Relative difference on Re[P_y]')
print('Py0_625_error', '-dpdf')

figure
subplot(2,1,1)
pk_Re_t = importdata('/home/liuhai/Downloads/z_P_Re(1).dat',',');
pk_Re = pk_Re_t(:,2:end);
pk1 = importdata('/home/liuhai/dynamics/Formal_V1_0703/run3/pk.dat',',');
plot(y, pk_Re(41,:))
hold on
plot(y, real(pk(:, 65)))
legend('Ben', 'Haiyi')
xlabel('y')
ylabel('Re[P_y]')
 head = sprintf(' Comparison on Re[P_y] at t = 2.5 ps', date);
title(head)
subplot(2,1,2)
plot(y, (pk_Re(41,:) - real(pk(:, 65)'))./max(real(pk(:, 65))))
xlabel('y')
ylabel('Relative difference on Re[P_y]')
print('Py_2_5', '-dpdf')




test1 = importdata('/home/liuhai/Downloads/more_ft_final(3).dat', ',');
figure
subplot(2,1,1)
semilogy(E001'*Ebind+4*Ebind, abs(E_fort_input(:,1)/Ebind - test1(:,4))./test1(:,4))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Relative error on Re[E(\omega)]')
head = sprintf('Comparison on $E_{freq}$', date);
title(head,  'Interpreter', 'Latex')
subplot(2,1,2)
semilogy(E001'*Ebind+4*Ebind, abs(E_fort_input(:,2)/Ebind - test1(:,5)))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Absolute error on Im[E(\omega)]')
print('E_freq_error', '-dpdf')

figure
subplot(2,1,1)
plot(E001'*Ebind+4*Ebind,(E_fort_input(:,1)/Ebind - test1(:,4))./test1(:,4))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Relative error on Re[E(\omega)]')
head = sprintf('Comparison on E_{freq}', date);
title(head)
subplot(2,1,2)
plot(E001'*Ebind+4*Ebind, abs(E_fort_input(:,2)/Ebind - test1(:,5)))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Absolute error on Im[E(\omega)]')

figure
subplot(4,1,1)
plot(E001'*Ebind+4*Ebind,  (p_fort_input(:,1) ))
hold on
plot(E001'*Ebind+4*Ebind,  test1(:,2))
legend('Haiyi', 'Ben')
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Re[$P_{freq}]$', 'Interpreter', 'Latex')
head = sprintf('Comparison on Re[$P_{freq}$]', date);
title(head, 'Interpreter', 'Latex')
subplot(4,1,2)
semilogy(E001'*Ebind+4*Ebind, abs((p_fort_input(:,1)-test1(:,2))./test1(:,2)))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Relative difference on Re[$P_{freq}$]', 'Interpreter', 'Latex')
head = sprintf('Comparison on Re[$P_{freq}$]', date);
title(head,  'Interpreter', 'Latex')
subplot(4,1,3)
plot(E001'*Ebind+4*Ebind,  (p_fort_input(:,2) ))
hold on
plot(E001'*Ebind+4*Ebind,  test1(:,3))
legend('Haiyi', 'Ben')
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('$P_{freq}$', 'Interpreter', 'Latex')
head = sprintf('Comparison on Im[$P_{freq}$]', date);
title(head, 'Interpreter', 'Latex')
subplot(4,1,4)
semilogy(E001'*Ebind+4*Ebind, abs(p_fort_input(:,2)-test1(:,3))./test1(:,3))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Relative difference on Im[$P_{freq}$]', 'Interpreter', 'Latex')
head = sprintf('Comparison on Im[$P_{freq}$]', date);
title(head,  'Interpreter', 'Latex')
print('P_freq_error', '-dpdf')

figure
subplot(4,1,1)
plot(E001'*Ebind+4*Ebind,  (p_fort_input(:,1) ))
hold on
plot(E001'*Ebind+4*Ebind,  test1(:,2))
legend('Haiyi', 'Ben')
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Re[$P_{freq}]$', 'Interpreter', 'Latex')
head = sprintf('Comparison on Re[$P_{freq}$]', date);
title(head, 'Interpreter', 'Latex')
subplot(4,1,2)
plot(E001'*Ebind+4*Ebind, abs((p_fort_input(:,1)-test1(:,2))./test1(:,2)))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Relative difference on Re[$P_{freq}$]', 'Interpreter', 'Latex')
head = sprintf('Comparison on Re[$P_{freq}$]', date);
title(head,  'Interpreter', 'Latex')
subplot(4,1,3)
plot(E001'*Ebind+4*Ebind,  (p_fort_input(:,2) ))
hold on
plot(E001'*Ebind+4*Ebind,  test1(:,3))
legend('Haiyi', 'Ben')
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('$P_{freq}$', 'Interpreter', 'Latex')
head = sprintf('Comparison on Im[$P_{freq}$]', date);
title(head, 'Interpreter', 'Latex')
subplot(4,1,4)
plot(E001'*Ebind+4*Ebind, abs(p_fort_input(:,2)-test1(:,3))./test1(:,3))
xlabel('h\omega -(E_g+E_{1s})[meV]')
ylabel('Relative difference on Im[$P_{freq}$]', 'Interpreter', 'Latex')
head = sprintf('Comparison on Im[$P_{freq}$]', date);
title(head,  'Interpreter', 'Latex')

figure
coul_mat1 = importdata('/home/liuhai/dynamics/with_coulomb_Fortran/coul_mat.dat', ',');
coul_mat2 =  importdata('/home/liuhai/dynamics/with_coulomb_Fortran/CoulombM.dat', ',');
coul_mat2 = coul_mat2(:,2:end)/2/pi;
subplot(3,1,1)
mesh(coul_mat1 - coul_mat2)
zlabel('absolute difference')
subplot(3,1,2)
plot(diag(coul_mat1 - coul_mat2))
title('Error on diagonal terms')
ylabel('absolute difference')
subplot(3,1,3)
mesh(coul_mat1)
title('Coulomb matrix')

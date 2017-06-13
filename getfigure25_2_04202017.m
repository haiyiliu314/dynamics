clear all
f1 = importdata('f1.dat');
psquared1 = importdata('psquared1.dat');
f2 = importdata('f2.dat');
psquared2 = importdata('psquared2.dat');
f3 = importdata('f3.dat');
psquared3 = importdata('psquared3.dat');
dt = importdata('dt.dat');
tgrid = importdata('tgrid.dat');
time = dt * (0:(tgrid)) * 1e12;
time1 = 0:0.01:1;
y1 = max(f1)*ones(1,101);
y2 = max(f2)*ones(1,101);
y3 = max(f3)*ones(1,101);
figure
subplot(3,1,1)
h = area(time1, y1);
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on
plot(time,f1)
hold on
plot(time,psquared1, '--')
axis([0 2.5 -max(f1)*0.1 max(f1)*1.1])

subplot(3,1,2)
h = area(time1, y2);
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on
plot(time,f2)
hold on
plot(time,psquared2, '--')
axis([0 2.5 -max(f2)*0.1 max(f2)*1.1])

subplot(3,1,3)
h = area(time1, y3);
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on
plot(time,f3)
hold on
plot(time,psquared3, '--')
axis([0 2.5 -max(f3)*0.1 max(f3)*1.1])
xlabel('Time(ps)')
hAxis = axes('visible', 'off');
h1 = text(-0.12, 0.5, 'Excitation');
set(h1, 'fontsize', 12, 'rotation' , 90, 'HorizontalAlignment', 'center')
% ylabel('Excitation')
f = importdata('f.dat');
psquared = importdata('psquared.dat');
dt = importdata('dt.dat');
tgrid = importdata('tgrid.dat');
time = dt * (0:(tgrid)) * 1e12;
time1 = 0:0.01:1;
y1 = max(f)*ones(1,101);
h = area(time1, y1);
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on
plot(time,f)
hold on
plot(time,psquared, '--')
axis([0 2.5 -max(f)*0.1 max(f)*1.1])
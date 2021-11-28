clear; clc;close all

N    = 10;
ne   = 200;
nt   = 10000;
isnap= 10

f = fopen("OUTPUT/snapshots_V.bin","r");
u = fread(f,"float64");
u = reshape(u,nt/isnap,[]);

f = fopen("OUTPUT/snapshots_analytical_V.bin","r");
v = fread(f,"float64");
v = reshape(v,nt/isnap,[]);


figure()
for i=1:nt/isnap
    subplot(211)
    plot(u(i,:));
    subplot(212)
    plot(v(i,:))
    pause(.0000009)
end
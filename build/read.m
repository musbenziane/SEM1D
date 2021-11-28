clear;
f = fopen("OUTPUT/source.bin","r");
data = fread(f,"float64");
plot(data);
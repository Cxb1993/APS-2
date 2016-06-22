i=1;
data = importdata('circle.dat');

for r = 1:200
    for c = 1: 200
        VolFrac(r,c) = data(i);
        i=i+1;
        
    end
end

dt = 0.005; dx = 0.02;
t=0;

while (t<2)
for r = 2:199
    for c= 2: 199
        VolFrac_n1(r,c) = dt*(VolFrac(r,c)+(VolFrac(r,c+1)-VolFrac(r,c-1))/(2*dx));
    end
end
t=t+1;

plot
end
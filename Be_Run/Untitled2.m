x=1:10;
y = 0.25 -  heaviside(x-3) + heaviside(x-7);
plot(y)

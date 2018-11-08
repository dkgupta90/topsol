jL = 310;
j0 = -0.006;
beta = 16.4;
max = 0.0;
maxV = 0.0;

for V = 0.5:0.001:0.58
    P = V*(jL + j0*exp(beta*V));
    if P > max
        max = P;
        maxV = V
    end
end
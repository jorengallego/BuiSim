TLow = zeros(35224,1);
for k = 0:366
    for i = 1:33
        TLow((k*96)+i)= 290.15;
    end
    for i = 34:93
        TLow((k*96)+i) = 293.15;
    end
    for i = 94:96
        TLow((k*96)+i) = 290.15;
    end
end
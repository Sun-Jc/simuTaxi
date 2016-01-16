function solution = binarySearch(func,rg)
eps = 0.00001;
a = rg(1);
b = rg(2);
while abs(a-b) > eps
    c = (a + b) /2;
    if func(a) * func(c) < 0
        b = c;
    else
        a = c;
    end
end
solution = b;
end
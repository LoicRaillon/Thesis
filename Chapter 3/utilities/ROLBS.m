function out = ROLBS(n, initReg, lambda)

if(n < 2 || n > 11)
    disp('n must be between 1 and 11')
    return
end

reg                 = de2bi(initReg);
N                   = 2^n - 1;
x                   = zeros(N,1);
for i = 1:N
    if (n <= 4 || n == 6)
        reg(n+1)    = xor(reg(1),reg(2));
    elseif (n == 5 || n == 11)
        reg(n+1)    = xor(reg(1),reg(3));
    elseif (n == 7 || n == 10)
        reg(n+1)    = xor(reg(1),reg(4));
    end
    x(i)            = round(reg(1));
    reg(1:n)        = reg(2:n+1);
end
out                 = repelem(x',1,lambda);

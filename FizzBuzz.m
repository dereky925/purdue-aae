clear;clc

for i = 1:15

    if mod(i,3) == 0
        fprintf("Fizz\n");
    elseif mod(i,5) == 0
        fprintf("Buzz\n");
    else
        fprintf("%f\n",i);
    end

    if mod(i,3) == 0 && mod(i,5) == 0
        fprintf("FizzBuzz\n");
    end

end



function k=fib(n)  
    if n == 0  
      k=0;  
      return;  
    elseif n == 1  
      k=1;  
    return;  
    else  
      k=fib(n - 1) + fib(n - 2);  
    return;  
    end  
end  
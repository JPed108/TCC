function V = MatrizV(x, tamanho)
    V(1,1,tamanho) = 1;
    for i=tamanho-1:-1:1
        V(i) = V(i+1)*x;
    end
end
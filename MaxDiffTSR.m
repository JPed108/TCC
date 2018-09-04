function MaxDiffTSR(Coeficientes, x, y, matriztempo)
lin = size(Coeficientes{1}, 1);
col = size(Coeficientes{1}, 2);
ind = sub2ind([col lin], x, y);
x = linspace(0, Coeficientes{end}.Tempo, Coeficientes{end}.nFrames);
if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
    x = log(x+1);
end

if(~Coeficientes{end}.Dividir)
    f = waitbar(1/8, 'Reconstruindo o vídeo térmico...', 'windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();
    frames(end).setAlwaysOnTop(1); 
    p = Coeficientes{1};
    p = permute(p, [2 1 3]);
    p = reshape(p, lin*col, []);  %Polinômios
    clear V;
    t = x';
    V(:, size(p, 2)) = ones(size(t,1),1);
    for j=size(V,2)-1:-1:1
        V(:,j)=V(:, j+1).*t;
    end
    mattmp = (V*p')';
    if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
        mattmp = exp(mattmp);
    end
    waitbar(1/4, f, 'Encontrando os máximos e mínimos...');
    mx = max(mattmp, [], 2);
    mn = min(mattmp, [], 2);
    waitbar(2/4, f, 'Normalizando...');
    mattmp = (mattmp - mn)./(mx-mn);
    matrizref = mattmp;
    waitbar(3/4, f, 'Subtraindo e encontrando a norma...');
    matrizref = abs(matrizref - matrizref(ind, :));
    waitbar(4/4, f, 'Encontrando os máximos das diferenças...');
    matrizref = max(matrizref, [], 2);
    figure; imshow(reshape(matrizref, col, lin)', []);
    close(f);
else
    matriztempo = reshape(matriztempo', lin*col, []);
    if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
        matind = round((exp(matriztempo)-1)*30);  %round por causa de erro de ponto flutuante
    else
        matind = round(matriztempo*30);
    end
    un = unique(matind);
    pAq = permute(Coeficientes{1}, [2 1 3]);
    pRes = permute(Coeficientes{2}, [2 1 3]);
    pAq = reshape(pAq, lin*col, []);  %Polinômios Aquecimento
    pRes = reshape(pRes, lin*col, []); %Polinômios Resfriamento
    imgAq = zeros(lin*col, 1);
    imgRes = zeros(lin*col, 1);
    %%%%%%%%%%%%%%%%%%REFERENCIA%%%%%%%%%%%%%%%%%%%%%
    f = waitbar(1/3, 'Reconstruindo a temperatura de referência...', 'windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();
    frames(end).setAlwaysOnTop(1); 
    matrizref = zeros(1, Coeficientes{end}.nFrames);
    indice = matind(ind);
    clear V;
    t = x(1:indice)';
    V(:, size(pAq, 2)) = ones(indice,1);
    for j=size(V,2)-1:-1:1
        V(:,j)=V(:, j+1).*t;
    end
    mattmp = (V*pAq(ind,:)')';
    if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
        mattmp = exp(mattmp);
    end
    waitbar(2/3, f, 'Encontrando os máximos e mínimos da referência(1)...');
    mx = max(mattmp, [], 2);
    mn = min(mattmp, [], 2);
    waitbar(3/3, f, 'Normalizando a referência(1)...');
    mattmp = (mattmp - mn)./(mx-mn);
    matrizref(1:indice) = mattmp;

    clear V;
    waitbar(1/3, f, 'Reconstruindo a temperatura de referência(2)...');
    t = x(indice+1:end)';
    V(:, size(pRes, 2)) = ones(size(t,1),1);
    for j=size(V,2)-1:-1:1
        V(:,j)=V(:, j+1).*t;
    end
    mattmp = (V*pRes(ind,:)')';
    if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
        mattmp = exp(mattmp);
    end
    waitbar(2/3, f, 'Encontrando os máximos e mínimos da referência(2)...');
    mx = max(mattmp, [], 2);
    mn = min(mattmp, [], 2);
    waitbar(3/3, f, 'Normalizando a referência(2)...');
    mattmp = (mattmp - mn)./(mx-mn);
    matrizref(indice+1:end) = mattmp;

    for i=1:length(un)
        waitbar(i/length(un), f, 'Reconstruindo vídeo térmico(1)...');
        indice = (un(i) == matind);
        matreconst =zeros(sum(indice), Coeficientes{end}.nFrames);
        clear V;
        t = x(1:un(i))';
        V(:, size(pAq, 2)) = ones(un(i),1);
        for j=size(V,2)-1:-1:1
            V(:,j)=V(:, j+1).*t;
        end
        mattmp = (V*pAq(indice,:)')';
        if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
            mattmp = exp(mattmp);
        end
        waitbar(i/length(un), f, 'Encontrando os máximos e mínimos (1)...');
        mx = max(mattmp, [], 2);
        mn = min(mattmp, [], 2);
        waitbar(i/length(un), f, 'Normalizando (1)...');
        mattmp = (mattmp - mn)./(mx-mn);
        matreconst(:, 1:un(i)) = mattmp;
        
        waitbar(i/length(un), f, 'Reconstruindo vídeo térmico(2)...');
        clear V;
        t = x(un(i)+1:end)';
        V(:, size(pRes, 2)) = ones(size(t,1),1);
        for j=size(V,2)-1:-1:1
            V(:,j)=V(:, j+1).*t;
        end
        mattmp = (V*pRes(indice,:)')';
        if(strcmp(Coeficientes{end}.Modelo, 'Polinômio Exponencial'))
            mattmp = exp(mattmp);
        end
        waitbar(i/length(un), f, 'Encontrando os máximos e mínimos (2)...');
        mx = max(mattmp, [], 2);
        mn = min(mattmp, [], 2);
        waitbar(i/length(un), f, 'Normalizando (2)...');
        mattmp = (mattmp - mn)./(mx-mn);
        matreconst(:, un(i)+1:end) = mattmp;
        waitbar(i/length(un), f, 'Encontrando o máximo dos módulos das diferenças(1)...');
        imgAq(indice) = max(abs(matreconst(:, 1:un(i)) - matrizref(:, 1:un(i))), [], 2);
        waitbar(i/length(un), f, 'Encontrando o máximo dos módulos das diferenças(2)...');
        imgRes(indice) = max(abs(matreconst(:, un(i)+1:end) - matrizref(:, un(i)+1:end)), [], 2);
    end
    figure; imshow(reshape(imgAq, col, lin)', []);
    title('Aquecimento');
    figure; imshow(reshape(imgRes, col, lin)', []);
    title('Resfriamento');
end
close(f);

end
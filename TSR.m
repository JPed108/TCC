function [coefs, errmat, mattempo] = TSR(seq, grau, modelo, dividir, filtrar, TamanhoFiltro)  
    Propriedades.Tempo = seq.Duration.TotalSeconds;
    Propriedades.nFrames = seq.Count;
    Propriedades.Dividir = dividir;
    Propriedades.Modelo = modelo;
    Propriedades.GrauDerivada = 0;
    Propriedades.FrameRate = seq.FrameRate;
    seq.First();
    h = seq.ThermalImage.Height;
    w = seq.ThermalImage.Width;
    err = zeros(h*w,1);
    ret = System.Drawing.Rectangle(0,0,w, h);
    matvid = zeros(h*w, seq.Count);
    f = waitbar(0, 'Convertendo vídeo em matriz...');
    indices = round(linspace(1, double(seq.Count), 20));
    for ind=1:seq.Count
        seq.SelectedIndex = ind - 1;
        if(any(ind == indices))
            waitbar(double(ind)/double(seq.Count), f);
        end
        matvid(:,ind) = double(seq.ThermalImage.GetValues(ret))';
        if(strcmp(modelo, 'Polinômio Exponencial'))
            matvid(:,ind) = log(matvid(:,ind));
        end
    end
    close(f);
    if(filtrar)
        f = waitbar(0, 'Aplicando filtro da mediana...');
        x = round(linspace(1,double(w*h),50));
        for i=2:length(x)
            waitbar(double(x(i))/double(w*h), f);
            matvid(x(i-1):x(i), :) = medfilt1(matvid(x(i-1):x(i), :), TamanhoFiltro, [], 2);
        end
        close(f);
    end
    
    if(strcmp(modelo, 'Polinômio'))
        x = linspace(0, seq.Duration.TotalSeconds, seq.Count)';
    else
        x = log(linspace(0, seq.Duration.TotalSeconds, seq.Count)+1)';
    end
    V_geral(:, grau+1) = ones(length(x), 1);
    for i=grau:-1:1
        V_geral(:, i) = V_geral(:, i+1).*x;
    end
    
    if(dividir)
        f = waitbar(0/1, 'Gerando matriz de índices de máximos...');
        coefs = cell(1,2);
        errmat = cell(1,2);
        [~, matind] = max(matvid, [], 2);
        mattempo = matind/double(seq.FrameRate);
        if(strcmp(modelo, 'Polinômio Exponencial'))
            mattempo = log(mattempo+1);
        end
        waitbar(1, f);
        close(f);
    else
        V = V_geral;
        [Q,R] = qr(V,0);
    end
   
    
    if(~dividir)
        f = waitbar(0, 'Gerando coeficientes do TSR...');
        p = (R\(Q'*matvid'))';
        waitbar(1, f);
        close(f);
        indices = round(linspace(1,size(p,1),10));
        f = waitbar(0, 'Calculando r^2...');
        for i=2:length(indices)
            mat = p(indices(i-1):indices(i),:);
            y = (V*mat')';
            med = mean(matvid(indices(i-1):indices(i),:), 2);
            SStot = sum((matvid(indices(i-1):indices(i),:)-med).^2,2);
            SSres = sum((matvid(indices(i-1):indices(i),:)-y).^2,2);
            err(indices(i-1):indices(i)) = 1-SSres./SStot;
            waitbar(i/length(indices), f);
        end
        errmat = {reshape(err, w, h, [])'};
        coefs = {permute(reshape(p, w, h, []), [2 1 3])};
        close(f);
    else
        un = unique(matind);
        for i=1:2
            p = zeros(w*h, grau+1);
            f = waitbar(0, strcat('Gerando coeficientes do TSR e calculando o r^2... (', int2str(i), ')'));
            for j=1:length(un)
                indice = (un(j) == matind);
                if(i==1)
                    mat = matvid(indice, 1:un(j));
                    t = x(1:un(j));
                else
                    mat = matvid(indice, un(j):end);
                    t = x(un(j):end);
                end
                V = [];
                V(:, grau+1) = ones(length(t), 1);
                for k=grau:-1:1
                    V(:, k) = V(:, k+1).*t;
                end
                [Q, R] = qr(V, 0);
                pol = (R\(Q'*mat'))';
                p(indice, :) = pol;
                y = (V*pol')';
                med = mean(mat, 2);
                SStot = sum((mat-med).^2,2);
                SSres = sum((mat-y).^2,2);
                err(indice) = 1-SSres./SStot;
                waitbar(j/length(un), f);
            end
            coefs{i} = permute(reshape(p, w, h, []), [2 1 3]);
            errmat{i} = reshape(err, w, h, [])';
            close(f);
        end
    end
    coefs{end+1} = Propriedades;
    if(dividir)
        mattempo = reshape(mattempo, w, h)';
    end
end

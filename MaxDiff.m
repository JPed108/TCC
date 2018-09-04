function [maxdiff] = MaxDiff(seq, x, y, dividir)
    seq.First();
    h = seq.ThermalImage.Height;
    w = seq.ThermalImage.Width;
    ind = sub2ind([w h], x, y);
    ret = System.Drawing.Rectangle(0,0,w, h);
    matvid = zeros(h*w, seq.Count);
    f = waitbar(0, 'Convertendo vídeo em matriz...');
    indices = round(linspace(1, double(seq.Count), 20));
    for i=1:seq.Count
        seq.SelectedIndex = i - 1;
        if(any(i == indices))
            waitbar(double(i)/double(seq.Count), f);
        end
        matvid(:,i) = double(seq.ThermalImage.GetValues(ret))';
    end
    close(f);
    [mx, indmx] = max(matvid, [], 2);
    if(~dividir)
        if(matvid(i, 1) > matvid(i, end))
            titulo = 'Resfriamento';
        else
            titulo = 'Aquecimento';
        end
        f = waitbar(1/5, 'Encontrando máximos e mínimos...');
        %Encontrando máximo e mínimo
        [mn] = min(matvid, [], 2);
        %Normalizando
        waitbar(2/5, f, 'Normalizando...');
        matvid = (matvid - mn)./(mx-mn);
        %diferença entre (x,y) e todos os pontos
        referencia = matvid(ind, :);
        waitbar(3/5, f, 'Subtraindo do ponto de referencia...');
        matvid = matvid - referencia;
        %O interesse é a maior diferença, então encontra-se o módulo
        waitbar(4/5, f, 'Encontrando o módulo das diferenças...');
        matvid = abs(matvid);
        %Encontra os máximos de cada píxel
        waitbar(5/5, f, 'Encontrando os máximos de cada diferença...');
        maxdiff = max(matvid, [], 2);
        close(f);
        figure; imshow(reshape(maxdiff, w, h)', []);
        title(titulo);
    else
        f = waitbar(1/8, 'Encontrando os mínimos...', 'windowstyle', 'modal');
        frames = java.awt.Frame.getFrames();
        frames(end).setAlwaysOnTop(1); 
        mn = min(matvid(:, 1:indmx-1), [], 2);
        waitbar(2/8, f, 'Normalizando...');
        matvid(:, 1:indmx-1) = (matvid(:, 1:indmx-1) - mn)./(mx-mn);
        waitbar(3/8, f, 'Subtraindo e encontrando os módulos...');
        referencia = matvid(ind, 1:indmx-1);
        matvid(:, 1:indmx-1) = abs(matvid(:, 1:indmx-1) - referencia);
        waitbar(4/8, f, 'Encontrando os máximos das diferenças...');
        maxdiff = max(matvid(:, 1:indmx-1), [], 2);
        figure; imshow(reshape(maxdiff, w, h)', []);
        title('Aquecimento');

        waitbar(5/8, f, 'Encontrando os mínimos...');
        mn = min(matvid(:, indmx:end), [], 2);
        waitbar(6/8, f, 'Normalizando...');
        matvid(:, indmx:end) = (matvid(:, indmx:end) - mn)./(mx-mn);
        referencia = matvid(ind, indmx:end);
        waitbar(7/8, f, 'Subtraindo e encontrando os módulos...');
        matvid(:, indmx:end) = abs(matvid(:, indmx:end) - referencia);
        waitbar(8/8, f, 'Encontrando os máximos das diferenças...');
        maxdiff = max(matvid(:, indmx:end), [], 2);
        figure; imshow(reshape(maxdiff, w, h)', []);
        title('Resfriamento');
        close(f);
    end
end
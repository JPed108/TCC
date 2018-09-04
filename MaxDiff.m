function [maxdiff] = MaxDiff(seq, x, y, dividir)
    seq.First();
    h = seq.ThermalImage.Height;
    w = seq.ThermalImage.Width;
    ind = sub2ind([w h], x, y);
    ret = System.Drawing.Rectangle(0,0,w, h);
    matvid = zeros(h*w, seq.Count);
    f = waitbar(0, 'Convertendo v�deo em matriz...');
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
        f = waitbar(1/5, 'Encontrando m�ximos e m�nimos...');
        %Encontrando m�ximo e m�nimo
        [mn] = min(matvid, [], 2);
        %Normalizando
        waitbar(2/5, f, 'Normalizando...');
        matvid = (matvid - mn)./(mx-mn);
        %diferen�a entre (x,y) e todos os pontos
        referencia = matvid(ind, :);
        waitbar(3/5, f, 'Subtraindo do ponto de referencia...');
        matvid = matvid - referencia;
        %O interesse � a maior diferen�a, ent�o encontra-se o m�dulo
        waitbar(4/5, f, 'Encontrando o m�dulo das diferen�as...');
        matvid = abs(matvid);
        %Encontra os m�ximos de cada p�xel
        waitbar(5/5, f, 'Encontrando os m�ximos de cada diferen�a...');
        maxdiff = max(matvid, [], 2);
        close(f);
        figure; imshow(reshape(maxdiff, w, h)', []);
        title(titulo);
    else
        f = waitbar(1/8, 'Encontrando os m�nimos...', 'windowstyle', 'modal');
        frames = java.awt.Frame.getFrames();
        frames(end).setAlwaysOnTop(1); 
        mn = min(matvid(:, 1:indmx-1), [], 2);
        waitbar(2/8, f, 'Normalizando...');
        matvid(:, 1:indmx-1) = (matvid(:, 1:indmx-1) - mn)./(mx-mn);
        waitbar(3/8, f, 'Subtraindo e encontrando os m�dulos...');
        referencia = matvid(ind, 1:indmx-1);
        matvid(:, 1:indmx-1) = abs(matvid(:, 1:indmx-1) - referencia);
        waitbar(4/8, f, 'Encontrando os m�ximos das diferen�as...');
        maxdiff = max(matvid(:, 1:indmx-1), [], 2);
        figure; imshow(reshape(maxdiff, w, h)', []);
        title('Aquecimento');

        waitbar(5/8, f, 'Encontrando os m�nimos...');
        mn = min(matvid(:, indmx:end), [], 2);
        waitbar(6/8, f, 'Normalizando...');
        matvid(:, indmx:end) = (matvid(:, indmx:end) - mn)./(mx-mn);
        referencia = matvid(ind, indmx:end);
        waitbar(7/8, f, 'Subtraindo e encontrando os m�dulos...');
        matvid(:, indmx:end) = abs(matvid(:, indmx:end) - referencia);
        waitbar(8/8, f, 'Encontrando os m�ximos das diferen�as...');
        maxdiff = max(matvid(:, indmx:end), [], 2);
        figure; imshow(reshape(maxdiff, w, h)', []);
        title('Resfriamento');
        close(f);
    end
end
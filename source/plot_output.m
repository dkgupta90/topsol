newmap = summer;
ncol = size(newmap, 1);
newmap(ncol+1,:) = [1 1 1];

%newmap = jet;
colormap(newmap);
nelx = 599;
nely = 599;
opath = 'output/output_penal_penalS/Output_';
%opath = 'output/output_maple_leaf_cluster/Output_';
opath = [opath, num2str(nelx), '_', num2str(nely), '_'];
figure(1);
% figure(2);
% figure(3);
opath1 = opath;
for i = 3:0.5:3
    for j = 3:0.5:3
        opath = [opath1, 'penal', num2str(i), '_penalS', num2str(j), '/'];
        if exist([opath, 'density.dat'])
        else
            continue;
        end
        load([opath, 'density.dat']);
        load([opath, 'current.dat']);
        load([opath, 'voltage.dat']);
        load([opath, 'workspace.mat']);
        %density(density > 0.49 & density < 0.51) = 3;
        figure(1)
        colormap(newmap);
        density(density == -1) = -0.1;
        imagesc(1 - density);
        hold on;
        str_title = (['Eff: ', num2str(Eff), '; p: ', num2str(i), '; r = ', num2str(j)]); 
        %title(str_title);
        
        hold off;
        axis off;
        colorbar;
        caxis([0, 1.1]);
%         figure(2)
%         imagesc(voltage);
%         figure(3)
%         imagesc(current);
    end
end
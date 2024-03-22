function TopOptVideo(ro,FEM,Rate,Animation,idx,ObjFun,Freq,ini)
% Video creation

% recover Symmetry
symm_mesh.conect = FEM.mesh.conect(:,1:8);
symm_mesh.conect(FEM.mesh.nel+1:2*FEM.mesh.nel,:) = ...
    FEM.mesh.conect(:,1:8)+FEM.mesh.nnod;
symm_mesh.nel = 2*FEM.mesh.nel;
symm_mesh.ncoord = FEM.mesh.ncoord;
symm_mesh.ncoord(FEM.mesh.nnod+1:2*FEM.mesh.nnod,:) = FEM.mesh.ncoord;
symm_axis = 2; % 1-x, 2-y, 3-z
symm_coordinate = max(FEM.mesh.ncoord(:,symm_axis));
symm_mesh.ncoord(FEM.mesh.nnod+1:2*FEM.mesh.nnod,symm_axis) = ...
    -FEM.mesh.ncoord(:,symm_axis) + 2*symm_coordinate;
ro(FEM.mesh.nel+1:2*FEM.mesh.nel,:) = ro;

% create the video writer
v = VideoWriter('TopOpt.avi');
% set the seconds per image
v.FrameRate = Rate;

% open the video writer
open(v);

switch Animation
case 1
figure('position', [0, 0, 800, 600]); ax = gca; ax.FontSize = 12;
for it = 1:size(ro,2)
    PlotTopology(symm_mesh,ro(:,it),'k','w');
    title(['Iteration: ' num2str(it)])
    frame = getframe(gcf);
    writeVideo(v, frame);
    pause(eps); ax = gca; cla(ax,'reset');
end

case 2 % Plot only eingefrequencies
figure('position', [0, 0, 1500, 600]); ax = gca; ax.FontSize = 18;
for loop = 1:length(idx)
    it = idx(loop);
    subplot(1,2,1);   ax = gca; cla(ax,'reset');
    PlotTopology(symm_mesh,ro(:,it),'k','w');
    title(['Iteration: ' num2str(it)],'FontSize',18)
    subplot(1,2,2); ax = gca; cla(ax,'reset');
    semilogy(Freq/1e3,ObjFun.FRF(:,it),'-k','LineWidth',1.5);
    ylim([min(min(ObjFun.FRF)) max(max(ObjFun.FRF))]);
    xlabel('Frequency (kHz)','FontSize',18)
    xline(ini.Target/1e3,'-k');
    if any(ObjFun.ResEigF(it,:))
        xline(nonzeros(ObjFun.ResEigF(it,:))/1e3,'--r','LineWidth',1.5);
    end
    frame = getframe(gcf); writeVideo(v, frame);
    pause(eps);
end

case 3 % Plot only Antieigenfreq
figure('position', [0, 0, 1500, 600]); ax = gca; ax.FontSize = 18;
for loop = 1:length(idx)
    it = idx(loop);
    subplot(1,2,1);   ax = gca; cla(ax,'reset');
    PlotTopology(symm_mesh,ro(:,it),'k','w'); 
    title(['Iteration: ' num2str(it)])
    subplot(1,2,2); ax = gca; cla(ax,'reset');
    semilogy(Freq/1e3,ObjFun.FRF(:,it),'-k','LineWidth',2);
    ylim([min(min(ObjFun.FRF)) max(max(ObjFun.FRF))]);
    xlabel('Frequency (kHz)','FontSize',18);
    xline(ini.Target/1e3,'-k');
    xline(nonzeros(ObjFun.AntiEigF(it)/1e3),'--b','LineWidth',1.5);
    xline(ObjFun.RefAntiFreq(it)/1e3,'--r');
    frame = getframe(gcf); writeVideo(v, frame);
    pause(eps);
end

case 4 % Plot both Resonances and Antieigenfreq
figure('position', [0, 0, 1500, 600]); ax = gca; ax.FontSize = 18;
for loop = 1:length(idx)
    it = idx(loop);
    subplot(1,2,1);   ax = gca; cla(ax,'reset');
    PlotTopology(symm_mesh,ro(:,it),'k','w'); 
    title(['Iteration: ' num2str(it)],'FontSize',18)

    subplot(1,2,2); ax = gca; cla(ax,'reset');
    semilogy(Freq/1e3,ObjFun.FRF(:,it),'-k','LineWidth',2);
    ylim([min(min(ObjFun.FRF)) max(max(ObjFun.FRF))]);
    xlim([min(Freq/1e3) max(Freq/1e3)])
    xlabel('Frequency (kHz)','FontSize',18);
    ax.FontSize = 18;
    ResEigF = nonzeros(ObjFun.ResEigF(it,:));
    AntiEigF = nonzeros(ObjFun.AntiEigF(it,:)); 
    xline(ini.Target/1e3,'--b','LineWidth',0.5);
    if any(ResEigF); xline(ResEigF/1e3,'--r','LineWidth',1.2); end
    xline(AntiEigF/1e3,'--b','LineWidth',1.2);
    frame = getframe(gcf); writeVideo(v, frame);
    pause(eps);
end

% end of switch
end

% close the writer object
close(v);
end
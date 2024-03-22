%% latex plots
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Clean data
idc = ~any(ObjFun.ObjFunVal,2);
ObjFun.ObjFunVal(idc,:) = [];
ObjFun.Szz(idc,:) = [];
ObjFun.Ux(idc,:) = [];
ObjFun.Uy(idc,:) = [];
ObjFun.Uz(idc,:) = [];
ObjFun.FRF1(:,idc) = [];
ObjFun.FRF2(:,idc) = [];
ObjFun.FRF3(:,idc) = [];
ObjFun.AntiEigF1(idc,:) = [];
ObjFun.AntiEigF2(idc,:) = [];
ObjFun.AntiEigF3(idc,:) = [];
ObjFun.RefAntiFreq1(idc,:) = [];
ObjFun.RefAntiFreq2(idc,:) = [];
ObjFun.RefAntiFreq3(idc,:) = [];
ObjFun.ResEigF1(idc,:) = [];
% ObjFun.ResEigF2(idc,:) = [];
ObjFun.ResEigF3(idc,:) = [];
% ObjFun.dFun(:,idc) = [];
rho(:,idc) = []; rho_h(:,idc) = [];
idx = 1:length(ObjFun.ObjFunVal);
clear idc; %save('TopOpt07')

%% Objective Function
figure('position', [0, 0, 800, 500])
plot(ObjFun.ObjFunVal,'-k','LineWidth',2); grid
xlabel('Iteration'); ylabel('Objective Function Value');
legend('Objective function')
ax = gca; ax.FontSize = 16;
% [minObj,idx] = min(ObjFun.ObjFunVal);

hold on
plot(ini.w1*ObjFun.Szz);
plot(ini.w2*ObjFun.Ux);
plot(ini.w3*ObjFun.Uy);
% plot(ini.w4*ObjFun.Uz);
% plot(ini.w5.*( (ObjFun.AntiEigF1 -ini.Target1)/ini.Target1 ).^2);
% plot(ini.w6.*( (ObjFun.AntiEigF2 -ini.Target2)/ini.Target2 ).^2);
% plot(ini.w7.*( (ObjFun.AntiEigF3 -ini.Target3)/ini.Target3 ).^2);

%% Volume evolution
TotalVolume  = sum(FEM.mesh.Ve); VolumenPerIt = (FEM.mesh.Ve')*rho;
VolumePercentage = (VolumenPerIt/TotalVolume)*100;
figure('position', [0, 0, 800, 500]);
plot(VolumePercentage,'-b','LineWidth',1.5); grid;
yline(ini.MaxVol*100,'--r'); yline(ini.MinVol*100,'--r')
ylim([0 ini.MaxVol*100+5]); % xlim([1 ini.MaxIter])
xlabel('Iteration'); ylabel('Volume percentage')
ax = gca; ax.FontSize = 16;

%% Topology animation
figure('position', [0, 0, 1500, 600]);
for loop = 1:length(idx)
    it = idx(loop);
    subplot(1,2,1); ax = gca; cla(ax,'reset'); grid on
    ro = rho_h(:,it); ro = MaterialModel(ini.MatModel,ro,ini.p);
    PlotTopology(FEM.mesh,ro,'k','w'); title(['Iteration: ' num2str(it)])
    
    subplot(1,2,2);
    if any(Freq)
    semilogy(Freq/1e3,ObjFun.FRF1(:,it),'-g','LineWidth',2); hold on
    semilogy(Freq/1e3,ObjFun.FRF3(:,it),'-b','LineWidth',2); hold off
    end
%     ylim([ min([min(min(ObjFun.FRF1)),min(min(ObjFun.FRF3))]) ...
%         max([max(max(ObjFun.FRF3)),max(max(ObjFun.FRF3))]) ]);
    if any(ObjFun.Ux)
        yline(ObjFun.Ux(it),'--k');
    end
    if any(ObjFun.Uz)
        yline(ObjFun.Uz(it),'--k');
    end
    if any(ObjFun.AntiEigF1(it))
        xline(ObjFun.RefAntiFreq1(it)/1e3,'--k');
        xline(nonzeros(ObjFun.AntiEigF1(it))/1e3,'--k');
    end
    if any(ObjFun.AntiEigF3(it))
        xline(ObjFun.RefAntiFreq3(it)/1e3,'--c');
        xline(nonzeros(ObjFun.AntiEigF3(it))/1e3,'--b');
    end
    if any(ObjFun.ResEigF1(it,:))
        xline(nonzeros(ObjFun.ResEigF1(it,:))/1e3,'-g');
    end
    if any(ObjFun.ResEigF3(it,:))
        xline(nonzeros(ObjFun.ResEigF3(it,:))/1e3,'-b');
    end
    xlabel('Frequency (kHz)');
    pause(eps);
end; grid on; ax = gca; ax.FontSize = 14;

%% Generate data for STL
ro = MaterialModel(ini.MatModel,rho_h(:,it),ini.p);
xPhys = TOPslicerData(FEM.mesh.Ce,ro); % Export data for TOPslicer
save('xPhys','xPhys'); TOPslicer

%% Objective funcion values
F1 = ini.w1*((ObjFun.AntiEigF(idx)-ini.Target)/ini.Target).^2;
F2 = sum(ini.w2*(sqrt((ObjFun.AntiEigF(idx)./...
    (ObjFun.ResEigF(idx,:)-ObjFun.AntiEigF(idx))).^2)) );

%% Plot Topology
figure('position', [0, 0, 800, 600])
ro = MaterialModel(ini.MatModel,rho_h(:,it),ini.p);
PlotTopology(FEM.mesh,ro,'k','w');
title(['Iteration = ' num2str(it)])

%% Remove disconnected memebers manually
Remove_Elements =  FEM.mesh.Ce(:,3) > 0.012;
rho_h(Remove_Elements,it) = 0;
Remove_Elements = FEM.mesh.Ce(:,2) < 0.002;
rho_h(Remove_Elements,it) = 0;

%% Topology with Symmetry recovered
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

figure('position', [0, 0, 800, 600])
ro = MaterialModel(ini.MatModel,rho_h(:,it),ini.p);
ro(FEM.mesh.nel+1:2*FEM.mesh.nel,:) = ro;
PlotTopology(symm_mesh,ro,'k','w'); grid

%% Plot FRF
figure('position', [0, 0, 1000, 600])
semilogy(Freq/1e3,ObjFun.FRF1(:,it),'-g','LineWidth',2); hold on
semilogy(Freq/1e3,ObjFun.FRF3(:,it),'-b','LineWidth',2); hold off
if any(ObjFun.AntiEigF1(it))
    xline(ObjFun.RefAntiFreq1(it)/1e3,'--g');
    xline(nonzeros(ObjFun.AntiEigF1(it))/1e3,'--k');
end
if any(ObjFun.AntiEigF3(it))
    xline(ObjFun.RefAntiFreq3(it)/1e3,'--c');
    xline(nonzeros(ObjFun.AntiEigF3(it))/1e3,'--b');
end
if any(ObjFun.ResEigF1(it,:))
    xline(nonzeros(ObjFun.ResEigF1(it,:))/1e3,'-k');
end
xlabel('Frequency (kHz)'); grid
xlim([Freq(1)/1e3 Freq(end)/1e3])
legend('abs(u)','abs(w)')
ax = gca; ax.FontSize = 16;

%% Plot density distrbution
ro = MaterialModel(ini.MatModel,rho_h(:,it),ini.p);
figure('position', [0, 0, 800, 600]);
plot(rho(:,it),'r','LineStyle','none','Marker','*'); grid; ylim([0 1])
hold on; plot(rho_h(:,it),'g','LineStyle','none','Marker','*')
plot(ro,'b','LineStyle','none','Marker','*')
ax = gca; ax.FontSize = 16; xlabel('Element number'); ylabel('Pseudo-density')

%% Topology Video
ro = MaterialModel(ini.MatModel,rho_h,ini.p);
FrameRate = 2; % Frames per second
Animation = 4;
TopOptVideo(ro,FEM,FrameRate,Animation,idx,ObjFun,Freq,ini)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Modal Eigenvectors solution
Mode = 1; E = EigV_h; % E = HarmDisp_h;
figure('position', [0, 0, 800, 600],'Color',[1 1 1])
PlotFEA(FEM.mesh,E,6,Mode); colorbar('off')

%% Generate interpolation models plots
ro = 0:0.01:1;
[roK,roM] = MaterialModel('SIMP',ro,3,1);
figure; plot(ro,roM,ro,roK); grid; legend('Mass','Stiffness')

%% Plot Heaviside curve
ro = 0:0.005:1;
eta = 0.5;
beta = 30;
rh(ro<=eta) = eta*( exp(-beta*(1-ro(ro<=eta)/eta)) - ...
              (1-ro(ro<=eta)/eta)*exp(-beta) );
rh(ro>eta) = (1-eta)*( 1 - exp( -beta*(ro(ro>eta)-eta)/(1-eta) ) ...
            + ( (ro(ro>eta)-eta)*exp(-beta) )/(1-eta) ) + eta;
figure; plot(ro,rh,'-b','LineWidth',2); grid;
xlabel('Input densitiy'); ylabel('Output density')
ax = gca; ax.FontSize = 16;
yline(eta,'--m'); %yline(0.7,'--k'); yline(0.5,'--k'); yline(0.1,'--k')
legend(['$\beta =$' num2str(beta)],'$\eta$')

%% Get antiresonance quality metrics
figure('position', [0, 0, 800, 500])
Qfactor = zeros(length(idx),1);
for loop = 1:length(idx)
    it = idx(loop);
    Uavg = ObjFun.FRF1(:,it);
    AntiUavg = 1./Uavg;
    [Peak,idP,Width,Prom] = findpeaks(AntiUavg);
    Prox = abs(ini.Target1 - Freq(idP));
    % Normalize data
    PeakN = Peak/max(Peak); 
    WidthN = Width/max(Width);
    PromN = Prom/max(Prom);
    ProxN = Prox/max(Prox);
    % Evaluate peak metrics
    factor = PeakN +PromN +WidthN -ProxN; [~,idf] = max(factor);
    % display( Freq(idP(idf)) )
    % display(Peak(idx))
    % display(Prox(idx))
    semilogy(Freq,Uavg); grid; pause(eps)
    
    % Q factor
    PeakRMS = Peak(idf)/sqrt(2);
    
    PrevAmp = Peak(idf)-PeakRMS;
    for pp = idP(idf):-1:1
        AmpDiff = abs( AntiUavg(pp)-PeakRMS );
        if AmpDiff > PrevAmp && AntiUavg(pp) < PeakRMS
            idc = pp+1;
            break
        else
            PrevAmp = AmpDiff;
        end
    end
    AmpDiff = AntiUavg(idc)-PeakRMS;
    if AmpDiff < 0
        f_l = interp1([AntiUavg(idc) AntiUavg(idc+1)],...
            [Freq(idc) Freq(idc+1)],PeakRMS);
    elseif AmpDiff >= 0
        f_l = interp1([AntiUavg(idc) AntiUavg(idc-1)],...
            [Freq(idc) Freq(idc-1)],PeakRMS);
    end
    
    PrevAmp = Peak(idf)-PeakRMS;
    for pp = idP(idf):1:length(Freq)
        AmpDiff = abs( AntiUavg(pp)-PeakRMS );
        if AmpDiff > PrevAmp && AntiUavg(pp) < PeakRMS
            idc = pp-1;
            break
        else
            PrevAmp = AmpDiff;
        end
    end
    AmpDiff = AntiUavg(idc)-PeakRMS;
    if AmpDiff >= 0
        f_h = interp1([AntiUavg(idc) AntiUavg(idc+1)],...
            [Freq(idc) Freq(idc+1)],PeakRMS);
    elseif AmpDiff < 0
        f_h = interp1([AntiUavg(idc) AntiUavg(idc-1)],...
            [Freq(idc) Freq(idc-1)],PeakRMS);
    end
    % HalfPower = 100*(f_h-f_l)/Freq(idP(idx)); % half-power bandwidth
    Qfactor(it) = Freq(idP(idf))/(f_h-f_l);
end
[B,I] = sort(Qfactor);
idQ = I(I>50);

%% Plot FEA only
FRF = FRF1; % FRF1, FRF2, or FRF3
semilogy(Freq/1e3,FRF,'-b','LineWidth',2); hold off
xlabel('Frequency (kHz)'); grid
xlim([Freq(1)/1e3 Freq(end)/1e3])
ax = gca; ax.FontSize = 16;
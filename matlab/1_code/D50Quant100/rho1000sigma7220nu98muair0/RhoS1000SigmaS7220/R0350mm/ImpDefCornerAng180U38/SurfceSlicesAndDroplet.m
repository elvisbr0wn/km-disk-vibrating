clc
clear
close all

load('U0.mat')
load('Ang.mat'); Cang = Ang * pi / 180;
load('vz.mat'); Vo = abs(vz(1));
files = dir(fullfile(pwd, "etaMatPer*.mat"));
N = length(files);
etaAux = [];
for i = 1:N
    load(files(i).name);
    etaAux = [etaAux, etaMatPer];
end
etaMatPer = etaAux;
load('z.mat')
load('etaOri.mat')
load('tvec.mat')
load('Fr.mat')

%%

cd ..
load('Ro.mat','Ro')%Sphere's radius in CGS

cd ..
load('rhoS.mat','rhoS')%Sphere density
%load('sigmaS.mat')%Sphere's surface tension

cd ..
load('rho.mat','rho')
%load('sigma.mat','sigma')
load('nu.mat','nu')
load('muair.mat')
load('g.mat','g') %gravitational constant

cd ..
%load('D.mat')%Domain diameter in units of droplet radii
%load('quant.mat')%number of dr's contained in an undeformed dropelt radius
load('nr.mat','nr')
load('dr.mat','dr')
load('r.mat')
load('zs.mat','zs')
%load('Deload('zs.mat','zs')lta.mat','Delta')
%load('IntMat.mat','IntMat')
%load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, D),'DTNnew345')
%DTN = DTNnew345;
%clear DTNnew345

%xplot = dr*(0:nr-1); save('xplot.mat','xplot')%I might remove or relocate this
load('xplot.mat')

%cd(['rho',num2str(1000*rho),'sigma',num2str(round(100*sigma)),'nu',num2str(round(10000*nu)),'muair',num2str(muair)])

%cd(['RhoS',num2str(rhoS*1000),'SigmaS',num2str(round(100*sigmaS))])
%load('Ma.mat','Ma')%Dimensionless mass of sphere
%load('Ra.mat','Ra')%Density ratio

%%
% cd ..
% load('Ro.mat');
% 
% cd ..
% load('rhoS.mat');

%cd ..
%load('r.mat')
%load('nr.mat')
%load('dr.mat')
%load('zs.mat','zs')
%load('Gamma.mat')
%load('wzero.mat')
%load('thetaZero.mat')
%load('xplot.mat')

%cd(['Rho',num2str(rhoS*1000)])

%cd(['ImpDefAng',num2str(round((pi-Cang)*180/pi)),'U',num2str(U0)])
%load('Rv.mat')
Rv = ones(size(etaMatPer, 2), 1);
%load('Fr.mat')


%Vo = abs(vz(1));

% load('etaMatPer1.mat')
% etaMatPer1 = etaMatPer;
% load('etaMatPer2.mat')
% etaMatPer2 = etaMatPer;
% load('etaMatPer3.mat')
% etaMatPer3 = etaMatPer;
% load('etaMatPer4.mat')
% etaMatPer4 = etaMatPer;
% load('etaMatPer5.mat')
% etaMatPer5 = etaMatPer;
% load('etaMatPer6.mat')
% etaMatPer6 = etaMatPer;
% load('etaMatPer7.mat')
% etaMatPer7 = etaMatPer;
% load('etaMatPer8.mat')
% etaMatPer8 = etaMatPer;
% load('etaMatPer9.mat')
% etaMatPer9 = etaMatPer;
% load('etaMatPer10.mat')
% etaMatPer10 = etaMatPer;
% load('etaMatPer11.mat')
% etaMatPer11 = etaMatPer;
% load('etaMatPer12.mat')
% etaMatPer12 = etaMatPer;
% load('etaMatPer13.mat')
% etaMatPer13 = etaMatPer;
% load('etaMatPer14.mat')
% etaMatPer14 = etaMatPer;
% load('etaMatPer15.mat')
% etaMatPer15 = etaMatPer;
% load('etaMatPer16.mat')
% etaMatPer16 = etaMatPer;
% 
% etaMatPer = [etaMatPer1,etaMatPer2,etaMatPer3,etaMatPer4,etaMatPer5,etaMatPer6,...
%     etaMatPer7,etaMatPer8,etaMatPer9,etaMatPer10,etaMatPer11,etaMatPer12,etaMatPer13,...
%     etaMatPer14,etaMatPer15,etaMatPer16];

ntimes = size(etaMatPer,2);

zplot = z(1:end-1);
etaOriplot = etaOri(1:end-1);

zmin = min(min(etaMatPer));
zmax = max(max(etaMatPer));


vidObj = VideoWriter('WavesAndDrop.mp4','MPEG-4');
set(vidObj,'FrameRate',10)
open(vidObj);

Gamma = 0; %%
wzero = 1; %%
thetaZero = 0; %%
zb = Gamma/(Fr*wzero^2)*cos(wzero*tvec+thetaZero); %Elevation of the pool
zbplot=zb(1:end-1);

for ii = 1:size(etaMatPer,2)
    RvCurr = 1; %Rv(ii); %% TEMP CHANGE
    Rh = sqrt(1/RvCurr);
    nlmax = floor(Rh/dr)+1;
    zs(1:nlmax) = RvCurr-RvCurr*sqrt(1-RvCurr*(0:dr:(nlmax-1)*dr).^2);
    zsplot = [(zplot(ii)-RvCurr)+[flipud(zs(2:nlmax));zs(1:nlmax)];flipud((z(ii)+RvCurr)...
        -[flipud(zs(2:nlmax));zs(1:nlmax)]);(z(ii)-RvCurr)+zs(nlmax)];
    xs = [xplot(nr-nlmax+2:nr+nlmax),fliplr(xplot(nr-nlmax+2:nr+nlmax)),xplot(nr-nlmax+2)];
    plot(xs,(zbplot(ii)+zsplot),'k','LineWidth',2)
    hold on
    axis equal

    plot([fliplr(-1*r),r(2:end)],(zbplot(ii)+[0;flipud(etaMatPer(:,ii));etaMatPer(2:end,ii);0]),...
        'color',[.4 .4 .4],'LineWidth',2)
    grid on
    set(gca,'xlim',[-5 5],'ylim',[-1.5 3],'Xtick',-5:5,'FontName','Times','FontSize',24);
    xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',24)
    ylabel('   $\frac{z}{R_o}\ \ \ $   ','interpreter','Latex','FontName','Times',...
        'FontSize',24,'rotation',0)
    t = (ii-1)/360;
    to = floor(t);
    t1 = floor(10*(t-to));
    t2 = round(100*(t-to-t1/10));
    if t2 == 10
        t2=0;
        t1=t1+1;
    end
    if t1 == 10
        t1 = 0;
        to = to+1;
    end
	ti = round(tvec(ii)*100);
    title(['$   tV_0/R_o =\ $',num2str(floor(ti/100)),'.',num2str(floor(mod(ti,100)/10)),num2str(mod(ti,10))],'FontSize',24,...
            'interpreter','latex','FontName','Times')    
    drawnow
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    hold off
end
close(vidObj);







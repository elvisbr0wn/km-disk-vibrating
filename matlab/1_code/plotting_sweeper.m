% This script will plot a number of variables relating such as coefficient
% of restitution, contact time and max deflection.

%D = 5;
Quant = [100, 200];
%rho = 1; % must multiply by x1000
%sigma = 72.20; % must multiply by x100
%nu = 9.78E-3; % Multiply by x10000
%muair = 0;
%RhoS = 1; % must multiply by x1000
%SigmaS = 72.20; % must multiply by x100
R = [0.035, 0.035]; % linspace(0.02, 0.05, 5)'; % must multiply by x10
%Ang = 180;
%U = 18; %linspace(28, 50, 5)';

data = {  };
files = dir("**/simulation_postprocessing.mat");
We = []; Bo = []; Oh = []; max_deflection = []; contact_time = []; 
coef_restitution = [];
plotting_data = table(We, Bo, Oh, max_deflection, contact_time, coef_restitution);
for ii = 1:length(files)
    
    if is_valid(files(ii).folder, data)
        
        if length(files) == 1; folder_name = files.folder; else folder_name = files(ii).folder; end
        cd(folder_name);
        load("U0.mat");
        cd ..
        load('Ro.mat','Ro')%Sphere's radius in CGS
        cd ..
        %load('rhoS.mat','rhoS')%Sphere density
        %load('sigmaS.mat')%Sphere's surface tension
        cd ..
        load('rho.mat','rho')
        load('sigma.mat','sigma')
        load('nu.mat','nu')
        load('muair.mat')
        load('g.mat','g') %gravitational constant
        cd(folder_name);
        We = rho * U0.^2 * Ro / sigma;
        Bo = rho * g * Ro.^2 / sigma;
        Oh = nu / sqrt(sigma * Ro * rho);

        load("simulation_postprocessing.mat");        
        max_deflection = abs(max_def); if isempty(max_deflection) == true; max_deflection = NaN; end
        coef_restitution = CRref; if isempty(coef_restitution) == true; coef_restitution = NaN; end
        contact_time = tcont; if isempty(contact_time) == true; contact_time = NaN; end

        plotting_data = [plotting_data; {We, Bo, Oh, max_deflection, contact_time, coef_restitution}];
   
    end
end

% Maximum deflection pot
Max_Def_Plot = scatter(plotting_data.We,plotting_data.coef_restitution,'MarkerEdgeColor',[ 0.4660    0.6740    0.1880],'LineWidth',4);
grid on;
%Center = plot(tvec(1:length(z)),z,'k','LineWidth',4);
set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
xlabel('   $We$   ','interpreter','LaTeX','FontSize',20); xlim([0, 2]);
ylabel('   $z / R_o \ \ \ $    ','interpreter','LaTeX','FontSize',20,'Rotation',0); ylim([0, 0.6]);

title(sprintf("Maximum Deflection for \n Bo = %.2g, Oh = %.2e", plotting_data.Bo(1), plotting_data.Oh(1)),'interpreter','LaTeX','FontSize',20);

function bool = is_valid(folder, data)
    bool = true;
end
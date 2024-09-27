function solve_motion(NameValueArgs)

arguments
    % Unless stated otherwise, all units are cgs. 
    % Default values amount to a water bath
    NameValueArgs.diskRadius (1, 1) double = 1 % Radius of the oscillating disk, in cm
    NameValueArgs.diskMass (1, 1) double = 5 % Mass in grams of the disk
    NameValueArgs.forceAmplitude (1, 1) double = 2000 % Amplitude of sinusoidal force applied to disk (in dynes)
    NameValueArgs.forceFrequency (1, 1) double = 20 % Frequency of sinusoidal force in Hz
    NameValueArgs.bathDensity (1, 1) double = 1 % Density of bath's fluid in g/cm^3
    NameValueArgs.bathSurfaceTension (1, 1) double = 72.20 % For water, in dynes/cm
    NameValueArgs.bathViscosity (1, 1) double = 0.978e-2 % Viscosity in Stokes (cgs)
    NameValueArgs.g (1, 1) double = 981 % Gravitational constant in cgs
    NameValueArgs.bathDiameter (1, 1) double = 100 % Diameter of the bath wrt to disk Radius
    NameValueArgs.spatialResolution (1, 1) double = 50 % Number of numerical radial intervals in one disk radius
    NameValueArgs.temporalResolution (1, 1) double = 20; % Number of temporal steps in one adimensional unit
    NameValueArgs.simulationTime (1, 1) double = 1; % Time to be simulated in seconds
    NameValueArgs.debug_flag (1, 1) logical = true; % To show some debugging info
end

datetimeMarker = datetime('now'); datetimeMarker.Format = 'yyyyMMddmmss';
NameValueArgs.datetimeMarker = datetimeMarker;
% Beware! I'm adding all these names into the current scope
cellfun(@(f) assignin('caller', f, NameValueArgs.(f)), fieldnames(NameValueArgs));
%Reset warning
lastwarn('', '');

close all

%tstart = tic;
%data in cgs

currfold = pwd;
addpath(currfold);
precomputedInverse = nan; 
cd ..
fold = fullfile(pwd, sprintf("D%dQuant%d", spatialResolution, bathDiameter));
try
    cd(fold)
    nr = ceil(spatialResolution*bathDiameter/2);
    load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, bathDiameter),'DTNnew345')
    DTN = DTNnew345;
    clear DTNnew345
    % Loading precomputed inverse if exists
    myfile = fullfile(pwd, sprintf("dtstep=%d.mat", forceFrequency*temporalResolution)); 
    if exist(myfile, "file"); load(myfile, "precomputedInverse"); end
    cd(currfold)
catch ME
    error("Could not load DTN for D=%d, Quant=%d. Please generate the matrix first", ...
        spatialResolution, bathDiameter);
end



% Loading some useful matrices
[dr, laplacian, pressureIntegral] = domainMaker(bathDiameter, spatialResolution);

% Dimensional parameters
%forceFrequency = 100; % Oscillations of theforce in Hertz.
%forceAmplitude = 1; % Force amplitude in cgs
%diskMass = 10; % Mass in grams of the object
%bathViscosity = 1;
%Characteristic Units

L_unit = diskRadius; 
M_unit = bathDensity * L_unit^3; % Mass unit. 
%T = sqrt(rhoS * Ro^3/sigmaS); % Characteristic time
T_unit = 1/forceFrequency;
V_unit = L_unit/T_unit;
F_unit = M_unit * L_unit/T_unit^2;

UNITS = struct('length', L_unit, 'mass', M_unit, 'time', T_unit, ...
    'velocity', V_unit, 'force', F_unit);
%Dimensionless numbers for equations
%Dr = rhoS/rho; %Sr = sigmaS/sigma;
Re = L_unit^2/(bathViscosity*T_unit);
Fr = L_unit/(g * T_unit^2) * inf; % Turning off gravity
We = bathDensity * L_unit.^3 / (bathSurfaceTension * T_unit^2); 

force_adim = forceAmplitude/diskMass * T_unit^2/L_unit;
freq_adim  = forceFrequency * T_unit;
obj_mass_adim = diskMass/M_unit;

%WeS  = rhoS*Ro^3/(sigmaS * T_unit^2); %This is for the bath/dropplet interaction.
%Westar = rhoS * U0.^2 * Ro / sigmaS; % velocity-based weber number (to compare with literature)
%Oh = nu*sqrt(rhoS/(sigmaS*Ro));

%Cang = (Ang/180)*pi; %contact angle to be imposed

%Physical parameters
%simulationTime = 20; %Earliest possible end of simulation in characteristic units


%Numerical Simulation parameters

dt = T_unit/temporalResolution; %
steps = ceil(simulationTime/dt); %estimated minimum number of timesteps
if steps * nr * 8 > 1e+9; warning('Spatial resolution and simulation times might be too big to store all matrices in memory'); end
%Inintial conditions for the fluid
%Zeroing result storing variables
%etaOri = zeros(1,steps+1);%height of the surface below the south pole
%obj_height = zeros(1,steps);%height of the centre of mass
%objVelocity = zeros(1,steps);%speed of the centre of mass
%numl = zeros(1,steps+1);%number of pressed mesh points at each time step
%saved_times = zeros(1, steps); %vector of times assuming no refinement has happened

%dt = tvec(2) - tvec(1); indexes_to_save = zeros(steps + 1, 1);
%current_to_save = 2; indexes_to_save(1) = 1;

%nlmax = zeros(1,steps+1);%Variable to store the number of nodes spanned by the deformed droplet
%zeroing variable that records each part of the sequence of surface states
%etaMatPer = zeros(length(etao),nsteps); 
etaInitial = zeros(nr,1); %initial surface elevation
phiInitial = zeros(nr,1); %initial surface potential


current_conditions = struct( ...
    "dt", dt(1), "time", 0, ...
    "center_of_mass", 0, "center_of_mass_velocity", 0, ...
    "bath_surface", etaInitial, "bath_potential", phiInitial, "pressure", zeros(spatialResolution+1, 1));
current_index = 1; %iteration counter
recordedConditions = cell(steps, 1);
recordedConditions{current_index} = current_conditions;


PROBLEM_CONSTANTS = struct("froude", Fr, "weber", We, ...
    "reynolds", Re, "dr", dr, "DEBUG_FLAG", debug_flag, ...
    "nr", nr, "contact_points", spatialResolution+1, ... 
    "force_amplitude", force_adim, "force_frequency", freq_adim, ...
    "DTN", DTN, "laplacian", laplacian, "obj_mass", obj_mass_adim, ...
    "pressure_integral", pressureIntegral(spatialResolution+1, :), ...
    "precomputedInverse", precomputedInverse);

fprintf("Starting simulation on %s\n", pwd);


% Names of the variables to be stored
savingvarNames = { ...
    getVarName(NameValueArgs), ...
    getVarName(PROBLEM_CONSTANTS), ...
    getVarName(recordedConditions), ...
    getVarName(UNITS) ...
};

variableValues = cell(size(savingvarNames));


%% Main Loop
try
    while (recordedConditions{current_index}.time * T_unit < simulationTime) 

        [recordedConditions{current_index+1}, PROBLEM_CONSTANTS] = ...
               advance_one_step(recordedConditions{current_index}, ...
                       PROBLEM_CONSTANTS);
        current_index = current_index + 1;

        if PROBLEM_CONSTANTS.DEBUG_FLAG == true
            width = 1.5; 
            width = min(nr, ceil(width*spatialResolution));
            eta = recordedConditions{current_index}.bath_surface;
            eta = [flipud(eta(2:width)); eta(1:width)];
            xplot = dr*(0:nr-1);
            plot([-fliplr(xplot(2:width)), xplot(1:width)], eta, 'b', 'Linewidth', 2);
            hold on
            
            % Filling the shape of the vibrating object
            x = [1, 1, -1, -1];
            z = recordedConditions{current_index}.center_of_mass;
            y = [z, z+1/10, z+1/10, z];
            fill(x, y, 'k');
            % xs = dr*(0:nlmax(tentative_index+1)-1);
            % zsplot = zs(1:nlmax(tentative_index+1))+RvTent+obj_height(tentative_index+1);
            % plot([-fliplr(xs(2:end)),xs],[flipud(zsplot(2:end));zsplot],'k','Linewidth',2);
            % hold on
            % thetaplot = linspace(0, thetaVec(end), 200);%-%-0:thetaVec(end)/400:thetaVec(end);
            % %-%-xsTop = xsoftheta(thetaplot,A2New,A3New);
            % %-%-zsTop = zsoftheta(thetaplot,A2New,A3New);
            % zsTop = zs_from_spherical(thetaplot, amplitudes_new);
            % xsTop = r_from_spherical(thetaplot, amplitudes_new); 
            % plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop]+zTent,'k','Linewidth',2);
            % width = min(nr, 200);
            % plot([-fliplr(xplot(2:width)),xplot(1:width)],[flipud(eta_accepted(2:width));eta_accepted(1:width)],'LineWidth',2);
            % hold off
            axis equal
            title(sprintf('   t = %0.3f s, z = %.2f', recordedConditions{current_index}.time*T_unit, z*L_unit),'FontSize',16);
            grid on
            hold off;
            %set(gca,'xlim',[-6 6])
            drawnow;
        end
    % 
  
     end % Outer while
    % 
    % 
    if ~isnan(PROBLEM_CONSTANTS.precomputedInverse)
        cd(fold)
        precomputedInverse = PROBLEM_CONSTANTS.precomputedInverse;
        save(myfile, "precomputedInverse")
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, "precomputedInverse");
    end

    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end

    results_saver("simulationResults", 1:(current_index-1), variableValues, savingvarNames, NameValueArgs);

catch ME

    if ~isnan(PROBLEM_CONSTANTS.precomputedInverse)
        cd(fold)
        precomputedInverse = PROBLEM_CONSTANTS.precomputedInverse;
        save(myfile, "precomputedInverse")
        PROBLEM_CONSTANTS = rmfield(PROBLEM_CONSTANTS, "precomputedInverse");
    end

    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end

    results_saver("errored_results", 1:(current_index-1), variableValues, savingvarNames, NameValueArgs);
       
    fprintf("Couldn't run simulation"); 
    
    save(sprintf("error_log%s.mat", datetimeMarker),'ME');
end % end while catch

%simul_time = toc(tstart);
%simul_time = simul_time - tstart;



%results_saver("", 1:(current_index-1), variableValues, savingvarNames);
%save('ProblemConditions.mat', "NameValueArgs", ...
%    "dt", "UNITS", "PROBLEM_CONSTANTS");
mypwd = split(pwd, "1_code"); mypwd = mypwd{2};
fprintf("Finished simulation on %s. Time elapsed: %0.2f minutes\n", mypwd, simul_time/60);
cd(currfold)


end

function results_saver(fileName, indexes, variables, variableNames, NameValueArgs)
    currfold = pwd;
    folders = { ...
        sprintf("rho%.2fgcm3-sigma%.2fdynecm-nu%%.4fSt", ...
        NameValueArgs.bathDensity, NameValueArgs.bathSurfaceTension, NameValueArgs.bathViscosity), ...
        sprintf("diskRadius%.2gcm-diskMass%.2gg", NameValueArgs.diskRadius, NameValueArgs.diskMass), ...
        sprintf("forceAmplitude%.2gdyne-forceFrequency%gHz", NameValueArgs.forceAmplitude, NameValueArgs.forceFrequency)
    };
    for ii = 1:length(folders)
        folder = folders{ii};
        if ~exist(folder, "dir"); mkdir(folder);  end
        cd(folder);
    end

    if indexes(2) == 1
        indexes = indexes(2:end); 
    end
    for ii = 1:length(variables)
       var = variables{ii};
       switch length(size(var))
           case 1
               if ~isstruct(var)
                   var = var(indexes);
               end
           case 2
               if iscell(var)
                   var = var{:, indexes};
               else
                   var = var(:, indexes);
               end
       end
       stru = struct(variableNames{ii}, var);
       save(sprintf('%s%s.mat', fileName, NameValueArgs.datetimeMarker), '-struct', 'stru', '-append');
    end
    cd(currfold)
end

function out = getVarName(var)
    out = inputname(1);
end

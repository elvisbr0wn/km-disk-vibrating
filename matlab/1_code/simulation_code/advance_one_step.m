function next_condition = advance_one_step(previous_conditions, PROBLEM_CONSTANTS)

    % If the matrix exists, import it

    % If it does not, create it and then save it
    dt = previous_conditions.dt;
    We = PROBLEM_CONSTANTS.weber;
    Re = PROBLEM_CONSTANTS.reynolds;
    Fr = PROBLEM_CONSTANTS.froude;
    Delta = PROBLEM_CONSTANTS.laplacian;
    DTN = PROBLEM_CONSTANTS.DTN;
    nr = PROBLEM_CONSTANTS.nr;
    cPoints = PROBLEM_CONSTANTS.contact_points;
    pIntegral = PROBLEM_CONSTANTS.pressure_integral;
    Ma = PROBLEM_CONSTANTS.obj_mass;
    F = PROBLEM_CONSTANTS.force_amplitude;
    t = previous_conditions.time + dt;
    w = PROBLEM_CONSTANTS.force_frequency;
    CoM = previous_conditions.center_of_mass;
    CoM_vel = previous_conditions.center_of_mass_velocity;
    % Preparing the matrix (2x2 block)
    Sist = [[eye(nr)-dt*2*Delta/Re,-dt*DTN];...
        [dt*(eye(nr)/Fr-Delta/We),eye(nr)-dt*2*Delta/Re]];
      
    % Completing the system
    Mat =  [[Sist(:,(cPoints+1):2*nr),...
        [zeros(nr,cPoints);dt*eye(cPoints);zeros(nr-cPoints,cPoints)],...
        zeros(2*nr,1),Sist(:,1:cPoints)*ones(cPoints,1)];
        [zeros(1,2*nr-cPoints),-dt*pIntegral(1:cPoints)/Ma, 1 ,0];
        [zeros(1,2*nr-cPoints),-zeros(1, cPoints)  ,-dt,1]];
    
    b = [previous_conditions.bath_surface; previous_conditions.bath_potential];
    %b = b-Sist(:,1:cPoints)*(CoM * ones(cPoints, 1)); %zs(1:cPoints)+Rv);
    indep = [b; CoM_vel - dt/Fr - dt*F*sin(w*(t+dt)); CoM];
    
    sol = Mat\indep;

    next_condition = previous_conditions;
    next_condition.bath_surface = [sol(end)* ones(cPoints, 1); ...
                                                sol(1:nr-cPoints)];
    next_condition.bath_potential = sol(nr-cPoints+1:2*nr-cPoints);
    next_condition.pressure = sol(2*nr-cPoints+1:2*nr);
    next_condition.center_of_mass_velocity = sol(end-1);
    next_condition.center_of_mass = sol(end);
    next_condition.time = next_condition.time + dt;
    
end
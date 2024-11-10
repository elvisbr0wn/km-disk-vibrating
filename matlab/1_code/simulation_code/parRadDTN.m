function DTNnew345 = parRadDTN(nr, D)
tic
maxtime = 10*60*60;
poolnum = 4; %Size of the parallel pool
workerload = 10;%Number of tasks that we expect worker to run before saving
%load('runNumber.mat','runNumber')
if true || runNumber == 0
    
    %load('nr.mat','nr')
    %load('D.mat','D')
    %load('rn.mat','rn')
    rn = 0:nr+1;
    dr = D/(2*nr);
    %load('dr.mat','dr')
    
    refp = 10;
    %save('refp.mat','refp')

    %Finding the DTN operator
    drp = dr/refp;
    numer = ceil(pi*D/drp);
    if mod(numer,2)==1
        numer = numer+1;
    end
    dtheta = 2*pi/numer;%use pi/even number
    DTNnew345=zeros(nr,nr);

    
    %% Integrating away from the singularity
    k = 1;
    for ii=2:(rn(k)+nr+1)
        idx1 = round(ii);
        if idx1<nr+1
            DTNnew345(k,idx1)   = DTNnew345(k,idx1)           - (( - ii^2/2 -  ii -1/3)*log((ii+1)/(ii)) +( ii^2/6 +ii/2 +1/3)     *(1-ii/(ii+1)) -(ii+.5)/6 +ii/2 +1/2);
        end
        if idx1<nr
            DTNnew345(k,idx1+1) = DTNnew345(k,idx1+1)         - (( 3*ii^2/2 +2*ii -1/2)*log((ii+1)/(ii)) +(-ii^2/2 -ii   +1/2 +1/ii)*(1-ii/(ii+1)) +(ii+.5)/2 -3*ii/2   -1);
            if idx1<nr-1
                DTNnew345(k,idx1+2) = DTNnew345(k,idx1+2)     - ((-3*ii^2/2 -  ii +1  )*log((ii+1)/(ii)) +( ii^2/2 +ii/2 - 1)      *(1-ii/(ii+1)) -(ii+.5)/2 +3*ii/2   +1/2);
                if idx1<nr-2
                    DTNnew345(k,idx1+3) = DTNnew345(k,idx1+3) - ((   ii^2/2      -1/6)*log((ii+1)/(ii)) +(-ii^2        +1)/6     *(1-ii/(ii+1)) +(ii+.5)/6 -ii/2);
                end
            end
        end
    end
    DTNnew345(1,1) = DTNnew345(1,1) + 1/2;
    DTNnew345(1,:) = DTNnew345(1,:)/dr;
    %Integrating the vincinity of the origin (= sing)
    DTNnew345(1,1)   = DTNnew345(k,1) +209/(54*dr);%- (-15+    4)/(6*dr);
    DTNnew345(1,2)   = DTNnew345(k,2) -29/(6*dr);%+ (-16+4/3*4)/(6*dr);
    DTNnew345(1,3)   = DTNnew345(k,3) + 7/(6*dr);%+ (  1-1/3*4)/(6*dr);
    DTNnew345(1,4)   = DTNnew345(k,4) - 11/(54*dr);%- (-15+    4)/(6*dr);

    
    %% Integrating the vincinity of the origin
    k = 2;
    for ii = 1:2*refp
        Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
            x1 = ii*cos(l)/refp;
            posr = radn-rn(k);
            DTNnew345(k,k)   = DTNnew345(k,k) - ((  -6- 3*posr+9*posr^2+ 3*posr^3- 3*posr^4)*posr+6*x1)/72   *Kern;
            DTNnew345(k,k-1) = DTNnew345(k,k-1) - ((-88+48*posr+62*posr^2-12*posr^3-10*posr^4)*posr+88*x1)/72   *Kern;
            DTNnew345(k,k)   = DTNnew345(k,k)   - ((  72-90*posr-90*posr^2+18*posr^3+18*posr^4)*posr- 72*x1)/72    *Kern;
            if k<nr-.5
                DTNnew345(k,k+1) = DTNnew345(k,k+1) - (( 24+48*posr+18*posr^2-12*posr^3-6*posr^4)*posr-24*x1)/72   *Kern;
                if k < nr-1.5
                    DTNnew345(k,k+2) = DTNnew345(k,k+2) - ((  -2- 3*posr+ posr^2+ 3*posr^3+ posr^4)*posr+ 2*x1)/72   *Kern;
                end
            end    
        end
    end
    %% Integrating away from the singularity
    for ii=2*refp+1:(rn(k)+nr)*refp
        Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
            idx1 = floor(radn);
            w1 = min(max(0,radn - idx1),1);
            if idx1 < .5
                    DTNnew345(k,idx1+1) = DTNnew345(k,idx1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                    DTNnew345(k,idx1+2) = DTNnew345(k,idx1+2) - (-w1+2)*w1^2*             Kern;
                    DTNnew345(k,idx1+3) = DTNnew345(k,idx1+3) - (w1-1)*w1^2/4*            Kern;
            elseif idx1<nr
                DTNnew345(k,idx1)   = DTNnew345(k,idx1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                DTNnew345(k,idx1+1) = DTNnew345(k,idx1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                if idx1<nr-1
                    DTNnew345(k,idx1+2) = DTNnew345(k,idx1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                    if idx1<nr-2
                        DTNnew345(k,idx1+3) = DTNnew345(k,idx1+3) - (w1^2-1)*w1/6*    Kern;
                    end
                end
            end
        end
    end
    DTNnew345(k,:) = dtheta/(2*pi*drp)*DTNnew345(k,:);
    DTNnew345(k,k) = DTNnew345(k,k) + 2/(4*dr+drp);   


    %% Integrating the vincinity of the origin
    k = 3;
    for ii = 1:2*refp
        Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
            x1 = ii*cos(l)/refp;
            posr = radn-rn(k);
            DTNnew345(k,k-2) = DTNnew345(k,k-2) - (( 124- 12*posr-149*posr^2+ 12*posr^3+ 25*posr^4)*posr-124*x1)/288   *Kern;
            DTNnew345(k,k-1) = DTNnew345(k,k-1) - ((-384+192*posr+288*posr^2- 48*posr^3- 48*posr^4)*posr+384*x1)/288   *Kern;
            DTNnew345(k,k)   = DTNnew345(k,k)   + ((-  4+ 10*posr+  5*posr^2-  2*posr^3-    posr^4)*posr+  4*x1)/8     *Kern;
            if k<nr-.5
                DTNnew345(k,k+1) = DTNnew345(k,k+1) - (( 128+192*posr+ 32*posr^2- 48*posr^3- 16*posr^4)*posr-128*x1)/288   *Kern;
                if k < nr-1.5
                    DTNnew345(k,k+2) = DTNnew345(k,k+2) - ((- 12- 12*posr+  9*posr^2+ 12*posr^3+ 3*posr^4)*posr+ 12*x1)/288   *Kern;
                end
            end    
        end
    end
    %% Integrating away from the singularity
    for ii=2*refp+1:(rn(k)+nr)*refp
        Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
        for l=dtheta/2:dtheta:pi-dtheta/4
            radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
            idx1 = floor(radn);
            w1 = min(max(0,radn - idx1),1);
            if idx1 < .5
                    DTNnew345(k,idx1+1) = DTNnew345(k,idx1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                    DTNnew345(k,idx1+2) = DTNnew345(k,idx1+2) - (-w1+2)*w1^2*             Kern;
                    DTNnew345(k,idx1+3) = DTNnew345(k,idx1+3) - (w1-1)*w1^2/4*            Kern;
            elseif idx1<nr
                DTNnew345(k,idx1)   = DTNnew345(k,idx1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                DTNnew345(k,idx1+1) = DTNnew345(k,idx1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                if idx1<nr-1
                    DTNnew345(k,idx1+2) = DTNnew345(k,idx1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                    if idx1<nr-2
                        DTNnew345(k,idx1+3) = DTNnew345(k,idx1+3) - (w1^2-1)*w1/6*    Kern;
                    end
                end
            end
        end
    end
    DTNnew345(k,:) = dtheta/(2*pi*drp)*DTNnew345(k,:);
    DTNnew345(k,k) = DTNnew345(k,k) + 2/(4*dr+drp); 
    
    % Shut down existing parallel pool.
    % if ~isempty(gcp('nocreate'))
    %     delete(gcp);
    % else
    %     disp('No pool exists.')
    % end
    % % Create a new parallel pool.
    % if isempty(gcp('nocreate'))
    %     parpool(poolnum)
    % else
    %     disp('A pool already exists.')
    % end
    for ll = 4:poolnum*workerload:nr 
        for k = ll:min(ll+poolnum*workerload-1,nr)
            k/nr;
            Line = zeros(1,nr);
            for ii = 1:2*refp
                Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
                for l=dtheta/2:dtheta:pi-dtheta/4
                    radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
                    x1 = ii*cos(l)/refp; 
                    posr = radn-rn(k);
                    Line(k-2) = Line(k-2) - (  ( 2-  posr-2*posr^2+  posr^3)*posr- 2*x1)/24   *Kern;
                    Line(k-1) = Line(k-1) - (4*(-4+4*posr+  posr^2-  posr^3)*posr+16*x1)/24   *Kern;
                    Line(k)   = Line(k)   - (posr^2-5)*posr^2/4                               *Kern;
                    if k<nr-.5
                        Line(k+1) = Line(k+1) -  (4*(4+4*posr- posr^2- posr^3)*posr-16*x1)/24 *Kern;
                        if k < nr-1.5
                            Line(k+2) = Line(k+2) - ((-2-posr+2*posr^2+posr^3)*posr+ 2*x1)/24 *Kern;
                        end
                    end    
                end
            end 
            
            for ii=2*refp+1:(rn(k)+nr)*refp
                Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
                for l=dtheta/2:dtheta:pi-dtheta/4
                    radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
                    idx1 = floor(radn); 
                    w1 = min(max(0,radn - idx1),1);
                    if idx1 < .5
                        Line(idx1+1) = Line(idx1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                        Line(idx1+2) = Line(idx1+2) - (-w1+2)*w1^2*             Kern;
                        Line(idx1+3) = Line(idx1+3) - (w1-1)*w1^2/4*            Kern;
                    elseif idx1<nr
                        Line(idx1)   = Line(idx1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                        Line(idx1+1) = Line(idx1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                        if idx1<nr-1
                            Line(idx1+2) = Line(idx1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                            if idx1<nr-2
                                Line(idx1+3) = Line(idx1+3) - (w1^2-1)*w1/6*    Kern;
                            end
                        end
                    end
                end
            end
            
            Line(:) = dtheta/(2*pi*drp)*Line(:);
            Line(k) = Line(k) + 2/(4*dr+drp);

            DTNnew345(k,:) = Line;
        end
        
        
        % save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
        % 
        % runtime = toc;
        % if runtime >= maxtime && false
        %     reStartAt = ll + poolnum*workerload;
        %     save('reStartAt.mat','reStartAt')
        %     runNumber = 1;
        %     save('runNumber.mat','runNumber')
        %     % Shut down the pool.
        %     if ~isempty(gcp('nocreate'))
        %         delete(gcp);
        %     else
        %         disp('No pool exists.')
        %     end
        %     return
        % end 
    end
    
    % save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
    % reStartAt = nr+1;
    % save('reStartAt.mat','reStartAt')
    % runNumber = 0;
    % save('runNumber.mat','runNumber')
    % % Shut down the pool.
    % if ~isempty(gcp('nocreate'))
    %     delete(gcp);
    % else
    %     disp('No pool exists.')
    % end  
    return
elseif runNumber > 0
    error("I should not be here!")
    load('nr.mat','nr')
    load('D.mat','D')
    load('rn.mat','rn')
    load('dr.mat','dr')
    load('refp.mat','refp')
    load('reStartAt.mat','reStartAt')
    
    drp = dr/refp;
    numer = ceil(pi*D/drp);
    if mod(numer,2)==1;
        numer = numer+1;
    end
    dtheta = 2*pi/numer;%use pi/even number
    
    if reStartAt > nr
        return
    else
        load(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
        for ll = reStartAt:poolnum*workerload:nr 
            parfor k = ll:min(ll+poolnum*workerload-1,nr)
                k/nr
                Line = zeros(1,nr);
                for ii = 1:2*refp
                    Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
                    for l=dtheta/2:dtheta:pi-dtheta/4
                        radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
                        x1 = ii*cos(l)/refp;
                        posr = radn-rn(k);
                        Line(k-2) = Line(k-2) - (  ( 2-  posr-2*posr^2+  posr^3)*posr- 2*x1)/24   *Kern;
                        Line(k-1) = Line(k-1) - (4*(-4+4*posr+  posr^2-  posr^3)*posr+16*x1)/24   *Kern;
                        Line(k)   = Line(k)   - (posr^2-5)*posr^2/4                               *Kern;
                        if k<nr-.5
                            Line(k+1) = Line(k+1) -  (4*(4+4*posr- posr^2- posr^3)*posr-16*x1)/24 *Kern;
                            if k < nr-1.5
                                Line(k+2) = Line(k+2) - ((-2-posr+2*posr^2+posr^3)*posr+ 2*x1)/24 *Kern;
                            end
                        end    
                    end
                end     
                for ii=2*refp+1:(rn(k)+nr)*refp
                    Kern = 2*(1/(ii-1/2)-1/(ii+1/2));
                    for l=dtheta/2:dtheta:pi-dtheta/4
                        radn = abs(sqrt((rn(k)+ii*cos(l)/refp)^2+(ii*sin(l)/refp)^2));
                        idx1 = floor(radn);
                        w1 = min(max(0,radn - idx1),1);
                        if idx1 < .5
                            Line(idx1+1) = Line(idx1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                            Line(idx1+2) = Line(idx1+2) - (-w1+2)*w1^2*             Kern;
                            Line(idx1+3) = Line(idx1+3) - (w1-1)*w1^2/4*            Kern;
                        elseif idx1<nr
                            Line(idx1)   = Line(idx1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
                            Line(idx1+1) = Line(idx1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
                            if idx1<nr-1
                                Line(idx1+2) = Line(idx1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                                if idx1<nr-2
                                    Line(idx1+3) = Line(idx1+3) - (w1^2-1)*w1/6*    Kern;
                                end
                            end
                        end
                    end
                end
                Line(:) = dtheta/(2*pi*drp)*Line(:);
                Line(k) = Line(k) + 2/(4*dr+drp);
                DTNnew345(k,:) = Line;
            end
            save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
            runtime = toc;
            if runtime >= maxtime
                reStartAt = ll + poolnum*workerload;
                save('reStartAt.mat','reStartAt')
                runNumber = runNumber+1;
                save('runNumber.mat','runNumber')
                % Shut down the pool.
                if ~isempty(gcp('nocreate'))
                    delete(gcp);
                else
                    disp('No pool exists.')
                end
                return
            end
        end
        save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
        reStartAt = nr+1;
        save('reStartAt.mat','reStartAt')
        runNumber = runNumber+1;
        save('runNumber.mat','runNumber')
        % Shut down the pool.
        if ~isempty(gcp('nocreate'))
            delete(gcp);
        else
            disp('No pool exists.')
        end  
        return
    end
end
end % End function
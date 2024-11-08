function DTNnew345 = DTNVectorized(nr, D)    

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
% Define the ii range as a vector and compute idx1 for all values of ii at once
ii_vals = 2:(rn(k) + nr + 1);
idx1_vals = round(ii_vals);

% Preallocate arrays for computed values for each term
term_0 = (-ii_vals.^2 / 2 - ii_vals - 1/3) .* log((ii_vals + 1) ./ ii_vals) ...
         + (ii_vals.^2 / 6 + ii_vals / 2 + 1/3) .* (1 - ii_vals ./ (ii_vals + 1)) ...
         - (ii_vals + 0.5) / 6 + ii_vals / 2 + 1/2;

term_1 = (3 * ii_vals.^2 / 2 + 2 * ii_vals - 1/2) .* log((ii_vals + 1) ./ ii_vals) ...
         + (-ii_vals.^2 / 2 - ii_vals + 1/2 + 1 ./ ii_vals) .* (1 - ii_vals ./ (ii_vals + 1)) ...
         + (ii_vals + 0.5) / 2 - 3 * ii_vals / 2 - 1;

term_2 = (-3 * ii_vals.^2 / 2 - ii_vals + 1) .* log((ii_vals + 1) ./ ii_vals) ...
         + (ii_vals.^2 / 2 + ii_vals / 2 - 1) .* (1 - ii_vals ./ (ii_vals + 1)) ...
         - (ii_vals + 0.5) / 2 + 3 * ii_vals / 2 + 1/2;

term_3 = (ii_vals.^2 / 2 - 1/6) .* log((ii_vals + 1) ./ ii_vals) ...
         + (-ii_vals.^2 + 1) / 6 .* (1 - ii_vals ./ (ii_vals + 1)) ...
         + (ii_vals + 0.5) / 6 - ii_vals / 2;

% Initialize a temporary matrix to accumulate updates
DTN_update = zeros(1, nr);

% Collect indices and values for accumarray
indices_0 = idx1_vals(idx1_vals < nr + 1);
values_0 = term_0(idx1_vals < nr + 1);

indices_1 = idx1_vals(idx1_vals < nr) + 1;
values_1 = term_1(idx1_vals < nr);

indices_2 = idx1_vals(idx1_vals < nr - 1) + 2;
values_2 = term_2(idx1_vals < nr - 1);

indices_3 = idx1_vals(idx1_vals < nr - 2) + 3;
values_3 = term_3(idx1_vals < nr - 2);

% Use accumarray to aggregate values at each index in DTN_update
DTN_update = DTN_update - accumarray(indices_0.', values_0.', [nr, 1]).' ...
                         - accumarray(indices_1.', values_1.', [nr, 1]).' ...
                         - accumarray(indices_2.', values_2.', [nr, 1]).' ...
                         - accumarray(indices_3.', values_3.', [nr, 1]).';

% Apply the update to DTNnew345 for row k
DTNnew345(k, :) = DTNnew345(k, :) + DTN_update;


DTNnew345(1,1) = DTNnew345(1,1) + 1/2;
DTNnew345(1,:) = DTNnew345(1,:)/dr;
%Integrating the vincinity of the origin (= sing)
DTNnew345(1,1)   = DTNnew345(k,1) +209/(54*dr);%- (-15+    4)/(6*dr);
DTNnew345(1,2)   = DTNnew345(k,2) -29/(6*dr);%+ (-16+4/3*4)/(6*dr);
DTNnew345(1,3)   = DTNnew345(k,3) + 7/(6*dr);%+ (  1-1/3*4)/(6*dr);
DTNnew345(1,4)   = DTNnew345(k,4) - 11/(54*dr);%- (-15+    4)/(6*dr);


%% Integrating the vincinity of the origin
k = 2;
% Define vector of i values and compute Kern for all values
i_vals = 1:2*refp;
Kern = 2 * (1 ./ (i_vals - 1/2) - 1 ./ (i_vals + 1/2));

% Define vector for l over the specified range and compute trigonometric terms
l_vals = dtheta/2 : dtheta : pi - dtheta/4;
cos_l = cos(l_vals);
sin_l = sin(l_vals);

% Compute radn, x1, and posr vectors
radn = abs(sqrt((rn(k) + i_vals' * cos_l / refp).^2 + (i_vals' * sin_l / refp).^2));
x1 = i_vals' * cos_l / refp;
posr = radn - rn(k);

% Vectorized updates for DTNnew345

% Update DTNnew345(k, k)
DTNnew345(k, k) = DTNnew345(k, k) - sum((( -6 - 3 * posr + 9 * posr.^2 + 3 * posr.^3 - 3 * posr.^4) .* posr + 6 * x1) / 72 .* Kern', 'all');
DTNnew345(k, k-1) = DTNnew345(k, k-1) - sum(((-88 + 48 * posr + 62 * posr.^2 - 12 * posr.^3 - 10 * posr.^4) .* posr + 88 * x1) / 72 .* Kern', 'all');
DTNnew345(k, k) = DTNnew345(k, k) - sum((( 72 - 90 * posr - 90 * posr.^2 + 18 * posr.^3 + 18 * posr.^4) .* posr - 72 * x1) / 72 .* Kern', 'all');

% Condition k < nr - 0.5
if k < nr - 0.5
    DTNnew345(k, k+1) = DTNnew345(k, k+1) - sum((( 24 + 48 * posr + 18 * posr.^2 - 12 * posr.^3 - 6 * posr.^4) .* posr - 24 * x1) / 72 .* Kern', 'all');
end

% Condition k < nr - 1.5
if k < nr - 1.5
    DTNnew345(k, k+2) = DTNnew345(k, k+2) - sum(((-2 - 3 * posr + posr.^2 + 3 * posr.^3 + posr.^4) .* posr + 2 * x1) / 72 .* Kern', 'all');
end
%% Integrating away from the singularity
for i=2*refp+1:(rn(k)+nr)*refp
    Kern = 2*(1/(i-1/2)-1/(i+1/2));
    for l=dtheta/2:dtheta:pi-dtheta/4
        radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
        i1 = floor(radn);
        w1 = min(max(0,radn - i1),1);
        if i1 < .5
                DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1+2)*w1^2*             Kern;
                DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1-1)*w1^2/4*            Kern;
        elseif i1<nr
            DTNnew345(k,i1)   = DTNnew345(k,i1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
            DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
            if i1<nr-1
                DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                if i1<nr-2
                    DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1^2-1)*w1/6*    Kern;
                end
            end
        end
    end
end
DTNnew345(k,:) = dtheta/(2*pi*drp)*DTNnew345(k,:);
DTNnew345(k,k) = DTNnew345(k,k) + 2/(4*dr+drp);   


%% Integrating the vincinity of the origin
k = 3;
for i = 1:2*refp
    Kern = 2*(1/(i-1/2)-1/(i+1/2));
    for l=dtheta/2:dtheta:pi-dtheta/4
        radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
        x1 = i*cos(l)/refp;
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
for i=2*refp+1:(rn(k)+nr)*refp
    Kern = 2*(1/(i-1/2)-1/(i+1/2));
    for l=dtheta/2:dtheta:pi-dtheta/4
        radn = abs(sqrt((rn(k)+i*cos(l)/refp)^2+(i*sin(l)/refp)^2));
        i1 = floor(radn);
        w1 = min(max(0,radn - i1),1);
        if i1 < .5
                DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (3*w1^3/4-7*w1^2/4+1)*    Kern;
                DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1+2)*w1^2*             Kern;
                DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1-1)*w1^2/4*            Kern;
        elseif i1<nr
            DTNnew345(k,i1)   = DTNnew345(k,i1)   - (-w1^2/6+w1/2-1/3)*w1*    Kern;
            DTNnew345(k,i1+1) = DTNnew345(k,i1+1) - (w1^3/2-w1^2-w1/2+1)*     Kern;
            if i1<nr-1
                DTNnew345(k,i1+2) = DTNnew345(k,i1+2) - (-w1^2/2+w1/2+1)*w1*  Kern;
                if i1<nr-2
                    DTNnew345(k,i1+3) = DTNnew345(k,i1+3) - (w1^2-1)*w1/6*    Kern;
                end
            end
        end
    end
end
DTNnew345(k,:) = dtheta/(2*pi*drp)*DTNnew345(k,:);
DTNnew345(k,k) = DTNnew345(k,k) + 2/(4*dr+drp); 

% % Shut down existing parallel pool.
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
checkpoints = (1:19)/20;
for pp = 1:length(checkpoints)
[~, checkpoints(pp)] = min(abs((4:nr) - nr*(checkpoints(pp) + .01))); 
end
parfor k = 4:nr
    if ismember(k, checkpoints); disp(k/nr); end
    Line = zeros(1,nr);
    
    % Define vectors for ii and l values
    ii_vals = (1:2*refp).';                       % Column vector for ii (transposed to column)
    l_vals = dtheta/2 : dtheta : pi - dtheta/4;   % Row vector for l
    
    % Calculate Kern for all ii values
    Kern_vals = 2 * (1 ./ (ii_vals - 1/2) - 1 ./ (ii_vals + 1/2));  % Column vector
    
    % Compute radn, x1, and posr for each (ii, l) combination
    cos_l = cos(l_vals);  % Row vector of cosines of l
    sin_l = sin(l_vals);  % Row vector of sines of l
    
    radn_vals = abs(sqrt((rn(k) + ii_vals * cos_l / refp).^2 + (ii_vals * sin_l / refp).^2));
    x1_vals = ii_vals * cos_l / refp;  % Matrix of x1 values
    posr_vals = radn_vals - rn(k);     % Matrix of posr values
    
    % Vectorized update for Line(k-2)
    Line(k-2) = Line(k-2) - sum(((2 - posr_vals - 2 * posr_vals.^2 + posr_vals.^3) .* posr_vals - 2 * x1_vals) ./ 24 .* Kern_vals, 'all');
    
    % Vectorized update for Line(k-1)
    Line(k-1) = Line(k-1) - sum((4 * (-4 + 4 * posr_vals + posr_vals.^2 - posr_vals.^3) .* posr_vals + 16 * x1_vals) ./ 24 .* Kern_vals, 'all');
    
    % Vectorized update for Line(k)
    Line(k) = Line(k) - sum((posr_vals.^2 - 5) .* posr_vals.^2 / 4 .* Kern_vals, 'all');
    
    % Conditional vectorized updates for Line(k+1) and Line(k+2)
    if k < nr - 0.5
        Line(k+1) = Line(k+1) - sum((4 * (4 + 4 * posr_vals - posr_vals.^2 - posr_vals.^3) .* posr_vals - 16 * x1_vals) ./ 24 .* Kern_vals, 'all');
        if k < nr - 1.5
            Line(k+2) = Line(k+2) - sum(((-2 - posr_vals + 2 * posr_vals.^2 + posr_vals.^3) .* posr_vals + 2 * x1_vals) ./ 24 .* Kern_vals, 'all');
        end
    end


    % Define ranges for ii and l
    ii = ((2*refp+1):((rn(k)+nr)*refp))';
    l_vals = (dtheta/2):dtheta:(pi-dtheta/4);
    
    % Compute Kernels and radn values for all combinations of ii and l
    Kern = (2 * (1 ./ (ii - 1/2) - 1 ./ (ii + 1/2)));
    %[ii_grid, l_grid] = ndgrid(ii, l_vals);
    radn = abs(sqrt((rn(k) + ii .* cos(l_vals) / refp).^2 + (ii .* sin(l_vals) / refp).^2));
    idxs = floor(radn);
    w1 = min(max(0, radn - idxs), 1);
    
    % Initialize Line vector update values
    LineUpdate = zeros(1, nr);
    sliceAndSum = @(exp, cond, idx, nr) accumarray(idx(cond), exp(cond), [nr 1])';
    
    % Accumulate updates for each index in Line based on conditions
    cond1 = idxs < 0.5; 
    LineUpdate = LineUpdate + sliceAndSum(-(3 * w1.^3 / 4 - 7 * w1.^2 / 4 + 1) .* Kern, cond1, idxs+1, nr);
    LineUpdate = LineUpdate + sliceAndSum(-(-w1 + 2) .* w1.^2 .* Kern, cond1, idxs+2, nr);
    LineUpdate = LineUpdate + sliceAndSum(-(w1 - 1) .* w1.^2 / 4 .* Kern, cond1, idxs+3, nr);
    
    
    cond2 = (idxs >= 0.5) & (idxs < nr); 
    
    LineUpdate = LineUpdate + sliceAndSum(-(-w1.^2 / 6 + w1 / 2 - 1/3) .* w1 .* Kern, cond2, idxs, nr);
    LineUpdate = LineUpdate + sliceAndSum(-(w1.^3 / 2 - w1.^2 - w1 / 2 + 1) .* Kern, cond2, idxs+1, nr);
    
    cond3 = (idxs < nr - 1) & (idxs >= 0.5); 
    cond4 = (idxs < nr - 2) & (idxs >= 0.5); 
    LineUpdate = LineUpdate + sliceAndSum(-(-w1.^2 / 2 + w1 / 2 + 1) .* w1 .* Kern, cond3, idxs+2, nr);
    LineUpdate = LineUpdate + sliceAndSum(-(w1.^2 - 1) .* w1 / 6 .* Kern, cond4, idxs+3, nr);
    
    
    Line = Line + LineUpdate;

    
    % Final adjustments to Line
    Line = dtheta / (2 * pi * drp) * Line;
    Line(k) = Line(k) + 2 / (4 * dr + drp);
    DTNnew345(k,:) = Line;
end
    
   

save(['DTNnew345nr',num2str(nr),'D',num2str(D),'refp',num2str(refp),'.mat'],'DTNnew345')
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

end

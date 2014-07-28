%maxthrough_generate is used to generate random locations and save them.
%For directional case there is a similar file in the other directory
%../maxthrough_directional_ga_multiP5  

clearvars -global;
global M N_f N B C SINRPU_lt SINRSU_lt pwr_pu f pwr_max L_pu L_su L_supu pmin pmax links alpha_su alpha_pu G_M G_S theta;

M =10;   %Number of secondary users
N_f =5; %Number of primary users
B = 6E06;  % Bandwidth
C = 3e08;  % light speed
pwr_max = 10^(-3/10);      % Maximum power budget 
f = 53e06;          % Frequency at the Start of band    
pwr_pu = 10^(6/10);      % Power transmitted by primary user
N = 1e-13;     % Noise level 
SINRPU_lt = 10^(20/10);    % primary user max
SINRSU_lt = 10^(10/10);    % secondary user max 
rng_location = 5000;% Order of location of SU and PU
conflict_count = 0;
%create the input.in file to be used by nsga2.c code which does the real
%multiobjective optimization
createinput(M,N_f);
while(conflict_count<2)
% PU emitter->reciever 1->N_f+1, 2->N_f+2, ......N_f->2_N_f
% PU location 
for(i= 1:N_f)
    pu(i,:) = [rng_location*rand(1,1) rng_location*rand(1,1)];
    flag = 0;
    while(flag ==0)
        pu(N_f+i,:) = [rng_location*rand(1,1) rng_location*rand(1,1)];
        % Pu distance calculation
        d_pu(i,N_f+i) = norm(pu(i,:)-pu(N_f+i,:));
        L_pu(i,N_f+i) = (C*C)/((4*pi*(f+B*(i-1)))^2 *(d_pu(i,N_f+i)^3.5));
        SINR_pu(i) = (pwr_pu*(L_pu(i,N_f+i)/N));
        if(SINR_pu(i) >(SINRPU_lt+1700))  % taking more than 100 to account for interference
            flag =1;
        end
    end
end

% SU location
% SU emitter-> reciever 1-> M+1, 2->M+2, 3->M+3,.......M->2M
for(i = 1:M)
    su(i,:) = [rng_location*rand(1,1),rng_location*rand(1,1)];
    flag = 0;
    while(flag==0)
        su(M+i,:) = [rng_location*rand(1,1),rng_location*rand(1,1)];
        % Distance between transmitter reciever
        d_su(i,M+i) = norm(su(i,:)-su(M+i,:));
        % Reciever is chosen such that effective communication even at
        % highest frequeny. This will take casre of lower frequencies as
        % well. 
        L_su(i,M+i,N_f) = (C*C)/((4*pi*(f+B*(N_f-1)))^2 *(d_su(i,M+i)^4)); 
        % Interference calculation (with all primary users only)
        Imax = -Inf;
        for(j=1:N_f)
            d_supu(M+i,j) = norm(su(M+i,:)-pu(j,:));
           I = N+ pwr_pu*((C*C)/((4*pi*(f+B*(j-1)))^2 *(d_supu(M+i,j)^4)));
            if(I>Imax)
                Imax = I;
            end
        end
        SINR_su_min(i) = (pwr_max*(L_su(i,M+i,N_f)/(N+I)));
        if(SINR_su_min(i)>(SINRSU_lt-5)) 
            flag = 1;
        end
    end
end

% L_su calculation 

for(i = 1:M)
    for(j = 1:M)
        % Distance from su i to su M+j (reciever)
        d_su(i,M+j) = norm(su(i,:)-su(M+j,:));
        for(k = 1:N_f) 
            % Loss due to secondary user i at reciever M+j at frequency k 
            L_su(i,M+j,k) = (C*C)/((4*pi*(f+B*(k-1)))^2*(d_su(i,M+j)^4));
        end
    end
end


% L_supu calculation

for(i = 1:M)
    for(j = 1:N_f)
        d_supu(i,N_f+j) = norm(su(i,:)-pu(N_f+j,:));
        % loss due to secondary user i at primary user N_f+j(only freq.
        % corresponding to primary user N_f+j)
        L_supu(i,N_f+j) = (C*C)/((4*pi*(f+B*(j-1)))^2*(d_supu(i,N_f+j)^4));
    end
end

for(i=1:M)
    for(j=1:N_f)
        freq = f+B*(j-1);
        L_supu(M+i,j) = (C*C)/(((4*pi*freq)^2) *(d_supu(M+i,j)^4));
    end
end

% pass all the L_su, L_pu, L_supu, N,N_F,M,B,C,pwr_max,f, pwr_pu and alpha
% to NSGA2.

input_values = [N N_f M B C pwr_max f pwr_pu SINRPU_lt SINRSU_lt] ;

[fval cons power cons_violation popnum] = main(su,pu,input_values);

'fval = ';
fval(1)
fval(2)
'cons_violations = ';
cons_violation
powerlevels = [];
j = 1;
for(i = power')
    if(i>0)
        powerlevels = [powerlevels j];
    end
    j = j+1;
end
powerlevels
%Save conflicting solution
conflict_count = 0;
if(popnum>1 && cons_violation> -1e-6 && size(powerlevels,2)>1)
    present = zeros([M 1]);
    for(i = 1:popnum)
        present(-fval(2,i)) = 1;
    end
    if(sum(present)>1)
        save(strcat('./conflict_location/su_con',int2str(conflict_count)),'su');
        save(strcat('./conflict_location/pu_con',int2str(conflict_count)),'pu');
        pareto_front = fval(:,1:popnum);
        save(strcat('./conflict_location/fval_con',int2str(conflict_count)),'pareto_front');
        conflict_count = conflict_count+1;
    end    
end

end


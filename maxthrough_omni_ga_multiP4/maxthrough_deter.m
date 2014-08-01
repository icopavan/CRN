%maxthrough_deter is the main function to be called to obtain
%the results for a given location and delta value stored in a mat file. 
%location of the input file is hardcoded and has to be changed at the run time. 
%There is a similar function - maxthrough_generate which is used to generate random
%locations and save them. This function can then be used to regenerate the
%results for those saved locations. 
%This method is for the omni directional case for problem 5.
%Problem 5 is maximizing throughput and maximizing the number of active
%users. For directional case there is a similar file in the other directory
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

%create the input.in file to be used by nsga2.c code which does the real
%multiobjective optimization
createinput(M,N_f);

%For M=2 , N_f = 1 case plotting function 

pu_all = load('../maxthrough_omni_ga_multiP5/conflict_location/pu0.mat','pu');
su_all = load('../maxthrough_omni_ga_multiP5/conflict_location/su0.mat','su');

pu = pu_all.pu;
su = su_all.su;

% PU emitter->reciever 1->N_f+1, 2->N_f+2, ......N_f->2_N_f
% PU location 
for(i= 1:N_f)
	% Pu distance calculation
        d_pu(i,N_f+i) = norm(pu(i,:)-pu(N_f+i,:));
        L_pu(i,N_f+i) = (C*C)/((4*pi*(f+B*(i-1)))^2 *(d_pu(i,N_f+i)^4));
        SINR_pu(i) = (pwr_pu*(L_pu(i,N_f+i)/N));
end

% SU location
% SU emitter-> reciever 1-> M+1, 2->M+2, 3->M+3,.......M->2M
for(i = 1:M)
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
% popnum is the number of rank 1 population members
[fval cons power cons_violation popnum] = main(su,pu,input_values);

'fval = ';
fval(1)
fval(2)
'cons_violations = ';
cons_violation
%Print power of one of the solutions
powerlevels = [];
j = 1;
for(i = power')
    if(i>0)
        powerlevels = [powerlevels j];
    end
    j = j+1;
end
powerlevels

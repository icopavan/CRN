%maxthrough_directional_deter generates the results for given input locations and 
%prints the result

clearvars -global;
global M N_f N B C SINRPU_lt SINRSU_lt pwr_pu f pwr_max L_pu L_su L_supu pmin pmax links alpha_su alpha_pu G_M G_S theta;
M = 10;   %Number of secondary users
N_f =5; %Number of primary users
B = 6E06;  % Bandwidth
C = 3e08;  % light speed
pwr_max = 10^(-3/10);      % Maximum power budget 
f = 53e06;          % Frequency at the Start of band    
pwr_pu = 10^(6/10);      % Power transmitted by primary user
N = 1e-13;     % Noise level 
SINRPU_lt = 10^(20/10);    % primary user max
SINRSU_lt = 10^(10/10);    % secondary user max 
createinput(M,N_f);

%Load locations
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



% G represents anteena gain for transmitter at a reciever.
% G_su(i,j) is gain of su transmitter i when reciever is also another su reciever j.
% G_supu(i,j) is the gain of su transmitter i at pu reciever j. 
% G_supu (i,j) for su reciever and pu transmitter is 1. 
% Defining G for all radios

% for su transmitters
for(i =1:M)
    % transmitter = i
    % reciever = j
    for(j = M+1:2*M)
        % check location to find j
        X = (su(j,1)-su(i,1));
        Y = (su(j,2)-su(i,2));
        % alpha is the angle made by line joining trans reciever with horizontal
        angle = atand(abs(Y)/abs(X));
        if(X>=0 && Y>=0)
            alpha_su(i,j) = angle;
        elseif(X>0 && Y<0)
            alpha_su(i,j) = 360-angle;
        elseif(X<0 && Y>0)
            alpha_su(i,j) = 180-angle;
        elseif(X<0&&Y<0)
            alpha_su(i,j) = 180+angle;
        end
    end
end

%for pu recievers
for(i = 1:M)
    % transmitr_su = i;
    %reciever_pu = j
    for(j = N_f+1:2*N_f)
        X = (pu(j,1)-su(i,1));
        Y = (pu(j,2)-su(i,2));
        angle = atand(abs(Y)/abs(X));
        if(X>=0 && Y>=0)
            alpha_pu(i,j) = angle;
        elseif(X>0 && Y<0)
            alpha_pu(i,j) = 360-angle;
        elseif(X<0 && Y>0)
            alpha_pu(i,j) = 180-angle;
        elseif(X<0&&Y<0)
            alpha_pu(i,j) = 180+angle;
        end
    end
end


theta = 60;
theta_dash = 30;
alph = (theta/360);
eta = (theta_dash/360);
G_M = 1/(alph+eta);
G_S = eta/((1-alph)*(alph+eta));
        
   

% pass all the L_su, L_pu, L_supu, N,N_F,M,B,C,pwr_max,f, pwr_pu and alpha
% to NSGA2. popnum is the number of rank 1 solutions
input_values = [N N_f M B C pwr_max f pwr_pu G_M G_S theta SINRPU_lt SINRSU_lt] ;
[fval cons power cons_violation popnum] = main(su,pu,input_values,alpha_su,alpha_pu);

'fval = ';
fval(1)
fval(2)
'cons_violations = ';
cons_violation
%Printing power levels for one of the rank 1 solution
powerlevels = [];
j = 1;
for(i = power')
    if(i>0)
        powerlevels = [powerlevels j];
    end
    j = j+1;
end
 powerlevels

%define phi
phi = power(M*N_f+1:end);
% plotting the results showing the directions:
% Plotting secondary users 
%h = figure; 
axis manual;
plot(su(:,1),su(:,2),'Marker','.','LineStyle','none');
axis([-1 2 -1 2]);
%set(h,'
for(i = 1:M)
    text(su(i,1),su(i,2),strcat('sut',num2str(i)));
end
for(i = M+1:2*M)
    text(su(i,1),su(i,2),strcat('sur',num2str(i-M)));
end
hold on;
h = plot(pu(:,1),pu(:,2),'r.');
for(i = 1:N_f)
    text(pu(i,1),pu(i,2),strcat('put',num2str(i)));
end

for(i = N_f+1:2*N_f)
    text(pu(i,1),pu(i,2),strcat('pur',num2str(i-N_f)));
end
hold on;
flag = 0;
%r defines the legth of line defining beam
r = 10000;
% Drawing lines showing transmission direction
for(i=1:M)
    flag = 0;
    for(j=1:N_f)
        if(power((i-1)*N_f+j) >0)
            flag = 1;
        end
    end
    if(flag)
        xst   =   [su(i,1) su(i,2)];
        xend1 = xst+[r*cosd(phi(i)-theta/2) r*sind(phi(i)-theta/2)];
        xend2 = xst+[r*cosd(phi(i)+theta/2) r*sind(phi(i)+theta/2)];
        xdata1 = [xst(1);xend1(1)];
        ydata1 = [xst(2);xend1(2)];
        xdata2 = [xst(1);xend2(1)];
        ydata2 = [xst(2);xend2(2)];
        plot(xdata1,ydata1,'-',xdata2,ydata2,'-');
        hold on;
        
    end
end
hold off;

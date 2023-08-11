clear all;
clc

%plug in the known inputs
d_o=0.01;
d_o=0.04;
r_i=0.005;      %m inner seed radius
r_o=0.02;       %m outer flesh radius

k_s=0.3;        %W/m-K
k_f=0.6;        %W/m-K

q_s_gen=50000;  %W/m2 flux from seed
q_f_gen=100000; %W/m2 flux from flesh

R_tc=0.02;      %m2-K/W contact resistance


%% Daytime

T_infd=15+273.15; % K
T_skyd=20+273.15; % K
h=50;             %W/m2-K

G_s=500;          %W/m2
alpha_s=0.75;     %absorptivity of skin
E=0.85;           %emissivity
alpha_emi=0.85;   %radiation aborptivity 


%% NightTime

T_infn=10+273.15; %K
T_skyn=5+273.15;  %K


%% node set up

dr=0.0005;          %difference between nodes
N_seed=(r_i/dr)+1;  %node numbers
N_flesh=(r_o/dr)+2;

delta_T=1;          %max error in temp of current vs previous iteration
tolerance_T=0.0001; %C- tolerance for exiting the loop

%set up coefficient, constant, and temperature matrix
T = ones(N_flesh,1);
c=1;                %iteration counter 


%% While loop for calculation of the derived equations in each node or nodes
while delta_T>tolerance_T   

c=c+1

%Seed radius boundary node-----1st part code
coeff_1 = q_s_gen*4/3*pi*(dr/2)^3;
coeff_2 = 4*pi*(dr/2)^2*k_s/dr;
T(1, c) = (coeff_1 + coeff_2*T(2,c-1))/coeff_2;


%Seed interior nodes-----2nd part code
    for i=2: (N_seed-1)
        r_inner = (i-(3/2))*dr;
        r_outer = (i-(1/2))*dr;
        A1 = (k_s/dr)*4*pi*r_inner^2;
        A2 = (k_s/dr)*4*pi*r_outer^2;
        V = q_s_gen*(4/3)*pi*(r_outer^3 - r_inner^3);
        T(i,c)= (A2*T(i+1,c-1) + A1*T(i-1,c-1) +V)/(A1 + A2);
    end   

    
%Seed side boundary and Flesh side boundary-----3rd part code N_seed
B_contact = (4*pi*r_i^2)/R_tc;
B_q_gen = q_s_gen*(4/3)*pi*(r_i^3 - (r_i - dr/2)^3 );
B_cond = (k_s/dr)*4*pi*r_i^2;
T(N_seed, c) = (B_q_gen + B_cond*T(N_seed -1, c-1) + B_contact*T(N_seed + 1, c-1) )/(B_contact + B_cond );


%Contact flesh side-----3rd part code N_seed+1
B_q_gen = q_f_gen*(4/3)*pi*((r_i + dr/2)^3 - r_i^3 );
B_cond2 = (k_f/dr)*4*pi*(r_i + dr/2)^2;
T(N_seed +1, c) = (B_q_gen + B_cond2*T(N_seed + 2, c-1) + B_contact*T(N_seed, c-1) )/(B_contact + B_cond2 );


%Flesh interior nodes-----4th part code
    for j=(N_seed+2):(N_flesh-1)
        r_inner = (j-(5/2))*dr; %the constant change for contact value, as this value should be larger comparing to previous one
        r_outer = (j-(3/2))*dr; 
        A1 = (k_f/dr)*4*pi*r_inner^2;
        A2 = (k_f/dr)*4*pi*r_outer^2;
        V = q_f_gen*(4/3)*pi*(r_outer^3 - r_inner^3);
        T(j,c)= (A2*T(j+1,c-1) + A1*T(j-1,c-1) +V)/(A1 + A2);
    end 


%Flesh outer boundary node-----5th part code
%F = @(T_M)(E_in + E_g - E_o) for the engergy balance
%the E_g term here is only for the LAST cv of thickness dr/2 
%F = S_cond + S_f_gen +S_radiation - S_conv(daytime) + S_emi(daytime)

M= N_flesh;
F = @(T_M)((k_f/dr)*4*pi*r_o^2*(T(N_flesh -1, c-1) - T_M)   + q_f_gen*(4/3)*pi*(r_o^3 - (r_o - dr/2)^3) + alpha_s*G_s*4*pi*r_o^2 - h*4*pi*r_o^2*(T_M-T_infd) + alpha_emi*(5.67*10^-8)*4*pi*r_o^2*(T_skyd^4 - T_M^4));
intermediate = fzero(F, [-10 1000]);
T(N_flesh, c)=intermediate;
delta_T=max(abs(T(:,c)-T(:, c-1)));

end
Temp_fin = T(:,end)-273.15;


%% plot the Temperature vs radius 

range = linspace(0, r_o, N_flesh);
plot(range, Temp_fin);
title('Plot of temperature vs radius at Day');
xlabel('Radius (m)');
ylabel('Temperature (Degree C)');
grid on;

%Therefore, from the output graph shows the heat is going out as it has a
%temperature drop from inner layers to the outer. 


%% Engergy Balance Calculation Base on Read Numbers for T_s

T_s=307.674284659641; %Surface temperature base on graph (DAY)

q_in= (k_f/dr)*4*pi*r_o^2*(T(N_flesh -1, c-1) - T_s)   + q_f_gen*(4/3)*pi*(r_o^3 - (r_o - dr/2)^3) + alpha_s*G_s*4*pi*r_o^2 
q_out= h*4*pi*r_o^2*(T_s-T_infd) - (alpha_emi*(5.67*10^-8)*4*pi*r_o^2*(T_skyd^4 - T_s^4))




clear;
clc;
close all;

% L_max=10; n_max=10;

Lambda_U=50;
Lambda_V=50;
Lambda_R=2*pi;

%% Main
times=0;
SU_O1_sim=zeros(10,1); % The total success delivery probability under Option 1
SE_O1_sim=zeros(10,1); % The total secrecy delivery probability under Option 1
SU_O2_sim=zeros(10,1); % The total success delivery probability under Option 2
SE_O2_sim=zeros(10,1); % The total secrecy delivery probability under Option 2

for n=1:10
    times=times+1;
    [SU_O1_sim(:,times),SE_O1_sim(:,times)]=Option1_sim(Lambda_V,Lambda_U,Lambda_R,n);
    [SU_O2_sim(:,times),SE_O2_sim(:,times)]=Option2_sim(Lambda_V,Lambda_U,Lambda_R,n);
end

%% Option 1 Simulation
function [SU_O1_sim,SE_O2_sim]=Option1_sim(Lambda_V,Lambda_U,Lambda_R,n)
R_c=1;
beta_c=2^R_c-1;
beta_f=-20;
beta_f=10^(beta_f/10);
R_e=R_c-0.5;
beta_e=2^R_e-1;
Lambda_R=Lambda_R/pi;

L=10;
n_tot=10;
if n<=5
    xi=(n*(3*n-1)+(n_tot-3*n+1)*(2*n-1))/(n_tot-n+1)^2;
else
    xi=1;
end

times_total=300000;
success_times_1=0;
success_times_2=0;
success_times_3=0;
success_times_4=0;
success_times_5=0;
success_times_6=0;
success_times_7=0;
success_times_8=0;
success_times_9=0;
success_times_10=0;

secrecy_times_1=times_total;
secrecy_times_2=times_total;
secrecy_times_3=times_total;
secrecy_times_4=times_total;
secrecy_times_5=times_total;
secrecy_times_6=times_total;
secrecy_times_7=times_total;
secrecy_times_8=times_total;
secrecy_times_9=times_total;
secrecy_times_10=times_total;

time_round_1=times_total;
time_round_2=times_total;
time_round_3=times_total;
time_round_4=times_total;
time_round_5=times_total;
time_round_6=times_total;
time_round_7=times_total;
time_round_8=times_total;
time_round_9=times_total;

for times=1:times_total
    SINR_EVA_2=0;SINR_EVA_3=0;SINR_EVA_4=0;SINR_EVA_5=0;SINR_EVA_6=0;SINR_EVA_7=0;SINR_EVA_8=0;SINR_EVA_9=0;SINR_EVA_10=0;
    [L_number,L_feture]=Road_Generation(Lambda_R);

    SINR_VAM_1=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
    SINR_EVA_1=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
    if SINR_VAM_1>=beta_c
        success_times_1=success_times_1+1;
    else
        SINR_NACK_1=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        if  SINR_NACK_1>=beta_f && L>=2
            time_round_1=time_round_1-1;
            SINR_VAM_2=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
            SINR_EVA_2=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
            if SINR_VAM_2>=beta_c/2
                success_times_2=success_times_2+1;
            else
                SINR_NACK_2=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                if SINR_NACK_2>=beta_f && L>=3
                    time_round_2=time_round_2-1;
                    SINR_VAM_3=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                    SINR_EVA_3=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                    if SINR_VAM_3>=beta_c/3
                        success_times_3=success_times_3+1;
                    else
                        SINR_NACK_3=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                        if SINR_NACK_3>=beta_f && L>=4
                            time_round_3=time_round_3-1;
                            SINR_VAM_4=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                            SINR_EVA_4=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                            if SINR_VAM_4>=beta_c/4
                                success_times_4=success_times_4+1;
                            else
                                SINR_NACK_4=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                if SINR_NACK_4>=beta_f && L>=5
                                    time_round_4=time_round_4-1;
                                    SINR_VAM_5=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                    SINR_EVA_5=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                                    if SINR_VAM_5>=beta_c/5
                                        success_times_5=success_times_5+1;
                                    else
                                        SINR_NACK_5=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                        if SINR_NACK_5>=beta_f && L>=6
                                            time_round_5=time_round_5-1;
                                            SINR_VAM_6=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                            SINR_EVA_6=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                                            if SINR_VAM_6>=beta_c/6
                                                success_times_6=success_times_6+1;
                                            else
                                                SINR_NACK_6=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                if SINR_NACK_6>=beta_f && L>=7
                                                    time_round_6=time_round_6-1;
                                                    SINR_VAM_7=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                    SINR_EVA_7=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                                                    if SINR_VAM_7>=beta_c/7
                                                        success_times_7=success_times_7+1;
                                                    else
                                                        SINR_NACK_7=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                        if SINR_NACK_7>=beta_f && L>=8
                                                            time_round_7=time_round_7-1;
                                                            SINR_VAM_8=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                            SINR_EVA_8=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                                                            if SINR_VAM_8>=beta_c/8
                                                                success_times_8=success_times_8+1;
                                                            else
                                                                SINR_NACK_8=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                                if SINR_NACK_8>=beta_f && L>=9
                                                                    time_round_8=time_round_8-1;
                                                                    SINR_VAM_9=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                                    SINR_EVA_9=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                                                                    if SINR_VAM_9>=beta_c/9
                                                                        success_times_9=success_times_9+1;
                                                                    else
                                                                        SINR_NACK_9=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                                        if SINR_NACK_9>=beta_f && L>=10
                                                                            time_round_9=time_round_9-1;
                                                                            SINR_VAM_10=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
                                                                            SINR_EVA_10=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
                                                                            if SINR_VAM_10>=beta_c/10
                                                                                success_times_10=success_times_10+1;
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if SINR_EVA_1>=beta_e
        secrecy_times_1=secrecy_times_1-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e
        secrecy_times_2=secrecy_times_2-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e
        secrecy_times_3=secrecy_times_3-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e
        secrecy_times_4=secrecy_times_4-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e
        secrecy_times_5=secrecy_times_5-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e
        secrecy_times_6=secrecy_times_6-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e
        secrecy_times_7=secrecy_times_7-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e || SINR_EVA_8>=beta_e
        secrecy_times_8=secrecy_times_8-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e || SINR_EVA_8>=beta_e || SINR_EVA_9>=beta_e
        secrecy_times_9=secrecy_times_9-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e || SINR_EVA_8>=beta_e || SINR_EVA_9>=beta_e || SINR_EVA_10>=beta_e
        secrecy_times_10=secrecy_times_10-1;
    end
    disp(['Option1: Lambda_V=',num2str(Lambda_V),'; Lambda_U=',num2str(Lambda_U),'; Lambda_R=',num2str(Lambda_R),'*pi; n=',num2str(n),'; simulation_times:',num2str(times)]);
end

S_sim_L_1=(success_times_1)/times_total;
S_sim_L_2=(success_times_1+success_times_2)/times_total;
S_sim_L_3=(success_times_1+success_times_2+success_times_3)/times_total;
S_sim_L_4=(success_times_1+success_times_2+success_times_3+success_times_4)/times_total;
S_sim_L_5=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5)/times_total;
S_sim_L_6=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6)/times_total;
S_sim_L_7=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7)/times_total;
S_sim_L_8=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7+success_times_8)/times_total;
S_sim_L_9=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7+success_times_8+success_times_9)/times_total;
S_sim_L_10=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7+success_times_8+success_times_9+success_times_10)/times_total;

SU_O1_sim=[S_sim_L_1;S_sim_L_2;S_sim_L_3;S_sim_L_4;S_sim_L_5;S_sim_L_6;S_sim_L_7;S_sim_L_8;S_sim_L_9;S_sim_L_10];

SE_sim_1=secrecy_times_1/times_total;
SE_sim_2=secrecy_times_2/times_total;
SE_sim_3=secrecy_times_3/times_total;
SE_sim_4=secrecy_times_4/times_total;
SE_sim_5=secrecy_times_5/times_total;
SE_sim_6=secrecy_times_6/times_total;
SE_sim_7=secrecy_times_7/times_total;
SE_sim_8=secrecy_times_8/times_total;
SE_sim_9=secrecy_times_9/times_total;
SE_sim_10=secrecy_times_10/times_total;

SE_O2_sim=[SE_sim_1;SE_sim_2;SE_sim_3;SE_sim_4;SE_sim_5;SE_sim_6;SE_sim_7;SE_sim_8;SE_sim_9;SE_sim_10];
end

%% Option 2 Simulation
function [SU_O2_sim,SE_O2_sim]=Option2_sim(Lambda_V,Lambda_U,Lambda_R,n)
R_c=1;
beta_c=2^R_c-1;
beta_f=-20;
beta_f=10^(beta_f/10);
R_e=R_c-0.5;
beta_e=2^R_e-1;
Lambda_R=Lambda_R/pi;

L=10;
n_tot=10;
if n<=5
    xi=(n*(3*n-1)+(n_tot-3*n+1)*(2*n-1))/(n_tot-n+1)^2;
else
    xi=1;
end

times_total=300000;
success_times_1=0;
success_times_2=0;
success_times_3=0;
success_times_4=0;
success_times_5=0;
success_times_6=0;
success_times_7=0;
success_times_8=0;
success_times_9=0;
success_times_10=0;

secrecy_times_1=times_total;
secrecy_times_2=times_total;
secrecy_times_3=times_total;
secrecy_times_4=times_total;
secrecy_times_5=times_total;
secrecy_times_6=times_total;
secrecy_times_7=times_total;
secrecy_times_8=times_total;
secrecy_times_9=times_total;
secrecy_times_10=times_total;

for times=1:times_total
    SINR_EVA_2=0;SINR_EVA_3=0;SINR_EVA_4=0;SINR_EVA_5=0;SINR_EVA_6=0;SINR_EVA_7=0;SINR_EVA_8=0;SINR_EVA_9=0;SINR_EVA_10=0;
    SINR_ACK_1=0;SINR_ACK_2=0;SINR_ACK_3=0;SINR_ACK_4=0;SINR_ACK_5=0;SINR_ACK_6=0;SINR_ACK_7=0;SINR_ACK_8=0;SINR_ACK_9=0;
    [L_number,L_feture]=Road_Generation(Lambda_R);

    SINR_VAM_1=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
    SINR_EVA_1=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
    if SINR_VAM_1>=beta_c
        success_times_1=success_times_1+1;
        SINR_ACK_1=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
    end
    if SINR_ACK_1<beta_f && L>=2
        SINR_VAM_2=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_2=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_2>=beta_c/2 && SINR_VAM_1<beta_c
            success_times_2=success_times_2+1;
            SINR_ACK_2=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && L>=3
        SINR_VAM_3=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_3=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_3>=beta_c/3 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2
            success_times_3=success_times_3+1;
            SINR_ACK_3=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && L>=4
        SINR_VAM_4=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_4=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_4>=beta_c/4 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2  && SINR_VAM_3<beta_c/3
            success_times_4=success_times_4+1;
            SINR_ACK_4=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && SINR_ACK_4<beta_f && L>=5
        SINR_VAM_5=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_5=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_5>=beta_c/5 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2 && SINR_VAM_3<beta_c/3 && SINR_VAM_4<beta_c/4
            success_times_5=success_times_5+1;
            SINR_ACK_5=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && SINR_ACK_4<beta_f && SINR_ACK_5<beta_f && L>=6
        SINR_VAM_6=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_6=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_6>=beta_c/6 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2 && SINR_VAM_3<beta_c/3 && SINR_VAM_4<beta_c/4 && SINR_VAM_5<beta_c/5
            success_times_6=success_times_6+1;
            SINR_ACK_6=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && SINR_ACK_4<beta_f && SINR_ACK_5<beta_f && SINR_ACK_6<beta_f && L>=7
        SINR_VAM_7=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_7=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_7>=beta_c/7 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2 && SINR_VAM_3<beta_c/3 && SINR_VAM_4<beta_c/4 && SINR_VAM_5<beta_c/5 && SINR_VAM_6<beta_c/6
            success_times_7=success_times_7+1;
            SINR_ACK_7=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && SINR_ACK_4<beta_f && SINR_ACK_5<beta_f && SINR_ACK_6<beta_f && SINR_ACK_7<beta_f && L>=8
        SINR_VAM_8=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_8=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_8>=beta_c/8 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2 && SINR_VAM_3<beta_c/3 && SINR_VAM_4<beta_c/4 && SINR_VAM_5<beta_c/5 && SINR_VAM_6<beta_c/6 && SINR_VAM_7<beta_c/7
            success_times_8=success_times_8+1;
            SINR_ACK_8=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && SINR_ACK_4<beta_f && SINR_ACK_5<beta_f && SINR_ACK_6<beta_f && SINR_ACK_7<beta_f && SINR_ACK_8<beta_f && L>=9
        SINR_VAM_9=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_9=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_9>=beta_c/9 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2 && SINR_VAM_3<beta_c/3 && SINR_VAM_4<beta_c/4 && SINR_VAM_5<beta_c/5 && SINR_VAM_6<beta_c/6 && SINR_VAM_7<beta_c/7 && SINR_VAM_8<beta_c/8
            success_times_9=success_times_9+1;
            SINR_ACK_9=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        end
    end
    if SINR_ACK_1<beta_f && SINR_ACK_2<beta_f && SINR_ACK_3<beta_f && SINR_ACK_4<beta_f && SINR_ACK_5<beta_f && SINR_ACK_6<beta_f && SINR_ACK_7<beta_f && SINR_ACK_8<beta_f && SINR_ACK_9<beta_f && L>=10
        SINR_VAM_10=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture);
        SINR_EVA_10=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi);
        if SINR_VAM_10>=beta_c/10 && SINR_VAM_1<beta_c && SINR_VAM_2<beta_c/2 && SINR_VAM_3<beta_c/3 && SINR_VAM_4<beta_c/4 && SINR_VAM_5<beta_c/5 && SINR_VAM_6<beta_c/6 && SINR_VAM_7<beta_c/7 && SINR_VAM_8<beta_c/8 && SINR_VAM_9<beta_c/9
            success_times_10=success_times_10+1;
        end
    end

    if SINR_EVA_1>=beta_e
        secrecy_times_1=secrecy_times_1-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e
        secrecy_times_2=secrecy_times_2-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e
        secrecy_times_3=secrecy_times_3-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e
        secrecy_times_4=secrecy_times_4-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e
        secrecy_times_5=secrecy_times_5-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e
        secrecy_times_6=secrecy_times_6-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e
        secrecy_times_7=secrecy_times_7-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e || SINR_EVA_8>=beta_e
        secrecy_times_8=secrecy_times_8-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e || SINR_EVA_8>=beta_e || SINR_EVA_9>=beta_e
        secrecy_times_9=secrecy_times_9-1;
    end
    if SINR_EVA_1>=beta_e || SINR_EVA_2>=beta_e || SINR_EVA_3>=beta_e || SINR_EVA_4>=beta_e || SINR_EVA_5>=beta_e || SINR_EVA_6>=beta_e || SINR_EVA_7>=beta_e || SINR_EVA_8>=beta_e || SINR_EVA_9>=beta_e || SINR_EVA_10>=beta_e
        secrecy_times_10=secrecy_times_10-1;
    end
    disp(['Option2: Lambda_V=',num2str(Lambda_V),'; Lambda_U=',num2str(Lambda_U),'; Lambda_R=',num2str(Lambda_R),'*pi; n=',num2str(n),'; simulation_times:',num2str(times)]);
end

S_sim_L_1=(success_times_1)/times_total;
S_sim_L_2=(success_times_1+success_times_2)/times_total;
S_sim_L_3=(success_times_1+success_times_2+success_times_3)/times_total;
S_sim_L_4=(success_times_1+success_times_2+success_times_3+success_times_4)/times_total;
S_sim_L_5=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5)/times_total;
S_sim_L_6=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6)/times_total;
S_sim_L_7=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7)/times_total;
S_sim_L_8=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7+success_times_8)/times_total;
S_sim_L_9=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7+success_times_8+success_times_9)/times_total;
S_sim_L_10=(success_times_1+success_times_2+success_times_3+success_times_4+success_times_5+success_times_6+success_times_7+success_times_8+success_times_9+success_times_10)/times_total;

SU_O2_sim=[S_sim_L_1;S_sim_L_2;S_sim_L_3;S_sim_L_4;S_sim_L_5;S_sim_L_6;S_sim_L_7;S_sim_L_8;S_sim_L_9;S_sim_L_10];

SE_sim_1=secrecy_times_1/times_total;
SE_sim_2=secrecy_times_2/times_total;
SE_sim_3=secrecy_times_3/times_total;
SE_sim_4=secrecy_times_4/times_total;
SE_sim_5=secrecy_times_5/times_total;
SE_sim_6=secrecy_times_6/times_total;
SE_sim_7=secrecy_times_7/times_total;
SE_sim_8=secrecy_times_8/times_total;
SE_sim_9=secrecy_times_9/times_total;
SE_sim_10=secrecy_times_10/times_total;

SE_O2_sim=[SE_sim_1;SE_sim_2;SE_sim_3;SE_sim_4;SE_sim_5;SE_sim_6;SE_sim_7;SE_sim_8;SE_sim_9;SE_sim_10];
end

%% Road Line Generation
function [L_number,L_feture]=Road_Generation(Lambda_R)
L_number=poissrnd(Lambda_R*2*pi)-1;
L_feture=zeros(2,1);
if L_number~=0
    for i=1:L_number
        L_feture(1,i)=rand(1,1);%Distance
        L_feture(2,i)=2*pi*rand(1,1);%Angle
    end
end
end

%% Simulation
function SINR_VAM=Monte_Carlo_VAM_sim(Lambda_V,Lambda_U,xi,L_number,L_feture)
P_u=25;
P_u=10^(P_u/10);
P_v=20;
P_v=10^(P_v/10);

alpha=4;
mu=1;

kappa_u=1/4;
kappa_v=1;

P_interfernce_u=0;
P_interfernce_v=0;
% Typical Road
% Generate RSU
R_number=poissrnd(Lambda_U);
while R_number==0
    R_number=poissrnd(Lambda_U);
end
R_position=zeros(2,R_number);
R_position(1,:)=-0.5+1*rand(1,R_number);
% Receiving RSU
Dis_min=zeros(1,1);
for R_times=1:R_number
    Dis_min(R_times)=sqrt(R_position(1,R_times)^2+R_position(2,R_times)^2);
end
[~,Receiving_R_p]=min(Dis_min);
Receiving_R_position=R_position(:,Receiving_R_p);
R_position(:,Receiving_R_p)=[];
% Interference of RSU on Typical Road
for R_times=1:R_number-1
    if rand(1,1)<=xi*kappa_u
        Dis_R=sqrt((R_position(1,R_times)-Receiving_R_position(1))^2+(R_position(2,R_times)-Receiving_R_position(2))^2);
        h=exprnd(mu);
        P_interfernce_u=P_interfernce_u+P_u*h*Dis_R^(-alpha);
    end
end
% Generate Vehicle
V_number=poissrnd(Lambda_V);
V_position=zeros(2,V_number);
V_position(1,:)=-0.5+1*rand(1,V_number);
% Interference of vehicle on typical road
for V_times=1:V_number
    if rand(1,1)<=xi*kappa_v
        Dis_V=sqrt((V_position(1,V_times)-Receiving_R_position(1))^2+(V_position(2,V_times)-Receiving_R_position(2))^2);
        h=exprnd(mu);
        P_interfernce_v=P_interfernce_v+P_v*h*Dis_V^(-alpha);
    end
end

if L_number~=0
    % Plot Interfering Road Line and Vehicle
    for L_point=1:L_number
        % Plot Road Line
        L_projection_x=L_feture(1,L_point)*cos(L_feture(2,L_point));
        L_projection_y=L_feture(1,L_point)*sin(L_feture(2,L_point));
        L_k=tan(pi/2+L_feture(2,L_point));
        % Road Line Length
        L_x_left=-0.5;
        L_y_left=L_k*(-0.5-L_projection_x)+L_projection_y;
        L_x_right=0.5;
        L_y_right=L_k*(0.5-L_projection_x)+L_projection_y;
        L_length=sqrt((L_x_right-L_x_left)^2+(L_y_right-L_y_left)^2);
        % Generate RSU
        R_number=poissrnd(Lambda_U*L_length);
        R_position=zeros(2,R_number);
        R_position(1,:)=-0.5+1*rand(1,R_number);
        for R_times=1:R_number
            R_position(2,R_times)=L_k*(R_position(1,R_times)-L_projection_x)+L_projection_y;
        end
        % Interference Power of RSU
        for R_times=1:R_number
            if rand(1,1)<=xi*kappa_u
                Dis_R=sqrt((R_position(1,R_times)-Receiving_R_position(1))^2+(R_position(2,R_times)-Receiving_R_position(2))^2);
                h=exprnd(mu);
                P_interfernce_u=P_interfernce_u+P_u*h*Dis_R^(-alpha);
            end
        end
        % Generate Vehicle
        V_number=poissrnd(Lambda_V*L_length);
        V_position=zeros(2,V_number);
        V_position(1,:)=-0.5+1*rand(1,V_number);
        for V_times=1:V_number
            V_position(2,V_times)=L_k*(V_position(1,V_times)-L_projection_x)+L_projection_y;
        end
        % Interference Power of Vehicle
        for V_times=1:V_number
            if rand(1,1)<=xi*kappa_v
                Dis_R=sqrt((V_position(1,V_times)-Receiving_R_position(1))^2+(V_position(2,V_times)-Receiving_R_position(2))^2);
                h=exprnd(mu);
                P_interfernce_v=P_interfernce_v+P_v*h*Dis_R^(-alpha);
            end
        end
    end
end

% VAM Distance
Dis_trans=sqrt(Receiving_R_position(1)^2+Receiving_R_position(2)^2);

% VAM Receiving Power
h=exprnd(mu);
P_receive=P_v*h*Dis_trans^(-alpha);

% SINR Calculation
SINR_VAM=P_receive/(P_interfernce_u+P_interfernce_v);
end

function SINR_NACK=Monte_Carlo_NACK_sim(Lambda_V,Lambda_U,xi,L_number,L_feture)
P_u=25;
P_u=10^(P_u/10);
P_v=20;
P_v=10^(P_v/10);

alpha=4;
mu=1;

kappa_u=1/4;
kappa_v=1;

P_interfernce_R=0;
P_interfernce_V=0;
% Typical Road
% Generate RSU
R_number=poissrnd(Lambda_U);
while R_number==0
    R_number=poissrnd(Lambda_U);
end
R_position=zeros(2,R_number);
R_position(1,:)=-0.5+1*rand(1,R_number);
% Receiving RSU
Dis_min=zeros(1,1);
for R_times=1:R_number
    Dis_min(R_times)=sqrt(R_position(1,R_times)^2+R_position(2,R_times)^2);
end
[~,Receiving_R_p]=min(Dis_min);
Receiving_R_position=R_position(:,Receiving_R_p);
% Interference of RSU on Typical Road
for R_times=1:R_number
    if rand(1,1)<=xi*kappa_u
        Dis_R=sqrt((R_position(1,R_times))^2+(R_position(2,R_times))^2);
        h=exprnd(mu);
        P_interfernce_R=P_interfernce_R+P_u*h*Dis_R^(-alpha);
    end
end
% Generate Vehicle
V_number=poissrnd(Lambda_V);
V_position=zeros(2,V_number);
V_position(1,:)=-0.5+1*rand(1,V_number);
% Interference of vehicle on typical road
for V_times=1:V_number
    if rand(1,1)<=xi*kappa_v
        Dis_V=sqrt((V_position(1,V_times))^2+(V_position(2,V_times))^2);
        h=exprnd(mu);
        P_interfernce_V=P_interfernce_V+P_v*h*Dis_V^(-alpha);
    end
end

% Generate Road Line
if L_number~=0
    % Plot Interfering Road Line and Vehicle
    for L_point=1:L_number
        % Plot Road Line
        L_projection_x=L_feture(1,L_point)*cos(L_feture(2,L_point));
        L_projection_y=L_feture(1,L_point)*sin(L_feture(2,L_point));
        L_k=tan(pi/2+L_feture(2,L_point));
        % Road Line Length
        L_x_left=-0.5;
        L_y_left=L_k*(-0.5-L_projection_x)+L_projection_y;
        L_x_right=0.5;
        L_y_right=L_k*(0.5-L_projection_x)+L_projection_y;
        L_length=sqrt((L_x_right-L_x_left)^2+(L_y_right-L_y_left)^2);
        % Generate RSU
        R_number=poissrnd(Lambda_U*L_length);
        R_position=zeros(2,R_number);
        R_position(1,:)=-0.5+1*rand(1,R_number);
        for R_times=1:R_number
            R_position(2,R_times)=L_k*(R_position(1,R_times)-L_projection_x)+L_projection_y;
        end
        % Interference Power of RSU
        for R_times=1:R_number
            if rand(1,1)<=xi*kappa_u
                Dis_R=sqrt((R_position(1,R_times))^2+(R_position(2,R_times))^2);
                h=exprnd(mu);
                P_interfernce_R=P_interfernce_R+P_u*h*Dis_R^(-alpha);
            end
        end
        % Generate Vehicle
        V_number=poissrnd(Lambda_V*L_length);
        V_position=zeros(2,V_number);
        V_position(1,:)=-0.5+1*rand(1,V_number);
        for V_times=1:V_number
            V_position(2,V_times)=L_k*(V_position(1,V_times)-L_projection_x)+L_projection_y;
        end
        % Interference Power of Vehicle
        for V_times=1:V_number
            if rand(1,1)<=xi*kappa_v
                Dis_R=sqrt((V_position(1,V_times))^2+(V_position(2,V_times))^2);
                h=exprnd(mu);
                P_interfernce_V=P_interfernce_V+P_v*h*Dis_R^(-alpha);
            end
        end
    end
end

% VAM Distance
Dis_trans=sqrt(Receiving_R_position(1)^2+Receiving_R_position(2)^2);

% VAM Receiving Power
h=exprnd(mu);
P_receive=P_u*h*Dis_trans^(-alpha);

% SINR Calculation
SINR_NACK=P_receive/(P_interfernce_R+P_interfernce_V);
end

function SINR_EVA=Monte_Carlo_EVA_sim(Lambda_V,Lambda_U,Lambda_R,xi)
W=0.3;

P_u=25;
P_u=10^(P_u/10);
P_v=20;
P_v=10^(P_v/10);

alpha=4;
mu=1;

kappa_u=1/4;
kappa_v=1;

P_interfernce_U=0;
P_interfernce_V=0;

% Generate Road Line
L_number=poissrnd(Lambda_R*2*pi);
while L_number==0
    L_number=poissrnd(Lambda_R*2*pi);
end
L_feture=zeros(2,L_number);
L_feture(1,:)=rand(1,L_number);%Distance
L_feture(2,:)=2*pi*rand(1,L_number);%Angle

% Plot Interfering Road Line and Vehicle
for L_point=1:L_number
    % Plot Road Line
    L_projection_x=L_feture(1,L_point)*cos(L_feture(2,L_point));
    L_projection_y=L_feture(1,L_point)*sin(L_feture(2,L_point));
    L_k=tan(pi/2+L_feture(2,L_point));
    % Road Line Length
    L_x_left=-0.5;
    L_y_left=L_k*(-0.5-L_projection_x)+L_projection_y;
    L_x_right=0.5;
    L_y_right=L_k*(0.5-L_projection_x)+L_projection_y;
    L_length=sqrt((L_x_right-L_x_left)^2+(L_y_right-L_y_left)^2);
    % Generate RSU
    R_number=poissrnd(xi*kappa_u*Lambda_U*L_length);
    R_position=zeros(2,R_number);
    R_position(1,:)=-0.5+1*rand(1,R_number);
    for R_times=1:R_number
        R_position(2,R_times)=L_k*(R_position(1,R_times)-L_projection_x)+L_projection_y;
    end
    % Interference Power of RSU
    for R_times=1:R_number
        Dis_R=sqrt((R_position(1,R_times))^2+(R_position(2,R_times))^2);
        h=exprnd(mu);
        P_interfernce_U=P_interfernce_U+P_u*h*Dis_R^(-alpha);
    end
    % Generate Vehicle
    V_number=poissrnd(xi*kappa_v*Lambda_V*L_length);
    V_position=zeros(2,V_number);
    V_position(1,:)=-0.5+1*rand(1,V_number);
    for V_times=1:V_number
        V_position(2,V_times)=L_k*(V_position(1,V_times)-L_projection_x)+L_projection_y;
    end
    % Interference Power of Vehicle
    for V_times=1:V_number
        Dis_R=sqrt((V_position(1,V_times))^2+(V_position(2,V_times))^2);
        h=exprnd(mu);
        P_interfernce_V=P_interfernce_V+P_v*h*Dis_R^(-alpha);
    end
end

% VAM Distance
Dis_trans=W*rand(1,1);

% VAM Receiving Power
h=exprnd(mu);
P_receive=P_v*h*Dis_trans^(-alpha);

% SINR Calculation
SINR_EVA=P_receive/(P_interfernce_U+P_interfernce_V);
end

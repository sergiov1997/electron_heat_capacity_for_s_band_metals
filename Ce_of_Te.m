clear variables;

%%% Uncomment this for performing calculations with polylogarithms

% z=-500:1:-0.01;
% z(length(z))=[];
% 
% ply_z=polylog(1.5,z);
% 
% figure(1)
% plot(z,ply_z,'k','linewidth',1.5)
% %hold on
% %plot(z,-(4/3)*pi^(-0.5)*(log(-z)).^1.5,':k','linewidth',1.5)
% set(gca,'FontSize',20) 
% xlabel('$z$','Interpreter','latex','FontSize',25);
% ylabel('$Li_{\frac{3}{2}} (z)$','Interpreter','latex','FontSize',25);
% title('Polilogarithm','Interpreter','latex','FontSize',30)
% 
% 
% figure(2)
% plot(ply_z,z,'k','linewidth',1.5)
% set(gca,'FontSize',20) 
% xlabel('$z$','Interpreter','latex','FontSize',25);
% ylabel('$Li_{\frac{3}{2}^{-1}} (z)$','Interpreter','latex','FontSize',25);
% title('Inverse Polilogarithm','Interpreter','latex','FontSize',30)
% 
% 
% figure(3)
% plot(ply_z,log(-z),'k','linewidth',1.5)
% set(gca,'FontSize',20) 
% xlabel('$z$','Interpreter','latex','FontSize',25);
% ylabel('$ln(-Li_{\frac{3}{2}}^{-1} (z))$','Interpreter','latex','FontSize',25);
% title('ln of Inverse Polilogarithm','Interpreter','latex','FontSize',30)
% 
% 
% figure(4)
% plot(z,ply_z+(4/3)*pi^(-0.5)*(log(-z)).^1.5,'k','linewidth',1.5)
% set(gca,'FontSize',20) 
% xlabel('$z$','Interpreter','latex','FontSize',25);
% ylabel('$Li_{\frac{3}{2}} (z)- \frac{-4*ln(-z)^{\frac{3}{2}}}{3*\sqrt(\pi)}$','Interpreter','latex','FontSize',25);
% title('Difference','Interpreter','latex','FontSize',30)
% 
% 
% figure(5)
% plot(z,ply_z,'k','linewidth',1.5)
% hold on
% plot(z,-(4/3)*pi^(-0.5)*(log(-z)).^1.5,'--k','linewidth',1.5)
% set(gca,'FontSize',20) 
% xlabel('$z$','Interpreter','latex','FontSize',25);
% ylabel('$Li_{\frac{3}{2}} (z) \ or \ \frac{-4*ln(-z)^{\frac{3}{2}}}{3*\sqrt(\pi)}$','Interpreter','latex','FontSize',25);
% legend('$Li_{\frac{3}{2}} (z)$','$\frac{-4*ln(-z)^{\frac{3}{2}}}{3*\sqrt(\pi)}$','Interpreter','latex','FontSize',30)
% title('Polilog and approxiation','Interpreter','latex','FontSize',30)
% 
% 
% myfun = @(x,K0)  polylog(3/2,x)-K0;  % parameterized function
% K0 = -2;                    % parameter
% fun = @(x) myfun(x,K0);    % function of x alone
% 
% x = fzero(fun,-5);
% 
% 
% %%% We plot the argument of the inverse polilogarithm for the range of
% %%% temperatures that we have
% 
% T_el=250:1:1e4;
% arg=-17538857.55./(T_el.^1.5);
% 
% figure(6)
% plot(T_el,arg,'k','linewidth',1.5)
% set(gca,'FontSize',25) 
% xlabel('$T_{el}$','Interpreter','latex','FontSize',25);
% ylabel('$Arg$','Interpreter','latex','FontSize',25);
% title('Argument of inverse polylog','Interpreter','latex','FontSize',25)


%%% We now calculate the electron heat capacity as a function of T_e for
%%% the chosen metal

%%% We first select the metal and option plots

prompt = "Enter a number to select the metal---> 1: Aluminium; 2: Copper; 3: Silver; 4:Gold ; 5: Silicon; 6:Ruthenium; 7:Pt; 8:Thungsten  ";
I_met = input(prompt);


promp_3="Do you want to plot comparison for published data for D-band metals? 1---> Yes; 2---> No; ";
comp_d=input(promp_3);

promp_1598="Do you want to plot plot functions F(x) and M(x) and their derivatives? (Note that Symbolic Math Toolbox needs to be installed): 1---> Yes; 2---> No; ";
option_F_M=input(promp_1598);

promp_2014="Do you want calculate the analytical solution for \mu and Ce? (Note that Symbolic Math Toolbox needs to be installed): 1---> Yes; 2---> No; ";
option_mu_analytical=input(promp_2014);

N_Points=600;

if I_met==1
    ne=18.1e28; %Atom density of Al, in m^{-3}
    Ae=2e7; % in K^{-2} s^{-1} 
    Bl=7e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Al.
    prompt2 = "Aluminium";

elseif I_met==2

    ne=8.45e28; %Atom density of Cu, in m^{-3}
    Ae=0.12e6; % in K^{-2} s^{-1} 
    Bl=3.35e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Cu.
    prompt2 = "Copper";

elseif I_met==3

    ne=5.85e28; %Atom density of Ag, in m^{-3}
    Ae=0.932e7; % in K^{-2} s^{-1} 
    Bl=1.02e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Ag.
    prompt2 = "Silver";

elseif I_met==4

    ne=5.90e28; %Atom density of Au, in m^{-3}
    Ae=1.2e7; % in K^{-2} s^{-1} 
    Bl=1.23e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Au.
    prompt2 = "Gold";

elseif I_met==5

    ne=9.9877e28; %Atom density of Si, in m^{-3}
    Ae=1.2e7; % in K^{-2} s^{-1} 
    Bl=1.23e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Au. 
    %%% THOSE ARE THE PARAMETERS OF AU, ONLY NEEDDED FOR THERMAL
    %%% CONDUCTIVITY, CHANGE THEM FOR THIS METAL IN CASE THE USER WANTS lambda_e  
    prompt2 = "Silicon";

elseif I_met==6

    ne=8.2223e28; %Atom density of Ru, in m^{-3}
    Ae=1.2e7; % in K^{-2} s^{-1} 
    Bl=1.23e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Au. 
    %%% THOSE ARE THE PARAMETERS OF AU, ONLY NEEDDED FOR THERMAL
    %%% CONDUCTIVITY, CHANGE THEM FOR THIS METAL IN CASE THE USER WANTS lambda_e      prompt2 = "Ruthenium ";

elseif I_met==7

    ne=6.6213e28; %Atom density of Ru, in m^{-3}
    Ae=1.2e7; % in K^{-2} s^{-1} 
    %%% THOSE ARE THE PARAMETERS OF AU, ONLY NEEDDED FOR THERMAL
    %%% CONDUCTIVITY, CHANGE THEM FOR THIS METAL IN CASE THE USER WANTS lambda_e      %%% THOSE ARE WRONG IN SI
    prompt2 = "Platinum ";

 elseif I_met==8

    ne=6.305e28; %Atom density of Ru, in m^{-3}
    Ae=1.2e7; % in K^{-2} s^{-1} 
    Bl=1.23e11; % in K^{-1} s^{-1}, parameters for calculating the electron scattering time for Au. 
    %%% THOSE ARE THE PARAMETERS OF AU, ONLY NEEDDED FOR THERMAL
    %%% CONDUCTIVITY, CHANGE THEM FOR THIS METAL IN CASE THE USER WANTS lambda_e      prompt2 = "Tungsten ";


end

h=6.63e-34; hbar=h/(2*pi); me=9.10938e-31; k_B=1.38065e-23;  %phisical constants, in SI units
E_F=(3*pi^2*ne)^(2/3)*hbar^2/(2*me);   %Fermi energy, in J
v_F=sqrt(2*E_F/me);  %Fermi velocity, in m/s
mu_min=-1.25*1.60217e-19+E_F; mu_max=E_F; mu=mu_min:(mu_max-mu_min)/N_Points:mu_max; %Ranges of potential to be used
gamma=pi^2*k_B^2*ne/(2*E_F);

T_00=(4/(3*sqrt(pi)))^(2/3)*E_F/k_B;   %%%Scaling temperature T_0, it appears in our analytical solutions

%%% We now calculate the chemical potential of the electron as a function
%%% of temperature using our aproximation and the quadratic model widely
%%% used

FIRST_TERM=-4*mu.^1.5/(3*sqrt(pi))+2*ne*h^3/(4*pi*sqrt(pi)*(2*me)^1.5);
T_e_ours=(1/k_B)*(FIRST_TERM.*mu.^0.5*(6*sqrt(pi)/pi^2)).^0.5;  %Our solution
T_e_polish=((E_F-mu)*12*E_F/(pi^2*k_B^2)).^0.5;                 %Solution from Majchrak et al paper[2]

figure(53)
plot(mu/(1.60217e-19),T_e_ours,'k','linewidth',1.5)
xlabel('$z$','Interpreter','latex','FontSize',25);
hold on
plot(mu/(1.60217e-19),T_e_polish,'--k','linewidth',1.5)
set(gca,'FontSize',30) 
ylabel('$T_e [K]$','Interpreter','latex','FontSize',30);
xlabel('$\mu [eV]$','Interpreter','latex','FontSize',30);
legend('Our result', 'Majchrzak & Poteralska paper')
title(['Chemical potential for ' prompt2])

% 
% %%% Uncomment for poltting the inverse of the polylogarithm
% 
% XXX=-20:0.05:-0.001;
% YYY=zeros(1,length(XXX));
% 
% for III = 1:length(XXX)
%     YYY(III)=inverting_polylogarithm_3_2(XXX(III));
%     COUNTER=III;
%     disp([' Progress is of ' num2str(100*III/length(XXX)) ' %'])
% end
% 
% figure(7)
% plot(XXX,YYY,'k','linewidth',1.5)
% set(gca,'FontSize',25) 
% xlabel('$x$','Interpreter','latex','FontSize',25);
% ylabel('$Li_{\frac{3}{2}}^{-1} (x)$','Interpreter','latex','FontSize',25);


TT_el=200:10:5e4;  %Ranges of electron temperatures, in K
mu_full=zeros(1,length(TT_el)); U_e=mu_full;
if option_mu_analytical==1
mu_full_dash=mu_full;
Ce_long_analytical=mu_full;
Ce_long_short_analytical=Ce_long_analytical;
end
C_e_alfull=zeros(1,length(TT_el)-2);
TTT_el=TT_el; TTT_el(1)=[]; TTT_el(length(TTT_el))=[]; 

%%% Here, we compute the Taylor expansion of the chemical potential and for
%%% the total internal electron energy from the analytical solution

TTe_Taylor_quadr=200:1000:5e4;
COEFFICIENT_SECOND_DERIVATIVE=-1.36088005; %%% Coefficient obtained by numerically calculating the second derivative in 0 of later defifined F(x)
COEFFICIENT_FOURTH_DERIVATIVE=-16.5661;     %%% Coefficient obtained by numerically calculating the fourth derivative in 0 of later defifined F(x)

mu_Taylor=zeros(1,length(TTe_Taylor_quadr));
mu_Taylor_o_4=mu_Taylor;

M_of_0=-0.483603;                         %%% Value of later defined function M at 0
COEFFICIENT_SECOND_DERIVATIVE_M=-2.72101; %%% Coefficient obtained by numerically calculating the second derivative in 0 of later defifined M(x)
COEFFICIENT_FOURTH_DERIVATIVE_M=33.11007;    %%% Coefficient obtained by numerically calculating the fourth derivative in 0 of later defifined M(x)


Ue_Taylor_O_2=zeros(1,length(TT_el));
Ue_Taylor_O_4=mu_Taylor_o_4;

Ce_Taylor_O_3=mu_Taylor_o_4;
Ce_Taylor_O_5=Ce_Taylor_O_3;

beta=4*pi*(2*me/h^2)^1.5*(-3*sqrt(pi)/4)*k_B^2.5*(1/6)*T_00^(-1.5)*COEFFICIENT_FOURTH_DERIVATIVE_M; %%% Coefficient beta froom the cubic approximaiion C_e = \gamma T_e + \beta T_e^3

for I=1:length(TTe_Taylor_quadr)
    mu_Taylor(I)=E_F+(1/2)*COEFFICIENT_SECOND_DERIVATIVE*(3*sqrt(pi)/4)^(2/3)*(k_B^2/E_F)*TTe_Taylor_quadr(I)^2;
    mu_Taylor_o_4(I)=mu_Taylor(I)+k_B/(24*T_00^3)*COEFFICIENT_FOURTH_DERIVATIVE*TTe_Taylor_quadr(I)^4;   %%%Taylor approximation to \mu at order two and four

    Ue_Taylor_O_4(I)=4*pi*(2*me/h^2)^1.5*(-3*sqrt(pi)/4)*k_B^2.5*(T_00^2.5*M_of_0+0.5*T_00^0.5*COEFFICIENT_SECOND_DERIVATIVE_M*TTe_Taylor_quadr(I)^2+(1/24)*(T_00^(-1.5))*COEFFICIENT_FOURTH_DERIVATIVE_M*TTe_Taylor_quadr(I)^4); %%%Taylor approximation to U_e at order four
    Ce_Taylor_O_3(I)=4*pi*(2*me/h^2)^1.5*(-3*sqrt(pi)/4)*k_B^2.5*(T_00^0.5*COEFFICIENT_SECOND_DERIVATIVE_M*TTe_Taylor_quadr(I)+(1/6)*T_00^(-1.5)*COEFFICIENT_FOURTH_DERIVATIVE_M*TTe_Taylor_quadr(I)^3);   %%% C_e from differenciating the fourth order Taylor expansion for U_e
    Ce_Taylor_O_5(I)=Ce_Taylor_O_3(I)+(9.92e-19)*TTe_Taylor_quadr(I)^5;
end

for III=1:length(TT_el)
   Ue_Taylor_O_2(III)=4*pi*(2*me/h^2)^1.5*(-3*sqrt(pi)/4)*k_B^2.5*(T_00^2.5*M_of_0+0.5*T_00^0.5*COEFFICIENT_SECOND_DERIVATIVE_M*TT_el(III)^2); %%%Taylor approximation to U_e at order two 

end


%%% We now calculate our exact, analytical solution for the chemical
%%% potential

for III=1:length(TT_el)


    Argument=-4*E_F^(1.5)/(3*sqrt(pi)*(k_B*TT_el(III))^1.5);

    mu_full(III)=k_B*TT_el(III)*log(-inverting_polylogarithm_3_2(Argument));

    if option_mu_analytical==1

    mu_full_dash(III)=k_B*log(-inverting_polylogarithm_3_2(-(T_00/TT_el(III))^1.5))+(3/2)*k_B*(T_00/TT_el(III))^1.5*1/(polylog(1/2,inverting_polylogarithm_3_2(-(T_00/TT_el(III))^1.5)));  %analytically found derivative of mu
    Ce_long_analytical(III)=-3*pi^1.5*(2*me/h^2)^1.5*k_B^2.5*(2.5*TT_el(III)^1.5*polylog(5/2,-exp(mu_full(III)/(k_B*TT_el(III))))+TT_el(III)^2.5*(mu_full_dash(III)*k_B*TT_el(III)-mu_full(III)*k_B)/(k_B*TT_el(III))^2*polylog(3/2,-exp(mu_full(III)/(k_B*TT_el(III)))));
    Ce_long_short_analytical(III)=-3*pi^1.5*(2*me/h^2)^1.5*k_B^2.5*(2.5*TT_el(III)^1.5*polylog(5/2,inverting_polylogarithm_3_2(-(T_00/TT_el(III))^1.5))-1.5*TT_el(III)^2.5*(T_00^3/(TT_el(III))^4)/polylog(1/2,inverting_polylogarithm_3_2(-(T_00/TT_el(III))^1.5)));
    %%% The solution Ce_long_short_analytical has mu(Te) substituted
    %%% already (also appears in the new version of the paper)
    end

    disp([' Progress is of ' num2str(100*III/length(TT_el)) ' %'])

end

if option_mu_analytical==1
    figure(2014)
    plot(TT_el,Ce_long_analytical,'k','linewidth',1.5)
    hold on
    plot(250:10:4e4, gamma*(250:10:4e4),'--k','linewidth',1.5)
    xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
    ylabel('$C_e [J m^{-3} K^{-1}]$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',30) 
    legend('Full analytical solution','First order approximation $\gamma \cdot T_e$','Interpreter','latex','FontSize',30)
    sgtitle(' Analytical solution of $ C_e (T_e) $','Interpreter','latex','FontSize',30);
    
    figure(2015)
    plot(TT_el,Ce_long_short_analytical,'k','linewidth',1.5)
    hold on
    plot(250:10:4e4, gamma*(250:10:4e4),'--k','linewidth',1.5)
    xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
    ylabel('$C_e [J m^{-3} K^{-1}]$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',30) 
    legend('Full analytical solution','First order approximation $\gamma \cdot T_e$','Interpreter','latex','FontSize',30)
    sgtitle(' Analytical solution of $ C_e (T_e) $','Interpreter','latex','FontSize',30);

    figure(2016)
    plot(TT_el/T_00,Ce_long_short_analytical/T_00^1.5,'k','linewidth',1.5)
    xlabel('$\frac{T_e}{T_0}$','Interpreter','latex','FontSize',30);
    ylabel('$\frac{C_e}{T_0^{\frac{3}{2}}} [J m^{-3} K^{-\frac{5}{2}}]$','Interpreter','latex','FontSize',30);
    set(gca,'FontSize',30) 
    %legend('Full analytical solution','First order approximation $\gamma \cdot T_e$','Interpreter','latex','FontSize',30)
    sgtitle(' Analytical solution of $ C_e (T_e) $ for all s-band metals','Interpreter','latex','FontSize',30);
end


figure(8)
plot(TT_el,mu_full/(1.60217e-19),'k','linewidth',1.5)
hold on
plot(T_e_ours,mu/(1.60217e-19),'--k','linewidth',1.5)
xlabel('$z$','Interpreter','latex','FontSize',25);
hold on
plot(T_e_polish,mu/(1.60217e-19),':k','linewidth',1.5)
set(gca,'FontSize',30) 
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$\mu [eV]$','Interpreter','latex','FontSize',30);
legend('Full analytical solution','Our expansion for low temperatures', 'Quadratic approximation')
title([' ' prompt2])

dTe_el=TT_el(2)-TT_el(1);

%%% We calculate the total internal energy of the electrons; the electron
%%% heat capacity is given by the derivative of this total internal energy
%%% wih respect to temperature

for III=1:length(TT_el)


    U_e(III)=4*pi*(2*me/h^2)^1.5*(-3*sqrt(pi)/4)*(k_B*TT_el(III))^2.5*incomplete_polylog(5/2,0,-exp(mu_full(III)/(k_B*TT_el(III))));

    disp([' Progress 2 is of ' num2str(100*III/length(TT_el)) ' %'])

end



%%% We now derive numerically the total internal energy in order to obtain
%%% the heat capacity

for II=2:length(TT_el)-1


    C_e_alfull(II-1)=(U_e(II+1)-U_e(II-1))/(2*dTe_el);

    

end



figure(9)
plot(TT_el,U_e,'k','linewidth',1.5)
hold on
plot(TT_el,Ue_Taylor_O_2,':k','linewidth',1.5,'MarkerSize', 6)
hold on
plot(TTe_Taylor_quadr,Ue_Taylor_O_4,'square k','linewidth',1.5,'MarkerSize', 6)
set(gca,'FontSize',30) 
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$U_e [J]$','Interpreter','latex','FontSize',30);
legend('Exact analytical solution for $U_e$', 'Taylor expansion for $U_e$ at $\mathcal{O} (T_e^2)$','Taylor expansion for $U_e$ at $\mathcal{O} (T_e^4)$','Interpreter','latex','FontSize',30)
title(['Total internal energy for ' prompt2])

Virginia_Ce_Data_Al=readmatrix('Datos_Ce_comparacion.xlsx','Range','D3:E1002');  %We read data for Al from excel file

figure(10)
plot(TTT_el, C_e_alfull,'k','linewidth',1.5)
hold on
plot(250:10:4e4, gamma*(250:10:4e4),'--k','linewidth',1.5)

%%% Uncomment for comparing with data obtained from https://compmat.org/electron-phonon-coupling/
if I_met==1
hold on  %%% Published data from Excel file in the case of aluminion
plot(Virginia_Ce_Data_Al(:,1), Virginia_Ce_Data_Al(:,2),':k','linewidth',1.5)
xlim([0 5.2e4])
end

set(gca,'FontSize',30) 
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$C_e [J m^{-3} K^{-1}]$','Interpreter','latex','FontSize',30);
legend('Full analytical solution','First order approximation $\gamma \cdot T_e$','Published Data','Interpreter','latex','FontSize',30)
title(['Electron specific heat of ' prompt2])

%% Uncomming for plotting the electron conductivity for the metal


TTT_la=250:10:5e3;
lambda_e=zeros(length(TTT_la),length(TTT_el)); lambda_e_approx=lambda_e;
Diff_per_labda=lambda_e; lambda_e_b_approx=lambda_e_approx;
Diff_per_labda_b=Diff_per_labda;

%%%We now calculate the electron thermal conductivity 

for I=1:length(TTT_el)
for J=1:length(TTT_la)
    
    tau_ee=1/(Ae*(TTT_el(I))^2+Bl*TTT_la(J)); %%% Electron relaxation time  
    lambda_e(J,I)=(1/3)*C_e_alfull(I)*v_F^2*tau_ee; %%% Using the exact solution for Ce
    lambda_e_approx(J,I)=lambda_e(1,1)*TTT_el(I)/TTT_la(J); %%% The approximation widely used for lambda e
    Diff_per_labda(J,I)=100*(lambda_e_approx(J,I)-lambda_e(J,I))/lambda_e(J,I);
    lambda_e_b_approx(J,I)=(1/3)*v_F^2*(gamma*TTT_el(I)+beta*TTT_el(I)^3)/(Ae*TTT_el(I)^2+Bl*TTT_la(J));  %%% Approximation by substituting cubic approx for C_e
    Diff_per_labda_b(J,I)=100*(lambda_e_b_approx(J,I)-lambda_e(J,I))/lambda_e(J,I);
end
end


%[T_el,T_la]=meshgrid(TTT_el,TTT_la);

figure(11)


%subplot(1,2,1)
surf(TTT_el,TTT_la,lambda_e,'EdgeColor','red')

hold on

surf(TTT_el,TTT_la,lambda_e_approx,'EdgeColor','blue')
%colorbar
shading interp
set(gca,'FontSize',30) 
ax = gca;
ax.YAxis.Exponent = 3;
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$T_l [K]$','Interpreter','latex','FontSize',30);
zlabel('$\lambda_e [W/(m \cdot K)]$','Interpreter','latex','FontSize',30);
%legend('Full analytical solution','First order approximation $\gamma \cdot T_e$', 'Interpreter','latex','FontSize',30)
title(['Thermal conductivity of '  prompt2])
% subplot(1,2,2)
% pcolor(TTT_el,TTT_la,lambda_e)
% colorbar
% shading interp
% set(gca,'FontSize',30) 
% xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
% ylabel('$T_l [K]$','Interpreter','latex','FontSize',30);
% zlabel('$\lambda_e [W/(m \cdot K)]$','Interpreter','latex','FontSize',30);
% % legend('Full analytical solution','First order approximation $\gamma \cdot T_e$', 'Interpreter','latex','FontSize',30)
% title('Thermal conductivity of Cu')

% figure(575)
% MAX_LAMBDA_E= max([max(max(lambda_e)) max(max(lambda_e_approx))]);
% 
% subplot(1,2,1)
% pcolor(TTT_el,TTT_la,lambda_e)
% shading interp
% clim([0 MAX_LAMBDA_E])
% set(gca,'FontSize',30) 
% ax = gca;
% ax.YAxis.Exponent = 3;
% colorbar
% xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
% ylabel('$T_l [K]$','Interpreter','latex','FontSize',30);
% zlabel('$\lambda_e [W/(m \cdot K)]$','Interpreter','latex','FontSize',30);
% title('Exact solution')
% hold on
% subplot(1,2,2)
% pcolor(TTT_el,TTT_la,lambda_e_approx)
% clim([0 MAX_LAMBDA_E])
% %colorbar
% shading interp
% set(gca,'FontSize',30) 
% ax = gca;
% ax.YAxis.Exponent = 3;
% xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
% ylabel('$T_l [K]$','Interpreter','latex','FontSize',30);
% zlabel('$\lambda_e [W/(m \cdot K)]$','Interpreter','latex','FontSize',30);
% %legend('Full analytical solution','First order approximation $\gamma \cdot T_e$', 'Interpreter','latex','FontSize',30)
% title(' $\lambda_e = \lambda_0 \frac{T_e}{T_l}$','Interpreter','latex','FontSize',30);
% colorbar
% %sgtitle(['Thermal conductivity of '  prompt2],'FontSize',30)
% sgtitle(['Thermal conductivity $[W/(m \cdot K)]$ of '  prompt2 ],'Interpreter','latex','FontSize',30)


figure(576)
pcolor(TTT_el,TTT_la,Diff_per_labda)
shading interp
clim([0 150])
set(gca,'FontSize',30) 
%ax = gca;
%ax.YAxis.Exponent = 3;
colorbar
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$T_l [K]$','Interpreter','latex','FontSize',30);
zlabel('$\lambda_e [W/(m \cdot K)]$','Interpreter','latex','FontSize',30);
title(['Percentaje error using the approximation  $\lambda_e = \lambda_0 \frac{T_e}{T_l}$' prompt2],'Interpreter','latex','FontSize',30)

figure(577)
pcolor(TTT_el,TTT_la,Diff_per_labda_b)
shading interp
%clim([0 150])
set(gca,'FontSize',30) 
ax = gca;
%ax.YAxis.Exponent = 3;
colorbar
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$T_l [K]$','Interpreter','latex','FontSize',30);
zlabel('$\lambda_e [W/(m \cdot K)]$','Interpreter','latex','FontSize',30);
title(['Percentaje error using the approximation  $\lambda_e = \frac{1}{3} v_F^2 \frac{\gamma T_e + \beta T_e^3}{A_e T_e^2 + B_l T_l }$' prompt2],'Interpreter','latex','FontSize',30)



% %%%Uncomment for plotting the function (k_B T_e)^1.5+Li_{3/2} (-exp(\mu / (k_B T_e)))
% RHS_fUNCTION=zeros(1,length(TT_el));
% for I=1:length(TT_el)
% RHS_fUNCTION(I)=(k_B*TT_el(I))^1.5*incomplete_polylog(1.5,0,-exp(E_F/(k_B*TT_el(I))));
% end
% figure(1234)
% plot(TT_el,((-3*sqrt(pi)/4)*RHS_fUNCTION)^(2/3),'k','linewidth',1.5)
% set(gca,'FontSize',30) 
% xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
% ylabel('$\left( k_B T_e \right)^{\frac{3}{2}} \cdot Li_{\frac{3}{2}}$','Interpreter','latex','FontSize',30);


%bb=inverting_polylogarithm_3_2(-2)


% function inv_Lit=inverting_polylogarithm_3_2_from_table(Y)
% if Y<-13300
%     inv_Lit=inverting_polylogarithm_3_2(Y);
% else
%     load yyy.mat y
% 
%     INDEX=floor((Y+13300)/0.1);
%     inv_Lit=y(INDEX);
% end
% end

if option_F_M==1
%%% Uncomment for plotting and operating with the function F(x)= x ln \left[ - Li_{\frac{3}{2}}^{-1} \left( - 1/x^(3/2) \right) \right]
%%% Note that for plottic M of x Symbolic Math Toolbox needs to be
%%% installed


x=0.002:0.002:1;  %%% x can be 0.002 if we´re only diff. numerically
dx=x(2)-x(1);

x_reduced=x(2:length(x)-1);
x_reduced_2=x(3:length(x)-2);
x_reduced_3=x(4:length(x)-3);
x_reduced_4=x(5:length(x)-4);



YY=zeros(1,length(x)); 
DF_OVER_DX=zeros(1,length(x)-2);
D2F_OVER_DX2=zeros(1,length(x)-4);
D3F_OVER_DX3=zeros(1,length(x)-6);
D4F_OVER_DX4=zeros(1,length(x)-8);


for I=1:length(x) %function F(x)
YY(I)=x(I)*log(-inverting_polylogarithm_3_2(-x(I)^(-3/2)));
if mod(I,10)==0
     disp([' Progress F(x) is of ' num2str(100*I/length(x)) ' %'])
end
end

for I=2:length(x)-1  %%% First derivative of F(x)

DF_OVER_DX(I-1)=(YY(I+1)-YY(I-1))/(2*dx);

end

for I=2:length(DF_OVER_DX)-1 %%% Second derivative of F(x)

D2F_OVER_DX2(I-1)=(DF_OVER_DX(I+1)-DF_OVER_DX(I-1))/(2*dx);

end

for I=2:length(D2F_OVER_DX2)-1 %%% Third derivative of F(x)

D3F_OVER_DX3(I-1)=(D2F_OVER_DX2(I+1)-D2F_OVER_DX2(I-1))/(2*dx);

end

for I=2:length(D3F_OVER_DX3)-1 %%% Fourth derivative of F(x)

D4F_OVER_DX4(I-1)=(D3F_OVER_DX3(I+1)-D3F_OVER_DX3(I-1))/(2*dx);

end

figure(1944)   %We plot the function F(x)
plot(x,YY,'k','linewidth',1.5)
set(gca,'FontSize',30) 
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ x \cdot ln \left[ - Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right]$','Interpreter','latex','FontSize',30);
title(' Adimensional function $ x \cdot ln \left[ - Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}}  \right) \right]$','Interpreter','latex','FontSize',30);

figure(1945)   %We plot the derivatives of F(X)

subplot(2,2,1)
plot(x_reduced,DF_OVER_DX,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d F(x)}{dx}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(2,2,2)
plot(x_reduced_2,D2F_OVER_DX2,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^2 F(x)}{dx^2}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(2,2,3)
plot(x_reduced_3,D3F_OVER_DX3,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^3 F(x)}{dx^3}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(2,2,4)
plot(x_reduced_4,D4F_OVER_DX4,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^4 F(x)}{dx^4}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
sgtitle(' Derivatives of $ x \cdot ln \left[ - Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}}  \right) \right]$','Interpreter','latex','FontSize',30);


%%% Uncomment for plotting and operating with the function M(x)= x^{\frac{5}{2}} Li_{\frac{5}{2}} \left(  Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right)

MM=YY;
DMM_OVER_DX=zeros(1,length(x)-2);
D2MM_OVER_DX2=zeros(1,length(x)-4);
D3MM_OVER_DX3=zeros(1,length(x)-6);
D4MM_OVER_DX4=zeros(1,length(x)-8);

%%% Uncomment for plotting the analytically calculated derivatives of M(x)
AN_DMM_OVER_DX=zeros(1,length(x));

% %%% Second derivative of M terms
% AN_D2MM_OVER_DX2=AN_DMM_OVER_DX;
% FIRST_TERM=AN_DMM_OVER_DX;
% SECOND_TERM=FIRST_TERM;
% THIRD_TERM=FIRST_TERM;
% TRial=AN_DMM_OVER_DX;
% 
% g_of_x=TRial; g_dash_of_x=TRial; FRIST_NUM=TRial; SECOND_NUM=TRial; 

for I=1:length(x)
    AN_DMM_OVER_DX(I)=(5/2)*x(I)^(3/2)*polylog(5/2,inverting_polylogarithm_3_2(-x(I)^(-3/2)))-(3/2)*(x(I)^(-3/2))*1/(polylog(1/2,inverting_polylogarithm_3_2(-x(I)^(-3/2))));
%     FIRST_TERM(I)=(15/4)*x(I)^(1/2)*polylog(5/2,inverting_polylogarithm_3_2(-x(I)^(-3/2)));
%     SECOND_TERM(I)=-(3/2)*x(I)^(-5/2)*(polylog(1/2,inverting_polylogarithm_3_2(-x(I)^(-3/2))))^(-1);
%     THIRD_TERM(I)=-(9/4)*x(I)^(-4)*polylog(-1/2,inverting_polylogarithm_3_2(-x(I)^(-3/2)))*(polylog(1/2,inverting_polylogarithm_3_2(-x(I)^(-3/2))))^(-3);
%     AN_D2MM_OVER_DX2(I)=FIRST_TERM(I)+SECOND_TERM(I)+THIRD_TERM(I);

%     g_of_x(I)=inverting_polylogarithm_3_2(-x(I)^(-3/2));
%     g_dash_of_x(I)=(3/2)*x(I)^(-5/2)*g_of_x(I)/polylog(1/2,g_of_x(I));
%     FRIST_NUM(I)=g_dash_of_x(I)*((5/2)*x(I)^4*polylog(3/2,g_of_x(I))+(3/2)*x(I)*polylog(-1/2,g_of_x(I))*(polylog(1/2,g_of_x(I)))^(-2));
%     SECOND_NUM(I)=(15/4)*x(I)^3*(polylog(5/2,g_of_x(I)))+9/(4*(polylog(1/2,g_of_x(I))));
% 
%     AN_D2MM_OVER_DX2(I)=FRIST_NUM(I)*x(I)^(-2.5)/g_dash_of_x(I)+SECOND_NUM(I)*x(I)^(-2.5);    


   % TRial(I)=(5/2)*x(I)^(3/2)*polylog(5/2,inverting_polylogarithm_3_2(-x(I)^(-3/2)));
    if mod(I,10)==0
     disp([' Progress M`(x) is of ' num2str(100*I/length(x)) ' %'])
    end
end

figure(13313)
plot(x,AN_DMM_OVER_DX,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d M(x)}{dx}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
sgtitle('Analytical Derivatives of $ M(x)= x^{\frac{5}{2}} Li_{\frac{5}{2}} \left(  Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right) $','Interpreter','latex','FontSize',30);


% %%% Uncomment for trial plots of the analytically calculated second
% %%% derivative of M
% figure(120614)
% plot(x,AN_D2MM_OVER_DX2,'k','linewidth',1.5)
% xlabel('$ x $','Interpreter','latex','FontSize',30);
% ylabel('$ Term$','Interpreter','latex','FontSize',30);
% set(gca,'FontSize',20) 


for I=1:length(x)
MM(I)=x(I)^(5/2)*polylog(5/2,inverting_polylogarithm_3_2(-x(I)^(-3/2)));
if mod(I,10)==0
     disp([' Progress M(x) is of ' num2str(100*I/length(x)) ' %'])
end
end


for I=2:length(x)-1  %%% First derivative of M(x)

DMM_OVER_DX(I-1)=(MM(I+1)-MM(I-1))/(2*dx);

end

for I=2:length(DMM_OVER_DX)-1 %%% Second derivative of M(x)

D2MM_OVER_DX2(I-1)=(DMM_OVER_DX(I+1)-DMM_OVER_DX(I-1))/(2*dx);

end

for I=2:length(D2MM_OVER_DX2)-1 %%% Third derivative of M(x)

D3MM_OVER_DX3(I-1)=(D2MM_OVER_DX2(I+1)-D2MM_OVER_DX2(I-1))/(2*dx);

end

for I=2:length(D3MM_OVER_DX3)-1 %%% Fourth derivative of M(x)

D4MM_OVER_DX4(I-1)=(D3MM_OVER_DX3(I+1)-D3MM_OVER_DX3(I-1))/(2*dx);

end

%%% Uncomment for calculating also the fiffth and sixth derivatives of M(x)

D5MM_OVER_DX5=zeros(1,length(x)-10);
D6MM_OVER_DX6=zeros(1,length(x)-12);
x_reduced_5=x(6:length(x)-5);
x_reduced_6=x(7:length(x)-6);


for I=2:length(D4MM_OVER_DX4)-1 %%% Fourth derivative of M(x)

D5MM_OVER_DX5(I-1)=(D4MM_OVER_DX4(I+1)-D4MM_OVER_DX4(I-1))/(2*dx);

end

for I=2:length(D5MM_OVER_DX5)-1 %%% Fourth derivative of M(x)

D6MM_OVER_DX6(I-1)=(D5MM_OVER_DX5(I+1)-D5MM_OVER_DX5(I-1))/(2*dx);

end



figure(666)
plot(x,MM,'k','linewidth',1.5)
set(gca,'FontSize',30) 
xlabel('$x$','Interpreter','latex','FontSize',30);
ylabel('$x^{\frac{5}{2}} Li_{\frac{5}{2}} \left(  Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right)$','Interpreter','latex','FontSize',30);
title('Function $x^{\frac{5}{2}} Li_{\frac{5}{2}} \left( Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right)$','Interpreter','latex','FontSize',30);

figure(6666)
subplot(2,2,1)
plot(x_reduced,DMM_OVER_DX,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d M(x)}{dx}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(2,2,2)
plot(x_reduced_2,D2MM_OVER_DX2,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^2 M(x)}{dx^2}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(2,2,3)
plot(x_reduced_3,D3MM_OVER_DX3,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^3 M(x)}{dx^3}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(2,2,4)
plot(x_reduced_4,D4MM_OVER_DX4,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^4 M(x)}{dx^4}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
sgtitle(' Derivatives of $ M(x)= x^{\frac{5}{2}} Li_{\frac{5}{2}} \left(  Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right) $','Interpreter','latex','FontSize',30);

figure(23323)
subplot(1,2,1)
plot(x_reduced_5,D5MM_OVER_DX5,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^5 M(x)}{dx^5}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
subplot(1,2,2)
plot(x_reduced_6,D6MM_OVER_DX6,'k','linewidth',1.5)
xlabel('$ x $','Interpreter','latex','FontSize',30);
ylabel('$ \frac{d^6 M(x)}{dx^6}$','Interpreter','latex','FontSize',30);
set(gca,'FontSize',20) 
sgtitle(' Derivatives of $ M(x)= x^{\frac{5}{2}} Li_{\frac{5}{2}} \left(  Li_{\frac{3}{2}}^{-1} \left( - \frac{1}{x^{\frac{3}{2}}} \right) \right) $','Interpreter','latex','FontSize',30);

end

% mu_full ;TT_el;  dTT_el=TT_el(2)-TT_el(1);
% 
% Dmu_OVER_DT=zeros(1,length(TT_el)-2);
% Dmu2_OVER_DT2=zeros(1,length(TT_el)-4);
% 
% TTel_reduced=TT_el(2:length(TT_el)-1);
% TTel_reduced_2=TT_el(3:length(TT_el)-2);
% mu_Expansion_from_full=zeros(1,length(TTel_reduced_2));
% 
% 
% 
% for I=2:length(TT_el)-1  %%% First derivative of F(x)
% 
% Dmu_OVER_DT(I-1)=(mu_full(I+1)-mu_full(I-1))/(2*dTT_el);
% 
% end
% 
% for I=2:length(Dmu_OVER_DT)-1 %%% Second derivative of F(x)
% 
% Dmu2_OVER_DT2(I-1)=(Dmu_OVER_DT(I+1)-Dmu_OVER_DT(I-1))/(2*dx);
% 
% end
% 
% for I=1:length(TTel_reduced_2)
% mu_Expansion_from_full(I)=E_F+Dmu2_OVER_DT2(1)*TTel_reduced_2(I)^2;
% end

figure(1776)
plot(TT_el,mu_full/(1.60217e-19),'k','linewidth',1.5)
hold on
plot(T_e_ours,mu/(1.60217e-19),'--k','linewidth',1.5)
xlabel('$z$','Interpreter','latex','FontSize',25);
hold on
plot(T_e_polish,mu/(1.60217e-19),':k','linewidth',1.5)
hold on
plot(TTe_Taylor_quadr,mu_Taylor_o_4/(1.60217e-19),'square k','linewidth',1.5,'MarkerSize', 6)
set(gca,'FontSize',30) 
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$\mu [eV]$','Interpreter','latex','FontSize',30);
legend('Full analytical solution','Our expansion for low temperatures', 'Quadratic approximation','Taylor Expansion from exact solution')
title([' ' prompt2])

figure(476)  %Here we plot the analytical solution of C_e of Aluminium as well as linear approximation, taylor expansion and published data
plot(TTT_el, C_e_alfull,'k','linewidth',1.5)
hold on
plot(250:10:4e4, gamma*(250:10:4e4),'--k','linewidth',1.5)
hold on
plot(TTe_Taylor_quadr,Ce_Taylor_O_3,'square k','linewidth',1.5,'MarkerSize', 6)
%% Uncomment for comparing with data obtained from https://compmat.org/electron-phonon-coupling/
if I_met==1
hold on  %%% Published data from Excel file in the case of aluminion
plot(Virginia_Ce_Data_Al(:,1), Virginia_Ce_Data_Al(:,2),':k','linewidth',1.5)
xlim([0 5.2e4])
end
set(gca,'FontSize',30) 
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$C_e [J m^{-3} K^{-1}]$','Interpreter','latex','FontSize',30);
legend('Full analytical solution','First order approximation $C_e = \gamma  T_e$','Third order taylor expansion $C_e = \gamma  T_e + \beta T_e^3 $','Published Data','Interpreter','latex','FontSize',24)
title(['Electron specific heat of ' prompt2])

%%% Uncomment for comparative plots with 5th order taylor expansion of Ce
figure(1588)  %Here we plot the analytical solution of C_e of Aluminium as well as linear approximation, taylor expansion and published data
plot(TTT_el, C_e_alfull,'k','linewidth',1.5)
% hold on  %%% Published data from Excel file in the case of aluminion
% plot(Virginia_Ce_Data_Al(:,1), Virginia_Ce_Data_Al(:,2),':k','linewidth',1.5)
hold on 
plot(TTT_el, gamma*TTT_el,'k--','linewidth',1.5)
hold on
plot(TTe_Taylor_quadr,Ce_Taylor_O_3,'+ k','linewidth',1.5,'MarkerSize', 7)
% hold on  
% plot(TTe_Taylor_quadr, Ce_Taylor_O_5,'*k','linewidth',1.5)
xlim([0 4.5e4])
set(gca,'FontSize',30) 
xlabel('$T_e [K]$','Interpreter','latex','FontSize',30);
ylabel('$C_e [J m^{-3} K^{-1}]$','Interpreter','latex','FontSize',30);
legend('Exact analytical solution','Linear approximation $C_e = \gamma  T_e$','Third order taylor expansion $C_e = \gamma  T_e + \beta T_e^3 $','Interpreter','latex','FontSize',24)
title( prompt2)


disp(['The maximum of the cubic approximation will be at ' num2str(sqrt(-gamma/(3*beta))) ' K for' ]); disp(prompt2);


if comp_d==1

Virginia_Ce_Data_Ag=readmatrix('Datos_Ce_comparacion.xlsx','Range','O3:P1002');
Virginia_Ce_Data_Cu=readmatrix('Datos_Ce_comparacion.xlsx','Range','Y3:Z1002');
Virginia_Ce_Data_Au=readmatrix('Datos_Ce_comparacion.xlsx','Range','BA3:BB1001');
Virginia_Ce_Data_W=readmatrix('Datos_Ce_comparacion.xlsx','Range','BY3:BZ1002');


Full_anal_Ce_Data_Ag=readmatrix('Datos_Ce_comparacion.xlsx','Range','R3:S3176');
Full_anal_Ce_Data_Cu=readmatrix('Datos_Ce_comparacion.xlsx','Range','AB3:AC3176');
Full_anal_Ce_Data_Au=readmatrix('Datos_Ce_comparacion.xlsx','Range','BD3:BE4981');
Full_anal_Ce_Data_W=readmatrix('Datos_Ce_comparacion.xlsx','Range','CB3:CC4981');

figure(661944)

subplot(2,2,1)
plot(Full_anal_Ce_Data_Ag(:,1),Full_anal_Ce_Data_Ag(:,2),'k','linewidth',1.5)
hold on
plot(Virginia_Ce_Data_Ag(:,1),Virginia_Ce_Data_Ag(:,2),':k','linewidth',1.5)
xlabel('$ T_e $','Interpreter','latex','FontSize',20);
ylabel('$C_e$','Interpreter','latex','FontSize',20);
legend('Analytical solution', 'Published data')
title('Silver')
set(gca,'FontSize',20) 

subplot(2,2,2)
plot(Full_anal_Ce_Data_Cu(:,1),Full_anal_Ce_Data_Cu(:,2),'k','linewidth',1.5)
hold on
plot(Virginia_Ce_Data_Cu(:,1),Virginia_Ce_Data_Cu(:,2),':k','linewidth',1.5)
xlabel('$ T_e $','Interpreter','latex','FontSize',20);
ylabel('$C_e$','Interpreter','latex','FontSize',20);
legend('Analytical solution', 'Published data')
title('Copper')
set(gca,'FontSize',20) 

subplot(2,2,3)
plot(Full_anal_Ce_Data_Au(:,1),Full_anal_Ce_Data_Au(:,2),'k','linewidth',1.5)
hold on
plot(Virginia_Ce_Data_Au(:,1),Virginia_Ce_Data_Au(:,2),':k','linewidth',1.5)
ylim([0 1.1*max(Virginia_Ce_Data_Au(:,2))])
xlabel('$ T_e $','Interpreter','latex','FontSize',20);
ylabel('$C_e$','Interpreter','latex','FontSize',20);
legend('Analytical solution', 'Published data')
title('Gold')
set(gca,'FontSize',20) 

subplot(2,2,4)
plot(Full_anal_Ce_Data_W(:,1),Full_anal_Ce_Data_W(:,2),'k','linewidth',1.5)
hold on
plot(Virginia_Ce_Data_W(:,1),Virginia_Ce_Data_W(:,2),':k','linewidth',1.5)
xlabel('$ T_e $','Interpreter','latex','FontSize',20);
ylabel('$C_e$','Interpreter','latex','FontSize',20);
legend('Analytical solution', 'Published data')
title('Tungsten')
set(gca,'FontSize',20) 



end



function inv_Li = inverting_polylogarithm_3_2(Y)

%%% This function calculates the inverse of the polylogarithm of order 3/2
%%% using Newton-Rapson method 

myfun = @(x,K0) incomplete_polylog(3/2,0,x)-K0;  % parameterized function
K0 = Y;                    % parameter
fun = @(x) myfun(x,K0);    % function of x alone

%%% We want to find the inverse of the polylogarithm, in other words, for
%%% each given Y we want to find the solution of the equation Li_{3/2}
%%% (x)=Y

if Y>-7.889   %%% We use a table of values to find the startinf points of the iterative methods
    
    load aaa.mat aaa
    
    [~,N]= min(abs(aaa-Y));
    
    START_POINT=-100+(N-1)*0.1;
else
    START_POINT=-exp((-3*sqrt(pi)*Y/4)^(2/3));  %%% For Y<< -1, the inverse polilogarithm approaxes this function [1]
end 

inv_Li = fzero(fun,START_POINT);  %%This finds the zero of the function Li_{3/2} (x) - Y=0.

end

function Li_inc = incomplete_polylog(s,b,z)

%%% This function calculates the incomplete polylogarithm of order s

%%% Input ---> s:order of the polylogarithm
%%%            b: first argument, ie, order limit
%%%            z:argument of the polylogarithm

fun = @(x,c) x.^c./(exp(x)/z-1);


Li_inc  = (1/gamma(1+s-1))*  integral(@(x) fun(x,s-1),b,Inf,'RelTol',1e-12);

end



%%% [1] Wood, David C., "The Computation of Polylogarithms." (1992).
%%% [2] Ewa Majchrzak, Jolanta Poteralska, Two-temperature microscale heat transfer model. Part I:
%%% Determination of electrons parameters, Scientic %%%  Research of the Institute of Mathematics and
%%% Computer Science, 2010, Volume 9, Issue 1, pages 99-108.
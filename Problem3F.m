function [] =Problem3F()
% % Use the propagation matrix to calculate the transmission coefficient T for
% % an array of N rectangular potential barriers, each habing width a height
% % V_0 and a barrier to barrier separation b.Calculate and plot how T depends
% % on the energy for an increasing number of barriers. Do you see any
% % oscillations in the transmission? Any regions where the transmision is
% % heavily reduced? What is the underlying physical explanation? Play with
% % the parameters( potential width, heigth, separation), and compare your
% % results wiwth those of theory. When / in what limits is the theoretical
% % model valid? Can you predict where the regions of high transmission (
% % bands) will apear? Also plot the damping -ln (T) and T.
% % we have a width a, distance between the barriers b and a hegth V_0


%% 
% Initialazing terms and vectors

% number of E/V0 values
n=5000;

% number of barriers

nbarrier = 110;

% we'll use the transfer matrix methode
% vector of values of E/V0

EoverV0vect = linspace(0,5,n);

% to store the calculated value for T

Tvect = zeros(1,n);
minuslnTvect = zeros(1,n);
kronigpennyvect = zeros(1,n);
deltakronigpennyvect = zeros(1,n);
abskronigpenneyvect = zeros(1,n);
deltaabskronigpennyvect = zeros(1,n);


% alfa = (V0*(a^2)*2*m)/(h^2); - a the length of the barriers
% We considere alfa dimensioneles and since k_1 is in the order of ~ 10^9
% and a ~ 10^-9 m  we consider alfa^2 to vary betwen 0-3 as for the E/V0
% gama = (V0*(b^2)*2*m)/(h^2); - b the distance between barriers
% Just as for alfa we consider gama dimensionelles since b will be in the order of 10^(-9) m.
% thus we will have Kb and Ka where we'ill have respectivilly b or c multiplied by gama or alfa ^1/2

% Alfa related to the length of the barriers
% Beta related to the distance between barriers

alfa = 0.01;
gama = 200;


%% 
% Here we will calculate the transmition using the  
%   - Matrix Approach
%   - Kronig - Penney model
%   - Delta - Kronnig - Penney model

for o=1:n
    
    
    % b = E/V0 and c = E/V0 - 1 
    
    b = sqrt(EoverV0vect(1,o));
    c = sqrt(EoverV0vect(1,o)-1);
    
    % (((E/V0)*V0*(b^2)*2*m)^(1/2))/(h^2)
    
    k_1gama = b*((gama)^(1/2)); % Related to E/V0 and to the distance between barriers
    
    % (((E/V0)*V0*(a^2)*2*m)^(1/2))/(h^2)
    
    k_1a=b*((alfa)^(1/2)); % Related to E/V0 and the length of the barrier
    k_2a=c*((alfa)^(1/2)); % Related to E/V0 -1 and the length of the barrier
    
    % We know that [A1;B1] = (1/(2*((b*c)^(1/2))))*D12*[A2;B2]
    
    D12 = (1/(2*((b*c)^(1/2)))).*[b+c,b-c;b-c,b+c];
    D21 = (1/(2*((b*c)^(1/2)))).*[b+c,c-b;c-b,b+c];
    
    
    % for x=a we can use the x' = 0 using the relation psi2(x)=psi2'(x-a)
    % in that way we find [A2;B2] = P2 * [A2';B2']
    
    P2 = [exp(-1i*k_2a),0;0,exp(1i*k_2a)];
    
    % Same thing as before but for the third part [A3;B3] = P1 * [A3';B3'], 
    %%we will consider where to pass from psi2'(x')=psi2''(x'-b)
    
    
    P1 = [exp(-1i*k_1gama),0;0,exp(1i*k_1gama)];
    
    % we'll define a variable responsable for representing the
    % passage of the particle by one barrier
    
    trasnfer1 = D12*P2*D21*P1;
    
    % we'll also define a matrix responsible for the passage
    % througth the last boundary
    
    transferbdr = D12*P2;
    
    % the matrix responsible to shift the wave function back to
    % it's original corrdinate system
    
    Pnbarrier = [exp(-1i*(k_1a*(1+(((nbarrier)-3)/2))+k_1gama*...
        (((nbarrier)-3)/2))),0;0,exp(1i*(k_1a*(1+(((nbarrier)-3)/2))+k_1gama*...
        (((nbarrier)-3)/2)))];
    
    
    
    % for  a nbarrier number of barriers we'll obtain a trasfer
    % matrix
    
    transfertot = (trasnfer1^(nbarrier-1))*transferbdr*D21*(Pnbarrier^(-1));
    
    % Trasmission 
    
    T = 1/(abs(transfertot(1))^2);
    
    Tvect(1,o)= T;
    
    % Damping
    
    minuslnTvect(1,o) = -log(T)/log(exp(1));
    
    
    %% 
    % Here I use the Kronig Penney model and the Delta Kronig Penney model
    
    
    % Kronig penney 
    
    kronigpennyvect(1,o) = cos(k_1gama)*cos(k_2a) -(((k_1a^2)+(k_2a^2))/(2*...
        k_1a*k_2a))*sin(k_1gama)*sin(k_2a);
    
    % Dirac Kronig penney 
    
    deltakronigpennyvect(1,o) = cos(k_1gama) + ((gama*alfa)^(1/2)/(2*k_1gama))*...
        sin(k_1gama);
    
    
    % Use if want to plot the absolut value of the Kronigpenneyvect
    % or the absolut value of the Dirac Kronig penney 
    
    abskronigpennyvect(1,o) = abs(kronigpennyvect(1,o));
    deltaabskronigpennyvect(1,o) = abs(deltakronigpennyvect(1,o));
    
end




%% ------------------------------------------------------------------------------------------------- 
% Plotting 

% Open a figure canvas

figure(1);

subplot(2,1,1)

hold on;

% Plot transmission probability as a function of energy
% Simply replace deltakronigpenneyvect by the deltaabskronigpenneyvect or
% the deltaallowedforbiddenvect and it will plot the desidered functions

plot(EoverV0vect,deltaabskronigpennyvect(1,:),'b-',EoverV0vect,  Tvect(1,:),'g-');

% Add labels to axes

xlabel('E/V_0');
ylabel('Translission');

% Add figure title

title(strcat('Transmission coefficient in function of E/V0 for gama = ',...
    num2str(gama),'alpha = ',num2str(alfa)));
legend('delta Kronig-Penney ' ,strcat('T for nbarriers = ',num2str(nbarrier)));

% Scale axes

xlim([0 max(EoverV0vect)]);
ylim([-0.2 1.5]);

% Turn on bounding box and background grid

box on;
grid on;

hold off
subplot(2,1,2)

hold on;

% Plot transmission probability as a function of energy

plot(EoverV0vect, abskronigpennyvect(1,:),'b-',EoverV0vect, Tvect(1,:),'r-');

% Add labels to axes

xlabel('E/V_0');
ylabel('Transmission');

% Add figure title

title(strcat('Transmission coefficient in function of E/V0 for gama = ',...
    num2str(gama),'alpha = ',num2str(alfa)));
legend('Kronig-Penney' ,strcat('T for nbarriers = ',num2str(nbarrier)));

% Scale axes

xlim([0 max(EoverV0vect)]);
ylim([-0.2 1.5]);

% turn on bounding box and background grid

box on;
grid on

figure(2);

hold on;

% Plot transmission probability as a function of energy

plot(EoverV0vect, minuslnTvect(1,:));

% Add labels to axes

xlabel('E/V_0');
ylabel('Transmission');

% Add figure title

title(strcat('Transmission coefficient in function of E/V0 for gama = ',...
    num2str(gama),'and alpha = ',num2str(alfa)));
legend(strcat('Damping') );

% Scale axes

xlim([0 max(EoverV0vect)]);
ylim([0 70]);

% Turn on bounding box and background grid

box on;
grid on

hold off




end


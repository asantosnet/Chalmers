
function [] =Problem2()

%%  Use the propagation matrix approach to calculate the transmission coefficient T of a 
%% rectangular potential well with width a and Height V0.Calculate and plot how T depends on energy.
%% When there are spikes in the trasmission coefficient, we are seeing the bound states



%           k_1                         k_2                      k_1
%                                                               
%                                        2                       
%               1                   3 
%                                                             
%           A1---->                   A2---->                
%           B1 <----                B2 <----                  A3---->
%                                                         B3 <----      
%     ----------------------|V_0                   |-------------------             
%                           |                      |
%                           |                      |
%                           |                      |
%                           |                      |
%                           |                      |
%                           |----------------------|  | energy  0
%                           0                      a
%*******************************************************************



% the constants

%m= 9.1*(10^-31);  %% Kg
%h= 6.6*(10^-16);  %% eV

% number of E/V0 values
n=500;

% number of L values

nl = 500;

% vector of values of E/V0 

EoverV0vect = linspace(0,1,n);

% to store the calculated value for T

Tvect = zeros(8,n);
tvect = zeros(8,n);
detmvect= zeros(8,n);
% b = (E/V0)^1/2
% c = (E/V0-1)^1/2

%L = 1;

%alfa = ((L^2)*2*m)/(h^2);
% we considere alfa dimensioneles and since k_1 is in the order of ~ 10^9
% and L ~ 10^-9 m  we consider alfa^2 to vary betwen 0-3 as for the E/V0 




% under a constant L of 10^-9, we will calculate for three diferent L points
alfa = [625,50,625];
for z=1:3
    for o=1:n

        c = sqrt(EoverV0vect(1,o));
        b = sqrt(EoverV0vect(1,o)-1);


        k_1L=b*((alfa(z))^(1/2));
        k_2L=c*((alfa(z))^(1/2));

        %% we know that [A1;B1] = (1/(2*((b*c)^(1/2))))*D12*[A2;B2]

        D12 = (1/(2*((b*c)^(1/2)))).*[b+c,b-c;b-c,b+c];
        D21 = (1/(2*((b*c)^(1/2)))).*[b+c,c-b;c-b,b+c];


        %% for x=a we can use the x' = 0 using the relation psi2(x)=psi2'(x-a)
        %% in that way we find [A2;B2] = P2 * [A2';B2']

        P2 = [exp(-1i*k_2L),0;0,exp(1i*k_2L)];

        %% Same thing as before but for the third part [A3;B3] = P1 * [A3';B3']

        P1 = [exp(-1i*k_1L),0;0,exp(1i*k_1L)];


        %% we'll find [A1;B1] = D12*P2*D21*P1^-1*[A3;B3]
        %% thus we can use t= D12*P2*D21*P1^-1*[A3;B3] = [t1*A3;t2*B3]
        %% we*re trying to calculate the transmission coefficient 
        %% T= out/in = (abs((A3*exp(1i*k_1*x)))^2)/(abs((A1*exp(1i*k_1*x)))^2) = abs(t1)^ 2

        t = D12*P2*D21/P1;


        T = 1/(abs(t(1))^2);

        Tvect(z,o)= T;
        tvect(z,o)= t(1);
        
        k = (-EoverV0vect(1,o)+1)*((alfa(z))^(1/2));
        
        M = [1,-1,-1,0;
            k/(1i*k_2L),-1,1,0;
            0,((1i*k_2L)/k)*exp(1i*k_2L),(-(1i*k_2L)/k)*exp(-1i*k_2L),exp(-k);
            0,exp(1i*k_2L),exp(-1i*k_2L),-exp(-k)];
        
              
        detmvect(z,o) = det(M);
%          
        
    end

end

%% for a fixed E/V0 value we will plot the transmition, we will now change the distance L, so the value of Alfa from 0 to 3 and for 3 different E/V0 points
% 
% alfa = alfasquare.^(1/2);
% Eoverv0 = [0.5,0.95,1.25,1.5,2];
% for z=1:5
%     for o=1:n
% 
%         c = sqrt(Eoverv0(z));
%         b = sqrt(Eoverv0(z)-1);
% 
% 
%         k_1L=b*((alfa(1,o))^(1/2));
%         k_2L=c*((alfa(1,o))^(1/2));
% 
%         %% we know that [A1;B1] = (1/(2*((b*c)^(1/2))))*D12*[A2;B2]
% 
%         D12 = (1/(2*((b*c)^(1/2)))).*[b+c,b-c;b-c,b+c];
%         D21 = (1/(2*((b*c)^(1/2)))).*[b+c,c-b;c-b,b+c];
% 
% 
%         %% for x=a we can use the x' = 0 using the relation psi2(x)=psi2'(x-a)
%         %% in that way we find [A2;B2] = P2 * [A2';B2']
% 
%         P2 = [exp(-1i*k_2L),0;0,exp(1i*k_2L)];
% 
%         %% Same thing as before but for the third part [A3;B3] = P1 * [A3';B3']
% 
%         P1 = [exp(-1i*k_1L),0;0,exp(1i*k_1L)];
% 
% 
%         %% we'll find [A1;B1] = D12*P2*D21*P1^-1*[A3;B3]
%         %% thus we can use t= D12*P2*D21*P1^-1*[A3;B3] = [t1*A3;t2*B3]
%         %% we*re trying to calculate the transmission coefficient 
%         %% T= out/in = (abs((A3*exp(1i*k_1*x)))^2)/(abs((A1*exp(1i*k_1*x)))^2) = abs(t1)^ 2
% 
%         t = D12*P2*D21/P1;
% 
%         T = 1/(abs(t(1))^2);
% 
%         Tvect(z+3,o)= T;
% 
%     end
% 
% end

% Open a figure canvas
figure(1);

subplot(1,2,1)
% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect, Tvect(1,:));
% Add labels to axes
xlabel('E/V_0');
ylabel('Translission');
% Add figure title
title('Transmission coefficient in function of E/V0');
legend(strcat('alpha = ',num2str(alfa(1))) ,strcat('alpha = ',...
    num2str(alfa(2))) ,strcat('alpha = ',num2str(alfa(3))) );
% Scale axes
xlim([0 max(EoverV0vect)]);
% Turn on bounding box and background grid
box on;
grid on;

figure(2);


% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect, tvect(1,:),'r-',EoverV0vect,zeros(1,n),'b-');
% Add labels to axes
xlabel('E/V_0');
ylabel('Translission');
% Add figure title
title('Transmission coefficient in function of E/V0');
legend(strcat('alpha = ',num2str(alfa(1))) );
% Scale axes
xlim([0 max(EoverV0vect)]);
ylim([-0.00000000004 0.00000000004])
% Turn on bounding box and background grid
box on;
grid on;


% Open a figure canvas
figure(3);

% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect, detmvect(1,:),'r--',EoverV0vect,zeros(1,n),'b-',...
    EoverV0vect, tvect(1,:),'g--');
% Add labels to axes
xlabel('E/V_0');
ylabel('Translission');
% Add figure title
title('det(M) in function of E/V0');
legend(strcat('alpha = ',num2str(alfa(1))));
% Scale axes
xlim([0 max(EoverV0vect)]);
ylim([-0.00000000004 0.00000000004])
% Turn on bounding box and background grid
box on;
grid on;



function [] =Problem1()

% Use the propagation matrix approach to calculate the transmission coefficient T of a 
% rectangular potential step with width L and Height V0.


%           k_1                         k_2                      k_1
%                                                               
%                                        2                       
%               1           |----------------------|  V_0       3 
%                           |                      |             
%           A1---->         |            A2---->   |             
%           B1 <----        |         B2 <----     |             A3---->
%                           |                      |        B3 <----      
%     ----------------------|                      |-------------------  0          | energy 
%                           0                      L
%*******************************************************************



% the constants

%m= 9.1*(10^-31);  %% Kg
%h= 6.6*(10^-16);  %% eV

% number of E/V0 values
n=50000;

% number of L values

nl = 50000;

% vector of values of E/V0 

EoverV0vect = linspace(0,6,n);
EoverV0vect3 = linspace(1,6,n);
EoverV0vect2 = linspace(0,1,n);

% to store the calculated value for T

Tvect = zeros(8,n);
transmissiontest = zeros(8,n);
transmissiontestprime = zeros(8,n);
% b = (E/V0)^1/2
% c = (E/V0-1)^1/2

%L = 1;

%alfa = ((L^2)*2*m)/(h^2);
% we considere alfa dimensioneles and since k_1 is in the order of ~ 10^9
% and L ~ 10^-9 m  we consider alfa^2 to vary betwen 0-3 as for the E/V0 

alfasquare = linspace(0,3000000,nl);


% under a constant L of 10^-9, we will calculate for three diferent L points
alfa = [10,50,100];

for z=1:3
    for o=1:n

        b = sqrt(EoverV0vect(1,o));
        c = sqrt(EoverV0vect(1,o)-1);


        k_1L=b*((alfa(z))^(1/2));
        k_2L=c*((alfa(z))^(1/2));

        % we know that [A1;B1] = (1/(2*((b*c)^(1/2))))*D12*[A2;B2]

        D12 = (1/(2*((b*c)^(1/2)))).*[b+c,b-c;b-c,b+c];
        D21 = (1/(2*((b*c)^(1/2)))).*[b+c,c-b;c-b,b+c];


        % for x=a we can use the x' = 0 using the relation psi2(x)=psi2'(x-a)
        % in that way we find [A2;B2] = P2 * [A2';B2']

        P2 = [exp(-1i*k_2L),0;0,exp(1i*k_2L)];

        % Same thing as before but for the third part [A3;B3] = P1 * [A3';B3']

        P1 = [exp(-1i*k_1L),0;0,exp(1i*k_1L)];


        % we'll find [A1;B1] = D12*P2*D21*P1^-1*[A3;B3]
        %thus we can use t= D12*P2*D21*P1^-1*[A3;B3] = [t1*A3;t2*B3]
        % we*re trying to calculate the transmission coefficient T= out/in = 
        %(abs((A3*exp(1i*k_1*x)))^2)/(abs((A1*exp(1i*k_1*x)))^2) = abs(t1)^ 2

        t = D12*P2*D21/P1;


               Tvect(z,o) = 1/(abs(t(1))^2);


        
        % In order to plot the function that will verify the propagation
        % matrix method

        cprime = sqrt(1-EoverV0vect2(1,o));
        bprime = sqrt(EoverV0vect2(1,o));
        
        cprimeprime = sqrt(EoverV0vect3(1,o)-1);
        bprimeprime = sqrt(EoverV0vect3(1,o));
        
        k_1Lprimeprime=bprimeprime*((alfa(z))^(1/2));
        k_2Lprimeprime=cprimeprime*((alfa(z))^(1/2));
        
        k_1Lprime=bprime*((alfa(z))^(1/2));
        k_2Lprime=cprime*((alfa(z))^(1/2));
        
        transmissiontest(z,o) = abs(((4*1i*bprime*cprime)/((-(cprime^2)+...
            (bprime^2) + 2*1i*cprime*bprime)*exp(k_2Lprime)+((cprime^2)-...
            (bprime^2)+2*1i*bprime*cprime)*exp(-k_2Lprime)))*exp(-1i*k_1Lprime))^2;

        transmissiontestprime(z,o) = (4*(k_1Lprimeprime^2)*...
            (k_2Lprimeprime^2))/(4*(k_1Lprimeprime^2)*(k_2Lprimeprime^2)+...
            (((k_1Lprimeprime^2-k_2Lprimeprime^2)^2)*((sin(k_2Lprimeprime))^2)));

    end

end

% for a fixed E/V0 value we will plot the transmition, we will now change the distance L, so the value of Alfa from 0 to 3 and for 3 different E/V0 points

alfa2 = alfasquare.^(1/2);
Eoverv0 = [0.5,0.95,1.25,1.5,2];
for z=1:5
    for o=1:n

        b = sqrt(Eoverv0(z));
        c = sqrt(Eoverv0(z)-1);


        k_1L=b*((alfa2(1,o))^(1/2));
        k_2L=c*((alfa2(1,o))^(1/2));

        % we know that [A1;B1] = (1/(2*((b*c)^(1/2))))*D12*[A2;B2]

        D12 = (1/(2*((b*c)^(1/2)))).*[b+c,b-c;b-c,b+c];
        D21 = (1/(2*((b*c)^(1/2)))).*[b+c,c-b;c-b,b+c];


        % for x=a we can use the x' = 0 using the relation psi2(x)=psi2'(x-a)
        % in that way we find [A2;B2] = P2 * [A2';B2']

        P2 = [exp(-1i*k_2L),0;0,exp(1i*k_2L)];

        % Same thing as before but for the third part [A3;B3] = P1 * [A3';B3']

        P1 = [exp(-1i*k_1L),0;0,exp(1i*k_1L)];


        % we'll find [A1;B1] = D12*P2*D21*P1^-1*[A3;B3]
        % thus we can use t= D12*P2*D21*P1^-1*[A3;B3] = [t1*A3;t2*B3]
        % we*re trying to calculate the transmission coefficient 
        %T= out/in = (abs((A3*exp(1i*k_1*x)))^2)/(abs((A1*exp(1i*k_1*x)))^2) = abs(t1)^ 2

        t = D12*P2*D21/P1;

        Tvect(z+3,o) = 1/(abs(t(1))^2);


    end

end
% Open a figure canvas
figure(1);


% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect, Tvect(1,:),'k-',EoverV0vect, Tvect(2,:),'g-',EoverV0vect, Tvect(3,:),'b-');
% Add labels to axes
xlabel('E/V_0');
ylabel('Transmission');
% Add figure title
title('Transmission coefficient in function of E/V0');
legend(strcat('alpha = ',num2str(alfa(1))) ,strcat('alpha = ',...
    num2str(alfa(2))) ,strcat('alpha = ',num2str(alfa(3))));
% Scale axes
xlim([0 max(EoverV0vect)]);

% Turn on bounding box and background grid
box on;
grid on;

figure(4);
subplot(2,1,1)
hold on;

plot(alfa2, Tvect(4,:),'r-',alfa2, Tvect(5,:),'g-');

xlabel('alpha');
ylabel('Transmission');

title('Transmission coefficient in function of alpha ( proportional to the length )');
legend('E/V0 = 0.5','E/V0 = 0.95');

xlim([0 max(alfa2)]);

box on;
grid on
subplot(2,1,2)
hold on;

plot(alfa2, Tvect(6,:),'b-',alfa2, Tvect(7,:),'k-',alfa2,Tvect(8,:),'m-');

xlabel('alpha');
ylabel('Transmission');

title('Transmission coefficient in function of alpha ( proportional to the length )');
legend('E/V0=1.25','E/V0=1.5','E/V0 =2');

xlim([0 max(alfa2)]);

box on;
grid on;

% Open a figure canvas
figure(2);

subplot(2,1,1)
% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect2, transmissiontest(1,:),'r-',EoverV0vect3,transmissiontestprime(1,:),'g-');
% Add labels to axes
xlabel('E/V_0');
ylabel('Transmission');
% Add figure title
title('Transmission coefficient in function of E/V0 - Manual Approach');
legend(strcat('alpha = ',num2str(alfa(1))) ,strcat('alpha = ',num2str(alfa(1)))) ;
% Scale axes
xlim([0 max(EoverV0vect)]);

% Turn on bounding box and background grid
box on;
grid on;
subplot(2,1,2)
% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect, Tvect(1,:),'r-');
% Add labels to axes
xlabel('E/V_0');
ylabel('Transmission');
% Add figure title
title('Transmission coefficient in function of E/V0 - Matrix Approach');
legend(strcat('alpha = ',num2str(alfa(1))));
% Scale axes
xlim([0 max(EoverV0vect)]);

% Turn on bounding box and background grid
box on;
grid on;
% Open a figure canvas
figure(3);

subplot(1,1,1)
% Clear the figure canvas
%clf(1);
hold on;
% Plot transmission probability as a function of energy
plot(EoverV0vect, Tvect(1,:),'b-',EoverV0vect2, transmissiontest(1,:),'r-',...
    EoverV0vect3,transmissiontestprime(1,:),'r-');
% Add labels to axes
xlabel('E/V_0');
ylabel('Tranmission');
% Add figure title
title('Transmission coefficient in function of E/V0 - Manual Approach superposed by Matrix Approach');
legend(' Matrix','Manual' ) ;
% Scale axes
xlim([0 max(EoverV0vect)]);

% Turn on bounding box and background grid
box on;
grid on;






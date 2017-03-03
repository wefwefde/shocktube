clear
N = 20;
step = 25000000000;
dt = 0.00000000002;

cp = 0.01;
R = 8.314;
%M = zeros(N);
%B = zeros(N);
% C = zeros(N-1);
x = linspace(0,N,N+1);
dx = 1;

%%%%%%%%%%%%%%%%% Initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(1,N+1);
E = [];
if mod(N,2) == 1
    for i = 1:((N+1)/2)
        E = [0.012 E 0.12];
    end
else
    for i = 1:(N/2)
        E = [0.012 E 0.12];
    end
    E = [E 0.12];
end

rou = [];
if mod(N,2) == 1
    for i = 1:((N+1)/2)
        rou = [0.125 rou 1];
    end
else
    for i = 1:(N/2)
        rou = [0.125 rou 1];
    end
    rou = [rou 1];
end
rou = rou*100;

%%%%%%%%%%%%%%%%%% Calculation vector generation %%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:step
    V1 = zeros(1,N+1);
    V2 = zeros(1,N+1);
    V3 = zeros(1,N+1);
    T = [];
    P = [];
    
    for i = 1:N
        V1(i) = V1(i)+(2/dx)*((2*u(i)*rou(i))/3-(u(i)*rou(i+1))/6-(u(i+1)*rou(i))/6-(u(i+1)*rou(i+1))/3);
        V1(i+1) = V1(i+1)+(2/dx)*((u(i)*rou(i))/3+(u(i)*rou(i+1))/6+(u(i+1)*rou(i))/6-(2*u(i+1)*rou(i+1))/3);
    end
    V1(1) = 0;
    V1(N+1) = 0;
    for i = 1:N
        T(i) = E(i)/(cp*rou(i))-(u(i))^2/(2*cp);
        T(i+1) = E(i+1)/(cp*rou(i+1))-(u(i+1))^2/(2*cp);
        
        V2(i) = V2(i)+(6/(2*dx*rou(i)+dx*rou(i+1)))*((2/3)*T(i)*rou(i)*R-(1/6)*T(i+1)*rou(i)*R-(1/6)*T(i)*rou(i+1)*R...
           -(1/3)*T(i+1)*rou(i+1)*R+(rou(i)*u(i)*u(i))/4+(rou(i)*u(i+1)*u(i))/12+(rou(i+1)*u(i)*u(i))/12+(rou(i+1)*u(i+1)*u(i))/12 ...
           -(rou(i)*u(i)*u(i+1))/4-(rou(i)*u(i+1)*u(i+1))/12-(rou(i+1)*u(i)*u(i+1))/12-(rou(i+1)*u(i+1)*u(i+1))/12);
        V2(i+1) = V2(i+1)+(6/(2*dx*rou(i+1)+dx*rou(i)))*((1/3)*T(i)*rou(i)*R+(1/6)*T(i+1)*rou(i)*R+(1/6)*T(i)*rou(i+1)*R...
           -(2/3)*T(i+1)*rou(i+1)*R+(rou(i)*u(i)*u(i))/12+(rou(i)*u(i+1)*u(i))/12+(rou(i+1)*u(i)*u(i))/12+(rou(i+1)*u(i+1)*u(i))/4 ...
           -(rou(i)*u(i)*u(i+1))/12-(rou(i)*u(i+1)*u(i+1))/12-(rou(i+1)*u(i)*u(i+1))/12-(rou(i+1)*u(i+1)*u(i+1))/4);
    end
    V2(1) = 0;
    V2(N+1) = 0;

    for i = 1:N
        P(i) = rou(i)*R*T(i);
        P(i+1) = rou(i)*R*T(i+1);
        
        V3(i) = V3(i)+(2/(dx))*(2*(u(i)*E(i)+P(i)*u(i))/3 - (u(i)*E(i+1)+u(i)*P(i+1))/6 - (u(i+1)*E(i)+u(i+1)*P(i))/6 - (u(i+1)*E(i+1)+u(i+1)*P(i+1))/3);
        V3(i+1) = V3(i+1)+(2/(dx))*((u(i)*E(i)+u(i)*P(i))/3 + (u(i)*E(i+1)+u(i)*P(i+1))/6 + (u(i+1)*E(i)+u(i+1)*P(i))/6 - 2*(u(i+1)*E(i+1)+u(i+1)*P(i+1))/3);
    end
    V3(1) = 0;
    V3(N+1) = 0;
    
    %%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     r = [0];
%     y1 = [0];
%     y2 = [0];
%     y3 = [0];
%     dr = dx/10;
%     for j = 1:N
%         for l = 1:10
%             r_new(l) = ((j-1)*10 + l)*dr;
%             x_cal = 0.1:0.1:1;
%             y1_new(l) = rou(j)*(1-x_cal(l)) + rou(j+1)*x_cal(l);
%             y2_new(l) = u(j)*(1-x_cal(l)) + u(j+1)*x_cal(l);
%             y3_new(l) = T(j)*(1-x_cal(l)) + T(j+1)*x_cal(l);
%         end
%         r = [r r_new];
%         y1 = [y1 y1_new];
%         y2 = [y2 y2_new];
%         y3 = [y3 y3_new];
%     end
%     count = count + 1;
%     subplot(3,1,1)
%     plot(x,rou)
%     title('Density');
%         
%     subplot(3,1,2)
%     plot(x,u)
%     title('Velocity');
%         
%     subplot(3,1,3)
%     plot(x,T)
%     title('Temperature');
%     count = 0;
    count = count+1;
    if count == 5000000000
       
        subplot(2,2,1)
        plot(x,rou)
        title('Density');
        
        subplot(2,2,2)
        plot(x,u)
        title('Velocity');
        
        subplot(2,2,3)
        plot(x,P)
        title('Pressure');
        
        subplot(2,2,4)
        plot(x,E)
        title('Energy');
        count = 0;
        
    end
%     plot(r,y1,y2,y3)
%     legend('Density','Velocity','Temperature') 
    rou = rou + V1.*dt;
    u = u + V2.*dt;
    E = E + V3.*dt;
end


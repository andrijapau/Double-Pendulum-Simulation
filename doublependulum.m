% Initialize Arrays
theta1 = zeros(1001,1);
theta2 = zeros(1001,1);
theta_p_1 = zeros(1001,1);
theta_p_2 = zeros(1001,1);
theta_pp_2 = zeros(1001,1);
theta_pp_1 =zeros(1001,1);
% Initialize Parameters
g = 9.81;
m1 = 5;
m2 = 10;
l1 = 1;
l2 = 1;
theta1(1) = pi/12;
theta2(1) = 0;
theta_p_1(1) = 0;
theta_p_2(1) = 0;
theta_pp_1(1) = 0;
theta_pp_2(1) = 0;

ti = 0;
tf =500;
N = 10000;
dt = (tf-ti)/N;

duration = tf;;
fps = 60;
movie = false;

nframes=duration*fps;
t = linspace(0,duration,nframes);
figure
subplot(1,2,1)

h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);

range=1.1*(l1+l2); 
axis([-range range -range range]); 
axis square;
set(gca,'nextplot','replacechildren');
for i=1:length(theta1)-1
    num1 = -g*(2*m1+m2)*sin(theta1(i));
    num2 = -m2*g*sin(theta1(i) - 2*theta2(i));
    num3 = -2*sin(theta1(i) - theta2(i))*m2*( ( (theta_p_2(i)*theta_p_2(i)) *l2 )+ ( (theta_p_1(i)*theta_p_1(i)) * l1 * cos(theta1(i)-theta2(i)) ));
    num4 = l1*(2*m1 + m2 - m2*cos(2*(theta1(i)-theta2(i))));
   
    numa = 2*sin(theta1(i)-theta2(i));
    numb = (theta_p_1(i)*theta_p_1(i))*l1*(m1+m2)+g*(m1+m2)*cos(theta1(i));
    numc = (theta_p_2(i)*theta_p_2(i)  )*l2*m2*cos(theta1(i)-theta2(i));
    numd = l2*(2*m1 + m2 - m2*cos(2*theta1(i)-2*theta2(i)));
     % Equations of motion
    theta_pp_1(i) = (num1 + num2 + num3)/num4;
    theta_pp_2(i) = (numa*(numb + numc))/numd;
    
    theta1(i+1) = theta1(i) + theta_p_1(i)*dt+theta_pp_1(i)*0.5*dt*dt;
    theta2(i+1) = theta2(i) + theta_p_2(i)*dt+theta_pp_2(i)*0.5*dt*dt;
    theta_p_1(i+1) = theta_p_1(i) + theta_pp_1(i)*dt;
    theta_p_2(i+1) = theta_p_2(i) + theta_pp_2(i)*dt;
   
    if (ishandle(h)==1)
        Xcoord=[0,l1*sin(theta1(i)),l1*sin(theta1(i))+l2*sin(theta2(i))];
        Ycoord=[0,-l1*cos(theta1(i)),-l1*cos(theta1(i))-l2*cos(theta2(i))];
        set(h,'XData',Xcoord,'YData',Ycoord);
        drawnow limitrate;
        F(i) = getframe;
        if movie==false
            pause(t(i+1)-t(i));
        end
    end
    subplot(1,2,2)
    plot(0,theta1(1));
    hold on;
    plot(i,theta1(i),'--x');
   
   
    
   
end

if movie==true
    movie2avi(F,'doublePendulumAnimation.avi','compression','Cinepak','fps',fps)
end

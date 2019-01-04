theta1 = zeros(1001,1);
theta2 = zeros(1001,1);
theta_p_1 = zeros(1001,1);
theta_p_2 = zeros(1001,1);

g = 9.81;
m1 = 1;
m2 = 1;
l1 = 1;
l2 = 2;
theta1(1) = pi/4;
theta2(1) = 0;
theta_p_1(1) = 0;
theta_p_2(1) = 0;
theta_pp_1(1) = 0;
theta_pp_2(1) = 0;

ti = 0;
tf = 1000;
N = 10000;
dt = (tf-ti)/N;
for t = 1:5000
    
    num1 = -g*((2*m1)+m2)*sin(theta1(t));
    num2 = -m2*g*sin(theta1(t) - (2*theta2(t)));
    num3 = -2*sin(theta1(t) - theta2(t))*m2*( ( (theta_p_2(t)^2) *l2 )+ ( (theta_p_1(t)^2) * l1 * cos(theta1(t)-theta2(t)) ));

    num4 = l1*(2*m1 + m2 - m2*cos(2*(theta1(t) -theta2(t))));

    theta_pp_1 = (num1 + num2 + num3)/num4;

    numa = 2*sin(theta1(t)-theta2(t));
    numb = (theta_p_1(t)^2)*l1*(m1+m2)+g*(m1+m2)*cos(theta1(t));
    numc = (theta_p_2(t)^2)*l2*m2*cos(theta1(t)-theta2(t));
    numd = l2*(2*m1 + m2 - m2*cos(2*theta1(t)-2*theta2(t)));

    theta_pp_2 = (numa*(numb + numc))/numd;
    
    theta1(t+1) = theta1(t) + theta_p_1(t)*dt;
    theta2(t+1) = theta2(t) + theta_p_2(t)*dt;
    theta_p_1(t+1) = theta_p_1(t) + theta_pp_1*dt;
    theta_p_2(t+1) = theta_p_2(t) + theta_pp_2*dt;
    
end
% x1 = l1*(sin(theta1));
% y1 = -1*l1*cos(theta1);
% x2 = x1 + l2*sin(theta2);
% y2 = y1 - l2*cos(theta2);

duration = tf;;
fps = 10;
movie = false;

nframes=duration*fps;
t = linspace(0,duration,nframes);
h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);

range=1.1*(l1+l2); axis([-range range -range range]); axis square;
set(gca,'nextplot','replacechildren');
for i=1:length(theta1)-1
    if (ishandle(h)==1)
        Xcoord=[0,l1*sin(theta1(i)),l1*sin(theta1(i))+l2*sin(theta2(i))];
        Ycoord=[0,-l1*cos(theta1(i)),-l1*cos(theta1(i))-l2*cos(theta2(i))];
        set(h,'XData',Xcoord,'YData',Ycoord);
        
        drawnow;
        F(i) = getframe;
        if movie==false
            pause(t(i+1)-t(i));
        end
    end
end
if movie==true
    movie2avi(F,'doublePendulumAnimation.avi','compression','Cinepak','fps',fps)
end


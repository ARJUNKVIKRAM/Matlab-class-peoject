function RotatingUnbalence(M,m,e,c,k,Omega,time)
%%
%       MECHANICAL VIBRATIONS COURSE PROJECT
%------------------------------------------------------------
%PROJECT DETAILS
%
%project Title: Simulation of Vibration of rotating unbalence system
%Created By, ARJUN K VIKRAM
%
% Format: RotatingUnbalence(Mass,mu,ecc,c,k,omega,time)
%
%Mass is fixed mass mu is roting masses in Kilogram
%ecc is eccentricity in metre
%Also, c is damping coefficient and k is stiffness(in N/m)
%omega is Angular speed (in rev/min)
%time is the total time of analysis in seconds
%
%PROBLEM DETAILS
%
%%
% Checking all the variables
if nargin<7
    M=10;m=1;e=0.3;c=5;k=500;Omega=50;time=4;
end
%Declaration of variables
%
t=0:0.01:time;%time in sec
%
%Calculation of secondary variables
%
w = (2*pi*Omega)/60;%Angular frequency (in rad/sec)
wn = sqrt(k/M);%Natural frequency
r = w/wn; %frequency ratio
del = c/(2*sqrt(k*M));%Damping factor
phi = atan((2*del*r)/(1-r^2));
th = w*t;
%field variables
X = 2*(m*(r^2)*(e/M))/sqrt(((1-r^2)^2)+(2*del*r)^2);
x = X*sin((w*t)-phi);
x_dot = X*w*cos(w*t-phi);
x_ddot = X*(w^2)*(-sin(w*t-phi));
%%
% Initialising the video recording process
%
v = VideoWriter('Record.avi');
v.Quality = 95;
open(v);
xline(0);yline(0); hold on
for i=1:length(t)-1
    s = line([t(1,i),t(1,i+1)],[x(1,i),x(1,i+1)],'Color','b');%Displacement
    vel = line([t(1,i),t(1,i+1)],[x_dot(1,i),x_dot(1,i+1)],'Color','r');%velocity
    a = line([t(1,i),t(1,i+1)],[x_ddot(1,i),x_ddot(1,i+1)],'Color','g');%Acceleration
    c1 = 0-7*e; %for mass_1 co-ordinates (along x axis)
    c2 = c1-4*e;%for mass_2 co-ordinates (along x axis)
    c3 = 0+4*e;%for mass co-ordinates  (along y axis)
    d1 = line([(c1+c2)/2,t(1,i)],[x(1,i+1),x(1,i+1)],'Color','k','LineStyle','--');%follower line
    point1 = plot((c1+c2)/2,x(1,i+1),'k.','Markersize',20,'MarkerFaceColor','k');%ref. point in the sysytem(to indicate the system is vibrating)
    point2 = plot(c1+x(1,i+1),c3+x(1,i+1),'k.','Markersize',5,'MarkerFaceColor','m');%Rotating mass center 1
    point3 = plot(c2+x(1,i+1),c3+x(1,i+1),'k.','Markersize',5,'MarkerFaceColor','m');%Rotating mass center 2
    p1 = plot(c1+e*cos(th(1,i+1))+x(1,i+1),c3+e*sin(th(1,i+1))+x(1,i+1),'k.','Markersize',10,'MarkerFaceColor','b');%Rotating mass 1
    p2 = plot(c2+e*cos(pi-th(1,i+1))+x(1,i+1),c3+e*sin(pi-th(1,i+1))+x(1,i+1),'k.','Markersize',10,'MarkerFaceColor','b');%Rotating mass 2
    l1 = line([c1+x(1,i+1), c1+e*cos(th(1,i+1))+x(1,i+1)],[c3+x(1,i+1), c3+e*sin(th(1,i+1))+x(1,i+1)],'Color','k');
    l2 = line([c2+x(1,i+1), c2+e*cos(pi-th(1,i+1))+x(1,i+1)],[c3+x(1,i+1), c3+e*sin(pi-th(1,i+1))+x(1,i+1)],'Color','k');
    %Plotting system geometry bh=> horizontal lines, bv=>vertical lines
    bh1 = line([c2-(2*e), c1+(2*e)],[c3+2*e+X+x(1,i+1), c3+2*e+X+x(1,i+1)],'Color','k');
    bh2 = line([c2-(5*e), c1+(5*e)],[c3-2*e-X+x(1,i+1), c3-2*e-X+x(1,i+1)],'Color','k');
    bh3 = line([c2-(5*e), c1+(5*e)],[-(c3/2)-2*e-X+x(1,i+1), -(c3/2)-2*e-X+x(1,i+1)],'Color','k');
    bh4 = line([c2-(5*e), c1+(5*e)],[c3-(14*e),c3-(14*e)],'Color','k','Markersize',10);%for fixed base line
    bv1 = line([c2-(2*e), c2-(2*e)],[c3-2*e-X+x(1,i+1), c3+2*e+X+x(1,i+1)],'Color','k');
    bv2 = line([c1+(2*e), c1+(2*e)],[c3-2*e-X+x(1,i+1), c3+2*e+X+x(1,i+1)],'Color','k');
    bv3 = line([c2-(5*e), c2-(5*e)],[-(c3/2)-2*e-X+x(1,i+1), c3-2*e-X+x(1,i+1)],'Color','k');
    bv4 = line([c1+(5*e), c1+(5*e)],[-(c3/2)-2*e-X+x(1,i+1), c3-2*e-X+x(1,i+1)],'Color','k');
    %for plotting spring
    sp1 = -(c3/2)-2*e-X+x(1,i+1); % y value of bottom line of body
    sp2 = c3-(14*e); % y value of base line
    spv = sp2:((sp1-sp2)/100):sp1;
    if x(1,i+1)>=0
        spr1 = plot(c1+(e/2)*sin(8*pi*spv),spv,'k');
        spr2 = plot(c2+(e/2)*sin(8*pi*spv),spv,'k');
    else
        spr1 = plot(c1+(e/2)*sin(10*pi*spv),spv,'k');
        spr2 = plot(c2+(e/2)*sin(10*pi*spv),spv,'k');
    end
    %Plot setting
    %
    legend([s vel a],'Displacement','Velocity','Acceleration','Location','Northeast','Orientation','Horizontal','Fontsize',5);
    pbaspect([2,1,1])
    box on
    grid on
    ylabel({'Displacement [m], Velocity [m/s],'; 'Acceleration [m/s^2]'},'Fontsize',8)
    xlabel('Time (in seconds)','Fontsize',8)
    xlim([-5 5]); ylim([-3.1 3]);
    set(gca,'XTick',0:0.5:5);
    set(gca,'YTick',-3:0.5:3);
    title('Vibration of Rotating Mass Unbalence','Fontsize',10)
    pause(0.01)
    frame = getframe(gcf);
    writeVideo(v,frame);
    % Deleting entities of the (i-1)th frame to make the simulation look clean
    if i==length(t)-1
        break;
    else
        delete(d1);delete(point1);delete(point2);delete(point3);delete(l1);delete(l2);delete(bh1);delete(bh2);delete(bh3);delete(bh4);
        delete(bv1);delete(bv2);delete(bv3);delete(bv4);delete(point2);delete(point3);delete(p1);delete(p2);delete(spr1);delete(spr2);
    end
end
close(v);
grid off
end

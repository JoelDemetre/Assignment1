%% Assignment 1
% Due: Feburary 4th 2018 \n
% Completed By: Joel Demetre (100943543)


%% Question 1A) Velocity
% Calculation of Velocity
%
% $$ vi = \frac{\KT}{mn}
%
T = 300;
Kbolt = 1.38064852E-23;
mo = 9.11e-31;
mn = mo * 0.26;
vx = (Kbolt*T/mn)^.5;
vy = (Kbolt*T/mn)^.5;
vth = (vx^2 + vy^2)^.5;

%% Question 1B) Mean Free Path
tmn = 0.2e-12;
MeanFreePath = tmn*vth;


%% Question 1C) Program


xlimits = [0, 2e-9];
ylimits = [0, 1e-9];
PlotHowMany = 10;
Timestep = 5e-17;
endtime = Timestep*10;
NumParticles = 1000;




%%SETUP INITIAL PARAMETERS
mycolors = hsv(PlotHowMany);
xprev = zeros(1, NumParticles);
yprev = zeros(1, NumParticles);
x = zeros(3,NumParticles);
y = zeros(2,NumParticles);
temp = zeros(2, NumParticles);

%Start the random distribution in x position
x(1,:) =  rand(1,NumParticles);
y(1,:) =  rand(1,NumParticles);
x(1,:) =  xlimits(1) + x(1,:).*(xlimits(2) - xlimits(1));
y(1,:) =  ylimits(1) + y(1,:).*(ylimits(2) - ylimits(1));

%Assign the random velocity and random angle
x(3,:) = rand(1, NumParticles)*2*pi;
x(2,:) = vth.*cos(x(3,:));
y(2,:) = vth.*sin(x(3,:));

%%LOOP
figure(1);
for i = 0:Timestep:endtime
   % pause(.1);
    xprev(1,:) = x(1,:);
    yprev(1,:) = y(1,:);
    for kt = 1:NumParticles
       if x(1,kt) +  x(2,kt)*Timestep < xlimits(1)
          x(1,kt) = xlimits(2);
          xprev(1,kt) = xlimits(2);
       elseif x(1,kt) +  x(2,kt)*Timestep > xlimits(2)
           x(1,kt) = xlimits(1);
           xprev(1,kt) = xlimits(1);
       end
       if y(1,kt) +  y(2,kt)*Timestep < ylimits(1) || y(1,kt) +  y(2,kt)*Timestep > ylimits(2)
           y(2,kt) = -y(2,kt);
       end
        x(1,kt) = x(1,kt) + x(2,kt).*Timestep;
        y(1,kt) = y(1,kt) + y(2,kt).*Timestep;
    if kt <= PlotHowMany
          plot([xprev(1,kt), x(1,kt)], [yprev(1,kt), y(1,kt)],'color', mycolors(kt,:) );
          hold on;
    end
   end
    hold on;
    VelSquared = mean((x(2,:).^2 + y(2,:).^2));
    CalcTemp = VelSquared*mn/2/Kbolt;
    title(['Average Temperature: ' num2str(CalcTemp)]);
    xlim([xlimits(1), xlimits(2)]);
    ylim([ylimits(1), ylimits(2)]);
end




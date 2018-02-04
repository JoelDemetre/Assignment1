clear all; %#ok<CLALL>
close all;
xlimits = [0, 2e-9];
ylimits = [0, 1e-9];
PlotHowMany = 10;
Timestep = 9e-17;
endtime = Timestep*1000;
NumParticles = 1000;



T = 300;
Kbolt = 1.38064852E-23;
mo = 9.11e-31;
mn = mo * 0.26;
vx = (Kbolt*T/mn)^.5;
vy = (Kbolt*T/mn)^.5;
vth = (vx^2 + vy^2)^.5;

tmn = 0.2e-12;
MeanFreePath = tmn*vth;

Pscat = 1-exp(-Timestep/tmn);
mycolors = hsv(PlotHowMany);
xprev = zeros(1, NumParticles);
yprev = zeros(1, NumParticles);
x = zeros(3,NumParticles);
y = zeros(2,NumParticles);
temp = zeros(2, NumParticles);
scatTime = zeros(1,NumParticles);

%Start the random distribution in x position
x(1,:) =  rand(1,NumParticles);
y(1,:) =  rand(1,NumParticles);
x(1,:) =  xlimits(1) + x(1,:).*(xlimits(2) - xlimits(1));
y(1,:) =  ylimits(1) + y(1,:).*(ylimits(2) - ylimits(1));

%Assign the random velocity and random angle
temp(1,:) = normrnd(vth, sqrt(mn/(Kbolt*T)), 1, NumParticles);
x(3,:) = rand(1, NumParticles)*2*pi;
x(2,:) = temp(1,:).*cos(x(3,:));
y(2,:) = temp(1,:).*sin(x(3,:));



figure(1);
%Update the Position
counter = 1;
for i = 0:Timestep:endtime
    pause(.1);
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
    %Scattering Check
    if Pscat>rand()
        temp = normrnd(vth, sqrt(mn/(Kbolt*T)));
        x(3,kt) = rand*2*pi;
        x(2,kt) = temp.*cos(x(3,kt));
        y(2,kt) = temp.*sin(x(3,kt));
        scatTime(1,kt) = 0;
    end
    end
    hold on;
    AvgScat = mean(scatTime(1,:));
    scatTime(1,:) = scatTime(1,:) + Timestep;
    VelSquared = mean((x(2,:).^2 + y(2,:).^2));
    CalcTemp = VelSquared*mn/2/Kbolt;
    title(['Average Temperature: ' num2str(CalcTemp) ' Average Scatter Time: ' num2str(AvgScat)]);
    xlim([xlimits(1), xlimits(2)]);
    ylim([ylimits(1), ylimits(2)]);
end

figure(2);
hist(sqrt(x(2,:).^2 + y(2,:).^2),100);



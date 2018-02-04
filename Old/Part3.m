%Initialize Parameters
%Area Limits
xlimits = [0, 2e-9];
ylimits = [0, 1e-9];
PlotHowMany = 10;
Kbolt = 1.38064852E-23;
Timestep = 3e-14;
endtime = 1e-10;
NumParticles = 2000;
Temp = 270;
mass = 9.11e-31;
SetVelAverage =900;
toumn = 0.2e-12;
DiffusionBarrierProbability = 0.2;



%Assign the Box Height and Width as Percentages of the Limits
xboxLim = [.4*(xlimits(2)-xlimits(1)), .6*(xlimits(2)-xlimits(1))];
yboxLim1 = [ylimits(1), .4*(ylimits(2)-ylimits(1))];
yboxLim2 = [.6*(ylimits(2)-ylimits(1)), ylimits(2)];
xbox = xboxLim([1 1 2 2 1]);
ybox1 = yboxLim1([1 2 2 1 1]);
ybox2 = yboxLim2([1 2 2 1 1]);





Pscat = 1-exp(-Timestep/toumn);
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

%Get the Particle Indexes that are out of bounds
IDX = uint32(1:NumParticles); 
ind = IDX((y(1,:) >= yboxLim2(1) | y(1,:) <= yboxLim1(2))& x(1,:) >= xboxLim(1) & x(1,:) <= xboxLim(2));

%Reassign Positions
counter = 1;
%size(ind, 2)
while counter <= size(ind,2)
    x(1,ind(counter)) = xlimits(1) + rand*(xlimits(2) - xlimits(1));
    y(1,ind(counter)) = ylimits(1) + rand*(ylimits(2) - ylimits(1));
    if ~((y(1,ind(counter)) >= yboxLim2(1) || y(1,ind(counter)) <= yboxLim1(2))&& x(1,ind(counter)) >= xboxLim(1) && x(1,ind(counter)) <= xboxLim(2))
       counter = counter + 1; 
    end
end

%Assign the random velocity and random angle
temp(1,:) = normrnd(SetVelAverage, sqrt(mass/(Kbolt*Temp)), 1, NumParticles);
x(3,:) = rand(1, NumParticles)*2*pi;
x(2,:) = temp(1,:).*cos(x(3,:));
y(2,:) = temp(1,:).*sin(x(3,:));



figure(1);
xlim([xlimits(1), xlimits(2)]);
ylim([ylimits(1), ylimits(2)]);
hold on;
rectangle('Position', [xboxLim(1), yboxLim1(1), xboxLim(2)-xboxLim(1), yboxLim1(2)]);
rectangle('Position', [xboxLim(1), yboxLim2(1), xboxLim(2)-xboxLim(1), yboxLim2(2) - yboxLim2(1)]);

counter = 1;
for i = 0:Timestep:endtime
   % pause(.1);
    xprev(1,:) = x(1,:);
    yprev(1,:) = y(1,:);
    for kt = 1:NumParticles
        
       %PARTICLE'S HITTING TOP AND BOTTOM
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
        
       %PARTICLE'S HITTING BOX
         if (y(1,kt) +  y(2,kt)*Timestep >= yboxLim2(1) && x(1,kt) +  x(2,kt)*Timestep >= xboxLim(1) && x(1,kt) +  x(2,kt)*Timestep <= xboxLim(2))
           [xinter1, xinter2, yinter1, yinter2]  = BoxIntercept( x(1,kt), y(1,kt), x(1,kt) +  x(2,kt)*Timestep, y(1,kt) +  y(2,kt)*Timestep, xboxLim(1), xboxLim(2),yboxLim2(1), ylimits(2));
           if rand > DiffusionBarrierProbability
            if xinter2 || xinter1
                 x(2,kt) = - x(2,kt);
             end
             if yinter2 || yinter1
                 y(2,kt) = - y(2,kt);
             end
           else
               temp = normrnd(SetVelAverage, sqrt(mass/(Kbolt*Temp)));
               if xinter1
                   x(3,kt) = pi/2 + rand*pi;
               elseif xinter2
                   x(3,kt) = pi/2 - rand*pi;
               elseif yinter2
                   x(3,kt) = rand*pi;
               elseif yinter1
                   x(3,kt) = -rand*pi;
               end
               x(2,kt) = temp.*cos(x(3,kt));
               y(2,kt) = temp.*sin(x(3,kt));
           end
         elseif y(1,kt) +  y(2,kt)*Timestep <= yboxLim1(2) && x(1,kt) +  x(2,kt)*Timestep >= xboxLim(1) && x(1,kt) +  x(2,kt)*Timestep <= xboxLim(2)
           [xinter1, xinter2, yinter1, yinter2]  = BoxIntercept( x(1,kt), y(1,kt), x(1,kt) +  x(2,kt)*Timestep, y(1,kt) +  y(2,kt)*Timestep, xboxLim(1), xboxLim(2), ylimits(1), yboxLim1(2));
           if rand > DiffusionBarrierProbability
           if xinter1 || xinter2
                x(2,kt) = - x(2,kt);
            end
            if yinter1 || yinter2
                y(2,kt) = - y(2,kt);
            end
           else
               temp = normrnd(SetVelAverage, sqrt(mass/(Kbolt*Temp)));
               if xinter1
                   x(3,kt) = pi/2 + rand*pi;
               elseif xinter2
                   x(3,kt) = pi/2 - rand*pi;
               elseif yinter2
                   x(3,kt) = rand*pi;
               elseif yinter1
                   x(3,kt) = -rand*pi;
               end
               x(2,kt) = temp.*cos(x(3,kt));
               y(2,kt) = temp.*sin(x(3,kt));
           end
         end
         
    %UPDATE POSITION   
    x(1,kt) = x(1,kt) + x(2,kt).*Timestep;
    y(1,kt) = y(1,kt) + y(2,kt).*Timestep;
    %PLOT
    if kt <= PlotHowMany
        plot([xprev(1,kt), x(1,kt)], [yprev(1,kt), y(1,kt)],'color', mycolors(kt,:) );
        hold on;
    end
    %Scattering Check
    if Pscat>rand()
        temp = normrnd(SetVelAverage, sqrt(mass/(Kbolt*Temp)));
        x(3,kt) = rand*2*pi;
        x(2,kt) = temp.*cos(x(3,kt));
        y(2,kt) = temp.*sin(x(3,kt));
        scatTime(1,kt) = 0;
    end
    end
    hold on;
    AvgScat = mean(scatTime(1,:));
    scatTime(1,:) = scatTime(1,:) + Timestep;
    averagevel = mean((x(2,:).^2 + y(2,:).^2).^.5);
    title(['Average Temperature: ' num2str(averagevel) ' Average Scatter Time: ' num2str(AvgScat)]);
    figure(2);
    hist(sqrt(x(2,:).^2 + y(2,:).^2),100);
    figure(1);
end










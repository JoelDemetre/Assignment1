%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlimits = [0, 2e-9];
ylimits = [0, 1e-9];
PlotHowMany = 1;
Timestep = 0.1;
endtime = 10000;
NumParticles = 1000;
%Temperature in K
T = 0;
velocity = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kbolt = 1.38064852E-23;
m = 9.11e-31;

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
temp(1,:) = velocity;
x(3,:) = rand(1, NumParticles)*2*pi;
x(2,:) = temp(1,:).*cos(x(3,:));
y(2,:) = temp(1,:).*sin(x(3,:));

transistionArray = [0,0,0,0;0,0,0,0;];

%Update the Position
counter = 1;
for i = 0:Timestep:endtime
    pause(.1);
    xprev(1,:) = x(1,:);
    yprev(1,:) = y(1,:);
    for kt = 1:NumParticles
       transistion = false;
       if y(1,kt) +  y(2,kt)*Timestep < ylimits(1) || y(1,kt) +  y(2,kt)*Timestep > ylimits(2)
           y(2,kt) = -y(2,kt);
       end
       if x(1,kt) +  x(2,kt)*Timestep < xlimits(1)
           x(1,kt) = x(1,kt) +  x(2,kt)*Timestep;
           x(1,kt) = xlimits(2)+ x(1,kt);
           y(1,kt) = y(1,kt) + y(2,kt).*Timestep;
           transistion = true;
           transistionArray(1,1) = xprev(1,kt);
           transistionArray(1,2) = xlimits(1);
           transistionArray(1,3) = xlimits(2);
           transistionArray(1,4) = x(1,kt);
           transistionArray(2,1) = yprev(1,kt);
           temp1 = transistionArray(1,1)-transistionArray(1,2);
           temp2 = transistionArray(1,3)-transistionArray(1,4);
           transistionArray(2,2) = yprev(1,kt) + (y(1,kt) -yprev(1,kt))*abs(temp1/(temp1 + temp2));
           transistionArray(2,3) = yprev(1,kt) + (y(1,kt) -yprev(1,kt))*abs(temp1/(temp1 + temp2));
           transistionArray(2,4) = y(1,kt);
       elseif x(1,kt) +  x(2,kt)*Timestep > xlimits(2)
           x(1,kt) = x(1,kt) +  x(2,kt)*Timestep;
           x(1,kt) = x(1,kt) - xlimits(2);
           y(1,kt) = y(1,kt) + y(2,kt).*Timestep;
           transistion = true;
           transistionArray(1,1) = xprev(1,kt);
           transistionArray(1,2) = xlimits(2);
           transistionArray(1,3) = xlimits(1);
           transistionArray(1,4) = x(1,kt);
           transistionArray(2,1) = yprev(1,kt);
           temp1 = transistionArray(1,1)-transistionArray(1,2);
           temp2 = transistionArray(1,3)-transistionArray(1,4);
           transistionArray(2,2) = yprev(1,kt) + (y(1,kt) -yprev(1,kt))*abs(temp1/(temp1 + temp2));
           transistionArray(2,3) = yprev(1,kt) + (y(1,kt) -yprev(1,kt))*abs(temp1/(temp1 + temp2));
           transistionArray(2,4) = y(1,kt);
       end
    if ~transistion
        x(1,kt) = x(1,kt) + x(2,kt).*Timestep;
        y(1,kt) = y(1,kt) + y(2,kt).*Timestep;
    end
    if kt <= PlotHowMany
        if transistion 
             plot([transistionArray(1,1),transistionArray(1,2)], [transistionArray(2,1),transistionArray(2,2)],'color', mycolors(kt,:));
             plot([transistionArray(1,3),transistionArray(1,4)], [transistionArray(2,3),transistionArray(2,4)],'color', mycolors(kt,:));
             hold on;
        else
          plot([xprev(1,kt), x(1,kt)], [yprev(1,kt), y(1,kt)],'color', mycolors(kt,:) );
          hold on;
        end
    end
   end
    hold on;
    averagevel = mean((x(2,:).^2 + y(2,:).^2).^.5);
    CalcTemp = 0;
    title(['Average Temperature: ' num2str(averagevel)]);
    xlim([0, 1000]);
    ylim([0, 1000]);
end





% Solve the 1-D wave equation subject to
% fixed boundary condition on right
% moving boundary condition on left
% -------------------------------------------

clear all
close all
clc

nFrames = 400;                % number of frames in the movie (time!)
tic
c      = 2;                   % set physical parameters
deltaT = 0.02;
deltaX = 0.1;
deltaY = deltaX;
lamb   = (c*deltaT/deltaX)^2;

xRight = 3;
yRight = 3;
x = linspace(0, xRight, (xRight/deltaX +1) ); % calculate x node locations
y = linspace(0, yRight, (yRight/deltaY +1) ); % calculate y node locations
n = length(x);
m = length(y);

%Initialize the matrices that handle the values of the height of the wave
%function u(x,y) 

uInit = zeros(n,m);         
uCur = uInit;
uFut = uInit;
%Create a mesh grid of linearly spaced x and y values that correspond with
%the values of the height of the wave function

[x,y] = meshgrid(0:(5/(xRight/deltaX)):5, 0:(5/(yRight/deltaY)):5);

%Plot the flat meshplot so that it does not have to be initialized over and
%over again in the for loop that handles the caculation for timesteps k>=2

meshc(x,y,uFut)
drawnow

%Initialize the 'wiggle' in our uCur matrix so that we have something to
%work off of when we iterate our finite difference method

uCur(:,1) = leftBoundary(1, deltaT); 

% calculate the first time-step, it is only necessary to go from indices
% 2:4 because the intiial wiggle only affects the first few calculated
% values of u(x,y), this saves time and speeds up the program.

for i = 2 : 4            
   for j = 2 : 4

    uCur(i,j) =  lamb*(uInit(i-1,j)+uInit(i+1,j)+uInit(i,j-1)+uInit(i, j+1)) +...
                 (2-4*lamb)*uInit(i,j) / 2;
    end
end

% and the rest of the time steps...
for k = 2 : nFrames
   if k < 100
    uFut(:,1) = leftBoundary(k, deltaT);
   else
       display(uCur(:,2));
       uFut(:,1) = uCur(:,2);
   end
    for i = 2 : n-1
            a = i*deltaX;
            for j = 2 : m-1
            [h,hx,hy] = depth2D(a, j*deltaY);
          
            uFut(i,j) = lamb*uCur(i+1,j)*((hx*deltaX/2)+h)+lamb*uCur(i,j+1)*((hy*deltaY/2)+h)+...
                      lamb*uCur(i-1,j)*(h-(hx*deltaX/2)) + lamb*uCur(i,j-1)*(h-((hy*deltaY/2)))+...
                      uCur(i,j)*(2-4*lamb*h) - uInit(i,j);
            end
    end
     
    
        meshc(x, y, uFut);
        set(gca, 'zlim', [-.5 .5])
        xlabel('X label'),
        ylabel('Y label'),
        zlabel('Z label'),
        grid on,
        title('Mesh plot and contour of peaks function'),
        hidden on             
        drawnow
        
      uInit = uCur;            % update u values
      uCur = uFut;
end
toc
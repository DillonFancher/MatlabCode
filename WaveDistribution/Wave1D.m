
% Solve the 1-D wave equation subject to
% fixed boundary condition on right
% moving boundary condition on left
% -------------------------------------------

clear all
close all
clc

nFrames = 400;                % number of frames in the movie (time!)

c      = 1;                   % set physical parameters
deltaT = 0.02;
deltaX = 0.1;
deltaY = deltaX;
lamb   = (c*deltaT/deltaX)^2;

xRight = 5;
yRight = 5;
x = linspace(0, xRight, (xRight/deltaX +1) ); 
y = linspace(0, yRight, (yRight/deltaY +1) );% calculate x node locations
n = length(x);
m = length(y);

uInit = zeros(n,m);         
uCur = uInit;
uFut = uCur;

[x,y] = meshgrid(0:(5/(xRight/deltaX)):5, 0:(5/(yRight/deltaY)):5);
drawnow

for j = 1:m
uCur(j,1) = leftBoundary(1, deltaT); 
end

% % set spatial boundary values
for j = 2 : m -1
for i = 2 : n -1                 % calculate the first time-step
    uCur(i,j) =  lamb*(uInit(i-1,j)+uInit(i+1,j)+uInit(i,j-1)+uInit(i, j+1)) +...
                 (2-4*lamb)*uInit(i,j) / 2;
end
end
display(uCur)
meshc(x,y,uFut)
xlabel('X label'),
   ylabel('Y label'),
   zlabel('Z label'),
   grid on,
   title('Mesh plot and contour of peaks function'),
   hidden on             
drawnow



uFut(1,1) = 0;
uFut(n,m) = 0;
                              % and the rest of the time steps...
for k = 2 : nFrames
   for j = 1:m
    uFut(j,1) = leftBoundary(k, deltaT);
    
   end
    for j = 2 : m-1
      for i = 2 : n-1
          uFut(i,j) = ( lamb*(uCur(i-1, j) + uCur(i+1,j) + uCur(i,j-1) + uCur(i,j+1)) +...
                      (2-4*lamb)*uCur(i,j) - uInit(i,j));
      end
  
    end 
   
    
        meshc(x, y, uFut);
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


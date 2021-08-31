function read_file %#codegen

clc
close all

data = readmatrix('data.txt');

%%%% define parameters for the pendulum %%%%%%
l1 = 1; l2 = 1; l3 = 0.5;

%Target x and y
x_ref = 1; y_ref = 1;


theta = [data(1,1) data(1,2) data(1,3); data(3,1) data(3,2) data(3,3)];


%%%%%% plot %%%%%%%%%%
for i=1:2
    %%%%%%%% prepping to get homogenous transformations %%%%
    c1 = cos(theta(i,1)); s1 = sin(theta(i,1));
    c2 = cos(theta(i,2)); s2 = sin(theta(i,2));
    c3 = cos(theta(i,3)); s3 = sin(theta(i,3));
    O01 = [0; 0]; O12 = [l1; 0]; O23 = [l2; 0];
    R01 = [c1 -s1; s1 c1];
    R12 = [c2 -s2; s2 c2];
    R23 = [c3 -s3; s3 c3];
    H01 = [R01, O01; 0, 0, 1];
    H12 = [R12, O12; 0, 0, 1];
    H23 = [R23, O23; 0, 0, 1];
    
    %%%%%%%% origin  in world frame  %%%%%%
    o = [0 0];
    
    %%%%% end of link1 in world frame %%%%
    P1 = [l1; 0; 1];
    P0 = H01*P1;
    p = P0(1:2); %same as p0
    
    %%%% end of link 2 in world frame  %%%%%%%
    Q2 = [l2; 0; 1];
    Q0 = H01*H12*Q2;
    q = Q0(1:2); %same as q0
    
    %%%% end of link 3 in world frame  %%%%%%%
    R3 = [l3; 0; 1];
    R0 = H01*H12*H23*R3;
    r = R0(1:2); %same as r0
    
    %%%%%% draw the curve on paper %%%
    plot(r(1), r(2), 'ko','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
    
    %Draw line from origin to end of link 1
    h1 = line([o(1) p(1)],[o(2) p(2)],'LineWidth',5,'Color','red');
    
    %Draw line from end of link 1 to end of link 2
    h2 = line([p(1) q(1)],[p(2) q(2)],'LineWidth',5,'Color','blue');
    
    %Draw line from end of link 2 to end of link 3
    h3 = line([q(1) r(1)],[q(2) r(2)],'LineWidth',5,'Color','green');
    
    text(data(2,1),data(2,2),'Start');
    text(x_ref,y_ref,'Goal');
    xlabel('x'); ylabel('y');
    grid on; %if you want the grid to show up.
    axis('equal'); %make the axis equal, to avoid scaling effect
    
    % These set the x and y limits for the axis (will need adjustment)
    xlim([-3 3]);
    ylim([-3 3]);
end  

end

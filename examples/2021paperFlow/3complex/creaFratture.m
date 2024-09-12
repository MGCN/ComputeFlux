clear all
close all

pt1=[0.0500 
0.0500
0.1500
0.1500
0.6500
0.7000
0.6000
0.3500
0.7500
0.1500];

pt2=[0.4160 0.2200
0.2750 0.2500
0.6300 0.4500
0.9167 0.4000
0.8333 0.849723
0.2350 0.849723
0.3800 0.8500
0.9714 0.8000
0.9574 0.9500
0.8363 0.4000];

pt3=[0.0624 0.1350 0.0900 0.5000 0.167625 0.167625 0.2675 0.7143 0.8155 0.9727]';

p=[pt1 pt2 pt3];

pa=p(:,1:2);
pb=p(:,3:4);

fx_string  = ['fx_string = '''];
fy_string  = ['fy_string = '''];
fd1_string = ['fd1_string = '''];
fd2_string = ['fd2_string = '''];
fa1_string = ['fa1_string = '''];

for i=1:10
    a=pa(i,:)';
    b=pb(i,:)';
    
    x=[a(1) b(1)];
    y=[a(2) b(2)];
    
    c=(a+b)/2;
    l=norm(b-a);
    
    aa=(a-c)/(l/2);
   
    if (aa(2)>=0.0)
        angolo=acos(aa(1));
    else
        angolo=2*pi-acos(aa(1));
    end
    
    angolo=angolo/pi*180;
    
    plot(x,y)
    hold on
    plot(a(1),a(2),'*')
    
    fx_string = [fx_string,' ',num2str(c(1))];
    fy_string = [fy_string,' ',num2str(c(2))];
    fd1_string = [fd1_string,' ',num2str(l)];
    fd2_string = [fd2_string,' ',num2str(1e-4)];
    fa1_string = [fa1_string,' ',num2str(angolo)];
    %angolo
    
end

    fx_string = [fx_string,''''];
    fy_string = [fy_string,''''];
    fd1_string = [fd1_string,''''];
    fd2_string = [fd2_string,''''];
    fa1_string = [fa1_string,''''];
    
    disp(fx_string)
    disp(fy_string)
    disp(fd1_string)
    disp(fd2_string)
    disp(fa1_string)
    
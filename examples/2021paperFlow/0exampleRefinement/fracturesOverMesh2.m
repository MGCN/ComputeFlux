clear all
close all

l=0.8/2;
a=0.05/2;

frac1=[-l  l l -l;
       -a -a a  a];
   
   theta=30/180*pi;
   rot=[cos(theta) -sin(theta); sin(theta) cos(theta)];
   
   frac2=rot*frac1;

   
   zoom=[0.2  0.4 0.4 0.2;
       0.1 0.1 0.3 0.3];
   
   
J=349;
I=749;

a=dir('*.jpeg');

for i=1:length(a)
    
    h=figure;
    inputName=a(i).name;
    outputName=[inputName(1:end-5),'.pdf'];
    outputName2=['zoom',outputName(5:end)];
    
    img = imread(inputName);
    img = img(1+J:end-J,1+I:end-I,:);
    
    imagesc([-0.5 0.5], [-0.5 0.5], flipdim(img,1) ); %flipud(img));

    
    hold on

    plot([frac1(1,:) frac1(1,1)],[frac1(2,:) frac1(2,1)],'LineWidth',2,'Color','r')
    plot([frac2(1,:) frac2(1,1)],[frac2(2,:) frac2(2,1)],'LineWidth',2,'Color','r')
    plot([zoom(1,:) zoom(1,1)],[zoom(2,:) zoom(2,1)],'LineWidth',1.5,'Color','b')
    
    set(gca,'ydir','normal');
    yticks([])
    xticks([])
    %axis equal
    set(gca,'FontSize', 30);  
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
    print(h,outputName,'-dpdf','-r0')

    xlim([0.2 0.4])
    ylim([0.1 0.3])
    print(h,outputName2,'-dpdf','-r0')
    
end
close all

!mv *pdf /Users/favinom/Desktop/wccm2021/figuresW/0intro

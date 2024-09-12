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

pa=[pa; 0, 0.5];
pb=[pb; 1, 0.9];


preName='realisticMesh';
postName='.jpg';

J=175;
I=375;

a=dir('*.png');

for i=1:length(a)
    
    h=figure;
    inputName=a(i).name;
    outputName=[inputName(1:end-4),'.pdf'];
    
    img = imread(inputName);
    img = img(1+J:end-J,1+I:end-I,:);
    
    imagesc([0 1], [0 1], flipdim(img,1) ); %flipud(img));

    
    hold on
    
    jstart=11;
    if(i==1)
        jstart=1;
    end
    
    for j=jstart:size(pa,1)
        
        color='r';
        color2='.r';
        lw=2;
        if (j==4 || j==5)
            color='b';
            color2='.b';
        end
        if (j==11)
            color='--k';
            color2='.k';
            lw=1;
        end
        
        p1=pa(j,:);
        p2=pb(j,:);
        
        x=[p1(1) p2(1)];
        y=[p1(2) p2(2)];
        plot(x,y,color,'LineWidth',lw)
        
        plot(p1(1),p1(2),color2,'MarkerSize',20)
        plot(p2(1),p2(2),color2,'MarkerSize',20)
    end
    %axis equal
    set(gca,'ydir','normal');
    yticks([0 0.5 1])
    %axis equal
    set(gca,'FontSize', 30);  
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
    print(h,outputName,'-dpdf','-r0')

end
close all

%l!mv *pdf /Users/favinom/Desktop/wccm2021/figuresW/3complex

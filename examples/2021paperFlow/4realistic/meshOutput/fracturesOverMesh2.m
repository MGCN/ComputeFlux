clear all
close all

p=csvread('network.csv',1);

p=p(:,2:end);

pa=p(:,1:2);
pb=p(:,3:4);

pa=[pa; 0, 500];
pa=[pa; 625, 0];

pb=[pb; 700, 500];
pb=[pb; 625, 600];

I=343;
J=209;


a=dir('*.png');

for i=1:length(a)
    
    h=figure;
    inputName=a(i).name;
    outputName=[inputName(1:end-4),'.pdf'];
    
    img = imread(inputName);
    
    img = img(1+J:end-J,1+I:end-I,:);
    
    imagesc([0 700], [0 600], flipdim(img,1) ); %flipud(img));

    hold on
    
    jstart=64;
    if(i==1)
        jstart=1;
    end
    
    for j=jstart:size(pa,1)
        
        color='r';
        color2='.r';
        lw=2;
        if (j>63)
            color='--k';
            color2='.k';
            lw=2;
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
    xticks([0 350 625 ])
    yticks([0 300 500  600])
    %axis equal
    set(gca,'FontSize', 30);  
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
    print(h,outputName,'-dpdf','-r0')

end
close all

!mv *pdf /Users/favinom/Desktop/wccm2021/figuresW/4realistic
%imshow(img)
%return
    
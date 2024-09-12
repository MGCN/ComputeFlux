clear all
close all

pa=[0 0.5;
    0.5 0;
    0.5 0.75;
    0.75 0.5;
    0.5 0.625;
    0.625 0.5;
    0 0.7;
    0.5 0];

pb=[1 0.5;
    0.5 1
    1 0.75;
    0.75 1;
    0.75 0.625;
    0.625 0.75;
    1 0.7;
    0.5 1];

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
    
    jstart=7;
    if(i==1)
        jstart=1;
    end
    
    for j=jstart:size(pa,1)
        
        color='r';
        color2='.r';
        lw=2.5;
        ms=25;
        if (j>6)
            color='--k';
            color2='.k';
            if (i>1)
                lw=2;
            else
                lw=1;
            end
            ms=15;
        end
        
        p1=pa(j,:);
        p2=pb(j,:);
        
        x=[p1(1) p2(1)];
        y=[p1(2) p2(2)];
        plot(x,y,color,'LineWidth',lw)
        
        plot(p1(1),p1(2),color2,'MarkerSize',ms)
        plot(p2(1),p2(2),color2,'MarkerSize',ms)
    end
    %axis equal
    set(gca,'ydir','normal');
    yticks([0 0.5 0.7 1])
    %axis equal
    set(gca,'FontSize', 30);  
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
    print(h,outputName,'-dpdf','-r0')

end
%close all

%!mv *pdf /Users/favinom/Desktop/wccm2021/figuresW/3complex

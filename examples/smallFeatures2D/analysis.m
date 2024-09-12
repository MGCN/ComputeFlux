clear all
close all

a=dir('*1.csv');

for i=1:length(a)
    
    name=a(i).name;
    newStr = split(name,'_');
    
    data{i}=csvread(name,1);
    
    pres{i}=data{i}(:,1);
    coor{i}=data{i}(:,3);
    
    res{i}=newStr{3};
    
end

start=1;

for i=start:length(a)
    
    plot(coor{i},pres{i});
    hold on
    
end
legend(res{start:end})
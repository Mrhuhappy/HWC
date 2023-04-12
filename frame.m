cmcc=importdata('./CMCC.txt');
tp1=cmcc;
mrna=importdata('./mrna_score.txt');
mirna=importdata('./mirna_score.txt');

Pro=importdata('./genename.txt');
a=mean(tp1);
b=mean(mrna);
c=mean(mirna);


number=10389;
s=1;


RWPOCPF2=zeros(number,1);
   for i=1:number
      for j=i+1:number
           if mrna(i,1)>mrna(j,1)-s/10 && mirna(i,1)>mirna(j,1)-s/10
              RWPOCPF2(i,1)=RWPOCPF2(i,1)+abs(mrna(i,1)-mrna(j,1))+abs(mirna(i,1)-mirna(j,1)); 
           end
           if mrna(i,1)-s/10 <mrna(j,1) && mirna(i,1)-s/10<mirna(j,1) 
              RWPOCPF2(j,1)=RWPOCPF2(j,1)+abs(mrna(j,1)-mrna(i,1))+abs(mirna(j,1)-mirna(i,1));  
           end
      end     
  end


value3 = zeros(number,1); 
inx3 = zeros(number,1);
RWPOCESnum3= zeros(6,1);
ESpro = importdata('./cancergene.txt');
filter =[50,100,150,200,250,300]; 
for i=1:1
  [value3(:,1),inx3(:,1)] = sort(RWPOCPF2(:,1),'descend'); 
end

for i = 1:6 
     RWPOCESnum3(i,1) = sum(ismember(Pro(inx3(1:filter(i),1)),ESpro));   
end
max2=max(RWPOCPF2);
tp2=(RWPOCPF2)/max2;

p=1;
RWPOCPF=zeros(number,11);
for k=0:10
   for i=1:number
      for j=i+1:number
           if tp1(i,1)>tp1(j,1)-p/1000 && tp2(i,1)>tp2(j,1)-p/1000
              RWPOCPF(i,k+1)=RWPOCPF(i,k+1)+(tp1(i,1)-tp1(j,1))*0.8+(tp2(i,1)-tp2(j,1))*(1-0.8); 
           end
           if tp1(i,1)-p/1000 <tp1(j,1) && tp2(i,1)-p/1000<tp2(j,1) 
              RWPOCPF(j,k+1)=RWPOCPF(j,k+1)+(tp1(j,1)-tp1(i,1))*0.8+(tp2(j,1)-tp2(i,1))*(1-0.8);  
           end
      end     
  end
end
value1 = zeros(number,11); 
inx1 = zeros(number,11);
RWPOCESnum = zeros(6,11);

for i=1:11
  [value1(:,i),inx1(:,i)] = sort(RWPOCPF(:,i),'descend'); 
end

filter =[50,100,150,200,250,300];
for j=1:11
  for i = 1:6 
     RWPOCESnum(i,j) = sum(ismember(Pro(inx1(1:filter(i),j)),ESpro)); 

  end
end

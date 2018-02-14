cd 'C:\Users\admins\Desktop\path\sim'

tic()
obj2=load('objmov5_extra.txt'); % the objects path points are now loaded
obj2=obj2';            
n=size(obj2,2);       %number of objects
m=size(obj2,1);       % number of points

%x=[ones(m,1),obj2(:,n)];
x=[ones(m,1),(obj2(:,n)./4)]; % Add a column of ones to x
theta = zeros(2, 1); % initialize fitting parameters
y = obj2(:,1:n-1);



stemp=zeros(m,n);
for i=1:(n)         % change of slopes are obtained
  for j=1:(m-1)
    stemp(j,i)=obj2(j+1,i)-obj2(j,i);
  end
end
    %to remember where the slopes are obtained
bp=zeros(m,n);
for i=1:(n)         
  for j=1:(m-1)
    bp(j,i)=round(abs(stemp(j+1,i)-stemp(j,i))*10);
    if(bp(j,i)!=0)                         % breakpoints are obtained where the bp value is 1
    bp(j,i)=1;
    endif
  end
end

count=0;                    %to remember the number of times the direction has changed
for(i=1:size(bp,1))
  if(bp(i,1)==1)
  count=count+1;
  endif
end
disp('number of times direction changes: ')
count
%number of segments = count+1
segnum=count+1;
disp('number of segments: ')
segnum
%segmenting and applying GD
k=0;
theta=zeros(2,1);
alpha=0.006
numit = 100000 %number of iterations
numit1= 850000;
h=x*theta; %X is Nx2 and theta is 2*1 hence, h is Nx1
ppath=[];
thetamatrix=[];

  theta=zeros(2,1);       %resetting theta to zero

  indexval=find(bp:1);
    xcord11=x(1:10,:);
    xcord1=repmat(xcord11,n-1,1);
    y11=y(1:10,:);
    m1=(size(y11,1)*size(y11,2));
      ycord1=y11(:);
    for n1=1:numit
      h=xcord1*theta;
      S=(h-ycord1).^2;
      J=(1/(2*m1))*sum(S);
      theta = theta - (alpha*(1/m1)*(xcord1'*(h-ycord1)));
    end
  ppath=[ppath;xcord1(1:10,:)*theta]; 
  ppath;
  thetamatrix=[thetamatrix;theta];
  
  
  theta=zeros(2,1);       %resetting theta to zero

  xcord12=x(10:15,:);
    xcord2=repmat(xcord12,n-1,1);
    y12=y(10:15,:);
    m2=(size(y12,1)*size(y12,2));
      ycord2=y12(:);
    for n1=1:numit
      h=xcord2*theta;
      S=(h-ycord2).^2;
      J=(1/(2*m2))*sum(S);
      theta = theta - (alpha*(1/m2)*(xcord2'*(h-ycord2)));
    end
  ppath=[ppath;xcord2(2:6,:)*theta]; 
  ppath;
  thetamatrix=[thetamatrix;theta];

    theta=zeros(2,1);       %resetting theta to zero

    xcord13=x(15:20,:);
    xcord3=repmat(xcord13,n-1,1);
    y13=y(15:20,:);
    m3=(size(y13,1)*size(y13,2));
      ycord3=y13(:);

    for n1=1:numit
      h=xcord3*theta;
      S=(h-ycord3).^2;
      J=(1/(2*m3))*sum(S);
      theta = theta - (alpha*(1/m3)*(xcord3'*(h-ycord3)));
    end
  ppath=[ppath;xcord3(2:6,:)*theta]; 
  thetamatrix=[thetamatrix;theta];
    theta=zeros(2,1);       %resetting theta to zero

  xcord14=x(20:30,:);
    xcord4=repmat(xcord14,n-1,1);
    y14=y(20:30,:);
    m4=(size(y14,1)*size(y14,2));
      ycord4=y14(:);
    for n1=1:numit
      h=xcord4*theta;
      S=(h-ycord4).^2;
      J=(1/(2*m4))*sum(S);
      theta = theta - (alpha*(1/m4)*(xcord4'*(h-ycord4)));
    end
  ppath=[ppath;xcord4(2:11,:)*theta]; 
  thetamatrix=[thetamatrix;theta];

    theta=zeros(2,1);       %resetting theta to zero

 

  xcord15=x(30:35,:);
    xcord5=repmat(xcord15,n-1,1);
    y15=y(30:35,:);
    m5=(size(y15,1)*size(y15,2));
      ycord5=y15(:);
    for n1=1:numit
      h=xcord5*theta;
      S=(h-ycord5).^2;
      J=(1/(2*m5))*sum(S);
      theta = theta - (alpha*(1/m5)*(xcord5'*(h-ycord5)));
    end
  ppath=[ppath;xcord5(2:6,:)*theta]; 
  thetamatrix=[thetamatrix;theta];
  
      theta=zeros(2,1);       %resetting theta to zero

  xcord16=x(35:40,:);
    xcord6=repmat(xcord16,n-1,1);
    y16=y(35:40,:);
    m6=(size(y16,1)*size(y16,2));
      ycord6=y16(:);
    for n1=1:numit
      h=xcord6*theta;
      S=(h-ycord6).^2;
      J=(1/(2*m6))*sum(S);
      theta = theta - (alpha*(1/m6)*(xcord6'*(h-ycord6)));
    end
  ppath=[ppath;xcord6(2:6,:)*theta]; 
  thetamatrix=[thetamatrix;theta];

      theta=zeros(2,1);       %resetting theta to zero

  xcord17=x(40:44,:);
    xcord7=repmat(xcord17,n-1,1);
    y17=y(40:44,:);
    m7=(size(y17,1)*size(y17,2));
      ycord7=y17(:);
    for n1=1:numit
      h=xcord7*theta;
      S=(h-ycord7).^2;
      J=(1/(2*m7))*sum(S);
      theta = theta - (alpha*(1/m7)*(xcord7'*(h-ycord7)));
    end
  ppath=[ppath;xcord7(2:5,:)*theta]; 
  thetamatrix=[thetamatrix;theta];


  thetamatrix

execution_time = toc()
number_of_samples = (n-1)
error_percentage=(1/m)*sum((obj2(:,1)-ppath).^2)*100
Accuracy_percentage = 100 - error_percentage

plot(obj2(:,1:(n-1)),'ko'); hold on; plot(ppath)   
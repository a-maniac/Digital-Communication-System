clc
clear all
%Carrier signal generation
Tb = 1; fc = 10;
t = 0:Tb/100:1;
c = sqrt(2/Tb)*sin(2*pi*fc*t);
%Message signal generation
N=8;
m = rand(1,N);
t1 = 0; t2 = Tb;
for i= 1:N
 t= [t1:0.01:t2]
 if(m(i)>0.5)
 m(i)=1;
 m_s = ones(1,length(t));
 else
 m(i)= 0;
 m_s = zeros(1,length(t));
 end
 message(i,:)= m_s;
 
%Product of message and carrier
 
ask_sig(i,:) = c.*m_s;
 t1 = t1+(Tb + 0.01);
 t2 = t2+(Tb + 0.01);
 
 %Plot message and ASK
 
subplot(5,1,2);
 axis([0 N -2 2]);
 plot(t,message(i,:),'r');
 title('Message signal');
 xlabel('m(t)');
 grid on
 hold on
 subplot(5,1,4);
 plot(t,ask_sig(i,:));
 title('ASK signal');
 xlabel('t--->');
 ylabel('s(t)');
 grid on
 hold on
end
hold off
%Plot the carrier signal and input binary bits
 subplot(5,1,3); plot(t,c);
 title('Carrier signal');
 xlabel('t--->');
 ylabel('c(t)');
 grid on
 subplot(5,1,1);
 stem(m)
 title('Binary data bits');
 xlabel('n--->');
 ylabel('b(n)');
 grid on
 
 %ASK demodulation
 t1=0; t2=Tb;
 for i= 1:N
 t=[t1:Tb/100:t2]
 %correlator
 x = sum(c.*ask_sig(i,:));
 if x>0
 demod(i)=1;
 else
 demod(i)=0;
 end
 end
 t1 = t1+(Tb + 0.01);
 t2 = t2+(Tb + 0.01);
 %Plot demodulated binary bits
 subplot(5,1,5);
 stem(demod)
 title('ASK demodulated signal');
 xlabel('n--->');
 ylabel('b(n)'); grid on

--------------------------------------------------------------------------------

%fsk


clc clear all
%Generate carrier signal tb=1;
fc1=2; fc2=5; t=0:tb/100:1;
c1=sqrt(2/tb)*sin(2*pi*fc1*t); c2=sqrt(2/tb)*sin(2*pi*fc2*t);
%Generate msg signal N=8;
m=rand(1,N); t1=0;
t2=tb; for i=1:N
t=[t1:(tb/100):t2] if m(i)>0.5
m(i)=1;
m_s=ones(1,length(t)); invm_s=zeros(1,length(t));
else
m(i)=0;
m_s=zeros(1,length(t)); invm_s=ones(1,length(t));
end message(i,:)=m_s;
%Product of carrier and msg fsk_sig1(i,:)=c1.*m_s; fsk_sig2(i,:)=c2.*invm_s; fsk=fsk_sig1+fsk_sig2;
 
%Plot msg and FSK subplot(3,2,1) axis([0 N -2 2]);
plot(t,message(i,:),'r');
title('Msg signal');xlabel('t--->');ylabel('m(t)'); grid on
hold on subplot(3,2,5); plot(t,fsk(i,:));
title('FSK signal');xlabel('t--->');ylabel('s(t)'); grid on
hold on t1=t1+(tb+0.1); t2=t2+(tb+0.1);
end hold off
%Plot carrier signal and input binary data subplot(3,2,3)
plot(t,c1);
title('Carrier signal 1');xlabel('t--->');ylabel('c1(t)'); grid on
subplot(3,2,4) plot(t,c2);
title('Carrier signal 2');xlabel('t--->');ylabel('c2(t)'); grid on
subplot(3,2,2) stem(m);
title('Binary data');xlabel('n--->');ylabel('b(t)'); grid on
%FSK Demodulation
 
t1=0;t2=tb; for i=1:N
t=[t1:tb/100:t2]
%Correlation x1=sum(c1.*fsk_sig1(i,:));
x2=sum(c2.*fsk_sig2(i,:)); x=x1-x2;
%decision device if x>0
demod(i)=1; else
demod(i)=0; end t1=t1+(tb+0.01); t2=t2+(tb+0.01);
end
%Plot demodulated output subplot(3,2,6) stem(demod);
title('FSK demodulated signal');xlabel('n--->');ylabel('b(t)'); grid on

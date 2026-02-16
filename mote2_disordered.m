% This code solves the gap equation matrix and search temperature (T_c) where he eigenvalue is 1.

format("long");

tic


% initializing squared frequency of insulator, omega_0^2
om0sqin=-1.0;
%om0sq=-0.0;        

%pl=parpool(18);

fileID=fopen('eigendata_Tvsm2_225_dis_kf1_mote2_alt35.dat','w');    % opening datafile to write
fprintf(fileID,'%1s %12s    %8s \n', '#',  'omega_0^2','Tc(K)');

for ilpomg=1:20      % loop for squared frequency of insulator omega_0^2

    disp(ilpomg);

om0sq=om0sqin+ilpomg*1.0;
%om0sq=0.01;  %meV



%setting the momentum integration and frequency sum
nsm=1000;
lambda=20;  % upper momentum cutoff
M=100;  % Matsubara frequencies


% initializing carrier density 
%nex=0;
%delnex=1;


%for ine=1:1       % carrier density loop
%nex=nex+delnex;
%nex=1;
%ned=nex*10^(16);    %carrier density in cm^-2

%disp(nex);

%parameters list1

kfermi=0.1*10^(8);  %Fermi momentum in cm^-1
alt=3.5;    % lattice constant in Angstrom         %change      
cbs=2*10^5*10^8;    %phonon velocity in  Angstrom/s  %change
hcut=6.582119*10^(-13);   %Planck constant in meV-s
kboltz=8.617333*10^(-5);    %Boltzmann constant in eV/K
melec=0.51099*10^6/(3*10^10)^2;      %mass of electron in eV/(cm/s)^2
%kfa=(3*pi^2*ned)^(1/3)*alt*10^(-8);  % 3D density      
kfa=kfermi*alt*10^(-8);  % 2D density      
melecef=2*melec;     % change  %%%%%%%%%%%%%%

%parameters list2
%kfermi=(2*pi*ned)^(1/2);  %Fermi momentum in cm^-1
efermi=hcut^2*10^(-6)*kfermi^2/(2.0*melecef); %Fermi energy in eV
efermi_red=efermi/kboltz;   %Fermi energy in K
vf=sqrt(2.0*efermi/melecef);   %Fermi velocity in cm/s
kavf=(vf/alt)*10^(8);          % v_F/a in 1/s
kbbyhct=(kboltz/hcut)*10^3;    % K_B/hcut in 1/(K-s)

%disp('Fermi energy in eV');
disp(efermi);

% DOS of STO
%ndos=0.524805*(kfa^1.33358)*10^(-3);  %meV^-1  for SrTiO3

ndos=kfermi*(alt*10^(-8))^2/(2*pi*vf*hcut);   %meV^-1    
%disp('ndos in meV-1');
%disp(ndos);

%melecef=(ndos*2*pi^2*hcut^2)*10^(-3)/((3*pi^2*ned)^(1/3)*(alt*10^(-8))^3); %eV/(cm/s)^2
%mratio=melecef/melec;  %effective mass


% defining bare coupling constant, gto
gbare=15;  %meV
D0=1;  %meV^-1
gbar=gbare^2*D0;        %meV



% defining phonon frequency at finite carrier density, 
%omegagsq=(om0sq+20.0*(ned/10^(20)));    %meV^2   %change
omegagsq=(om0sq);    %meV^2  

%  phonon frequency in the ordered state
if omegagsq < 0
omegagsq=2*abs(omegagsq);    % phonon frequency in the FE state in meV
end
omegag=sqrt(omegagsq);   % phonon frequency in the PE state in meV


% dimensionless "r" in the bosonic propagator
ar1=omegag^2*(alt/(hcut*cbs))^2;
rc=ar1;
disp(sqrt(rc));
%% effective Coupling Constant
ggm=gbar*ndos/(kfa);

%disp("effective coupling constant");
%disp(ggm);

%effective coupling constant in gap equation
fac=ggm;
facpi=ggm*4*kfa;



% initializing temperature
temp=0.0;      % K
  nt=10^10;
  delnt=0.000001;
ncnt=0;
ncnt1=0;
ncnt2=0;
ncnt3=0;
ncnt4=0;

  for ll=1:nt   % temperature Loop

%  AI reducing number of temperature iterations and simultaneously increasing the precision of the T_c value	  
	  
if ncnt==0
    delnt=delnt*10;
end
      
      if ncnt==1
          if ncnt1==0
      delnt = delnt/10.0;
    temp=0;
    ncnt1=ncnt1+1;
          end
     end
if ncnt==2
    if ncnt2==0
    temp=temp-delnt;
    delnt=delnt/10; 
    ncnt2=ncnt2+1;
    end
end
if ncnt==3
    if ncnt3==0
    temp=temp-delnt;
    delnt=delnt/10;
    ncnt3=ncnt3+1;
    end
end

if ncnt==4
    if ncnt4==0
    temp=temp-delnt;
    delnt=delnt/10;
    ncnt4=ncnt4+1;
    end
end

if ncnt==5
xfrmt=' %10.5f  %18.13f \n ';
fprintf(fileID,xfrmt,om0sq,temp-delnt);    % writing T_c vs. carrier density in datafile
disp(temp-delnt);
%disp((temp-delnt)*kboltz*10^3);
%disp(gbar^2*ndos/(2*pi));




break
end

tempold=temp;
  
if ncnt==0
    temp=delnt;
else
temp=temp+delnt;
end


% frequency rescaling

   fener=efermi_red;
mxkl=fix(0.5*(efermi_red/(pi*temp)-1));
ds=0.1;
if mxkl>M
   imll=fix(0.1*M);
        Mq=M;
       dimq=M+1;
aal=(log(mxkl-Mq)-log(ds))/((Mq-imll)*1.0);
abl=(imll*log(mxkl-Mq)-Mq*log(ds))/((Mq-imll)*1.0);
else
  imll=mxkl+1;
        Mq=mxkl;
       dimq=mxkl+1;
        aal=0;    
        abl=0;    
end



   % Initialize Matrix
   mat=zeros(dimq,dimq);


% Matrix formation
for ii=1:dimq      % frequency loop

 nopii=0;
 if ii >= (1+imll)
    nopii=1;
end
fk=(2*(ii-1)+1)*pi*temp+nopii*2*exp(aal*(ii-1)-abl)*pi*temp;    

sumq0=0.0;
for iq=1:dimq      % frequency loop
 nopiq=0;

if iq >= (1+imll)
    nopiq=1;
end

fq=(2*(iq-1)+1)*pi*temp+nopiq*2*exp(aal*(iq-1)-abl)*pi*temp;


if ne(ii,iq)==1
dfnqp=0.0;
dfnqm=0.0;
dh=1.0/nsm;
qlm=0.000001;

for irr=1:nsm+1  % momentum integration of bosonic function

%bpropm=lambda/(rc+qlm^2*lambda^(2)+(fq-fk)^2*(kbbyhct*alt/cbs)^2+facpi*(abs(fq-fk)/qlm)*(kbbyhct/kavf)*lambda^(-1));
%bpropp=lambda/(rc+qlm^2*lambda^(2)+(fq+fk)^2*(kbbyhct*alt/cbs)^2+facpi*(abs(fq+fk)/qlm)*(kbbyhct/kavf)*lambda^(-1));

bpropm=lambda/(rc+qlm^2*lambda^(2)+facpi*(abs(fq-fk)/qlm)*(kbbyhct/kavf)*lambda^(-1));
bpropp=lambda/(rc+qlm^2*lambda^(2)+facpi*(abs(fq+fk)/qlm)*(kbbyhct/kavf)*lambda^(-1));

if irr==1
    dfnqp=dfnqp+bpropp;
    dfnqm=dfnqm+bpropm;
elseif irr==nsm+1
    dfnqp=dfnqp+bpropp;
    dfnqm=dfnqm+bpropm;
else
if mod(irr,2)==0
    dfnqp=dfnqp+bpropp*4.0;
    dfnqm=dfnqm+bpropm*4.0;
else
    dfnqp=dfnqp+bpropp*2.0;
    dfnqm=dfnqm+bpropm*2.0;
end
end
qlm=qlm+dh;
end    % momentum  integration ends

if  iq < (imll+1)
sumq0=sumq0+(-dfnqp+dfnqm)*dh/3.0;
else
sumq0=sumq0+((-dfnqp+dfnqm)*dh/3.0)*(1+exp(aal*(iq-1)-abl)*(exp(aal)-1));
end

end

end
sumq0=1+fac*temp*sumq0/fk;        %denominator in the gap equation    


for j=1:dimq       % frequency loop


nopj=0;
if j >= (1+imll)
    nopj=1;
end
fp=(2*(j-1)+1)*pi*temp+nopj*2*exp(aal*(j-1)-abl)*pi*temp;

dfnp=0.0;
dfnm=0.0;
dh=1.0/nsm;
plm=0.000001;

for ir=1:nsm+1  % momentum integration bosonic function

%bpropm=lambda/(rc+plm^2*lambda^(2)+(fp-fk)^2*(kbbyhct*alt/cbs)^2+facpi*(abs(fp-fk)/plm)*(kbbyhct/kavf)*lambda^(-1));
%bpropp=lambda/(rc+plm^2*lambda^(2)+(fp+fk)^2*(kbbyhct*alt/cbs)^2+facpi*(abs(fp+fk)/plm)*(kbbyhct/kavf)*lambda^(-1));

bpropm=lambda/(rc+plm^2*lambda^(2)+facpi*(abs(fp-fk)/plm)*(kbbyhct/kavf)*lambda^(-1));
bpropp=lambda/(rc+plm^2*lambda^(2)+facpi*(abs(fp+fk)/plm)*(kbbyhct/kavf)*lambda^(-1));


if ir==1
    dfnp=dfnp+bpropp;
    dfnm=dfnm+bpropm;
elseif ir==nsm+1
    dfnp=dfnp+bpropp;
    dfnm=dfnm+bpropm;
else
if mod(ir,2)==0
    dfnp=dfnp+bpropp*4.0;
    dfnm=dfnm+bpropm*4.0;
else
    dfnp=dfnp+bpropp*2.0;
    dfnm=dfnm+bpropm*2.0;
end
end
plm=plm+dh;
end    % momentum integration ends

dfnp=dfnp*dh/3.0;
dfnm=dfnm*dh/3.0;


% matrix elements
matelem=0.0;
if ne(ii,j)==1
matelem=matelem+((dfnm+dfnp)/fp)/sumq0;
end

      if  j < (imll+1)
               mat(ii,j)=fac*matelem*temp;
              else
              mat(ii,j)=fac*matelem*temp*(1+exp(aal*(j-1)-abl)*(exp(aal)-1));    %matrix(k0 x p0) elements
       end


        end
        end
%mat;
[V,D]=eigs(mat,1);
Y=[temp,D];
if D<1
    ncnt=ncnt+1;
end
%disp(Y);
    

  
  
  end     %  temperature loop end
 
  
%end   % density loop end  
  
  
  

end
fclose(fileID);

%delete(pl)

toc

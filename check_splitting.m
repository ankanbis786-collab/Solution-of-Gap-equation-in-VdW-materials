% This code solves the gap equation matrix and search temperature (T_c) where he eigenvalue is 1.

format("long");

tic


% initializing squared frequency of insulator, omega_0^2
om0sqin=1;
%om0sq=-10.0;        

%pl=parpool(18);

fileID=fopen('eigendata_loop1_tb_kf_59_b05_check1_Tcb_chek_ef.dat','w');    % opening datafile to write
%fprintf(fileID,'%1s %12s   %15s  %8s \n', '#',  'omega_0^2', 'delta/efermi0','Tc(K)');
fprintf(fileID,'%1s %12s   %15s  %15s  %10s \n', '#', 'omega_0^2', 'delta/efermi0', 'Tc(K)', 'rc2');

for ilpomg=1:20     % loop for squared frequency of insulator omega_0^2

    disp(ilpomg);

om0sq=om0sqin-ilpomg*1.0;
%om0sq=-1;  %meV



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

alt=3.9;    % lattice constant in Angstrom         %change      
cbs=7.5*10^5*10^8;    %phonon velocity in  Angstrom/s  %change
hcut=6.582119*10^(-13);   %Planck constant in meV-s
kboltz=8.617333*10^(-5);    %Boltzmann constant in eV/K
melec=0.51099*10^6/(3*10^10)^2;      %mass of electron in eV/(cm/s)^2
melecef=2*melec;     % change  %%%%%%%%%%%%%%
kbbyhct=(kboltz/hcut)*10^3;    % K_B/hcut in 1/(K-s)

omegagsq=(om0sq);    %meV^2  


if omegagsq < 0
omegagsq=2*abs(omegagsq);    % phonon frequency in the FE state in meV
end
omegag=sqrt(omegagsq);   % phonon frequency in the PE state in meV


% dimensionless "r" in the bosonic propagator
ar1=omegag^2*(alt/(hcut*cbs))^2;
rc=ar1;
rc2=rc^2;
%disp(rc2);

% defining bare coupling constant, gto
gbare=25;  %meV
D0=1;  %meV^-1
gbar=gbare^2*D0;        %meV

bq=0.05; % dimensionless
delta=gbare*sqrt(rc/bq)*10^(-3);      %eV
disp(delta);
kfermi0=0.59*10^(8);  %Fermi momentum in cm^-1    
efermi0=hcut^2*10^(-6)*kfermi0^2/(2.0*melecef); %Fermi energy in eV
disp(efermi0);
efermi_red=efermi0/kboltz;   %Fermi energy in K

%efermi=efermi-delta;   % eV
%kfermi=0.59*10^(8);  %Fermi momentum in cm^-1
%kfa=kfermi*alt*10^(-8);  % 2D density      
%efermi=hcut^2*10^(-6)*kfermi^2/(2.0*melecef); %Fermi energy in eV


efermi=efermi0-delta/2;   % eV
kfermi=sqrt(2.0*melecef*efermi)/(hcut*10^(-3));  %Fermi momentum in cm^-1
kfa=kfermi*alt*10^(-8);  % 2D density 
%efermi_red=efermi/kboltz;   %Fermi energy in K
%disp("cutoff");
%disp(efermi_red);
vf=sqrt(2.0*efermi/melecef);   %Fermi velocity in cm/s
kavf=(vf/alt)*10^(8);          % v_F/a in 1/s
ndos=kfermi*(alt*10^(-8))^2/(2*pi*vf*hcut);   %meV^-1    

%delta=0.035*10^(8);  %Fermi momentum in cm^-1
%bq=0.5; % dimensionless
%delta=(gbare/hcut)*sqrt(rc/bq)/(cbs*10^(-8));
%disp(delta);

efermi2=efermi0+delta/2;   % eV
kfermi2=sqrt(2.0*melecef*efermi2)/(hcut*10^(-3)); 
%kfermi2=kfermi-delta;  %Fermi momentum in cm^-1
kfa2=kfermi2*alt*10^(-8);  % 2D density
%efermi2=hcut^2*10^(-6)*kfermi2^2/(2.0*melecef); %Fermi energy in eV
vf2=sqrt(2.0*efermi2/melecef);   %Fermi velocity in cm/s
kavf2=(vf2/alt)*10^(8);          % v_F/a in 1/s
ndos2=kfermi2*(alt*10^(-8))^2/(2*pi*vf2*hcut);   %meV^-1    

%disp([ndos,ndos2]);
%disp([kfermi,kfermi2]);




%% effective Coupling Constant
ggm=gbar*ndos/(kfa);

%disp("effective coupling constant");
%disp(ggm);

%effective coupling constant in gap equation
fac=ggm;
facpi=gbar*(ndos/kavf+ndos2/kavf2);

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
%xfrmt=' %10.5f     %18.10f   %18.13f \n ';
%fprintf(fileID,xfrmt,om0sq,delta/efermi0,temp-delnt);    % writing T_c vs. carrier density in datafile
xfrmt = ' %10.5f   %18.10f   %18.13f   %12.6f \n';
fprintf(fileID, xfrmt, om0sq, delta/efermi0, temp-delnt, rc2);
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
%disp(mxkl);
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

bpropm=lambda/(rc+qlm^2*lambda^(2)+2*facpi*(abs(fq-fk)/qlm)*(kbbyhct)*lambda^(-1));
bpropp=lambda/(rc+qlm^2*lambda^(2)+2*facpi*(abs(fq+fk)/qlm)*(kbbyhct)*lambda^(-1));

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

bpropm=lambda/(rc+plm^2*lambda^(2)+2*facpi*(abs(fp-fk)/plm)*(kbbyhct)*lambda^(-1));
bpropp=lambda/(rc+plm^2*lambda^(2)+2*facpi*(abs(fp+fk)/plm)*(kbbyhct)*lambda^(-1));


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

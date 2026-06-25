program SWs_CMC42_Yin_Yang_GSP04_3rd
!********************************************************************************
! This program was written by Rasoul Mirzaei Shiri.
!
! It solves the shallow water equations on a Yin-Yang grid (MYYIC) using a
! fourth-order MacCormack scheme with Runge-Kutta time integration
! and the third stabilization strategy.
!
! The code is developed and employed to solve the Unstable Jet (Gallewsky et al., 2004) test case
! for validation and assessment of the numerical method.
!********************************************************************************
implicit none
integer,parameter::nx=256,ny=(nx/2)+1,kx=512,ky=(kx/2)+1,n_epsilon_x=nx/8,n_epsilon_y=(ny-1)/8&
,mx=(3*nx/4)+2*n_epsilon_x,my=((ny-1)/2)+1+2*n_epsilon_y,kky=2*(ny-1),kkx=nx/2,days=8
integer::i,l,n,j,knex,kney,k_cof
REAL*8,parameter::delt=160.0d0
REAL*8::a2,a3,a4,t,lx,ly,pi,twopi,dlanda,dphi,top_x_limit,top_y_limit&
,aa,a12,a13,dtime,bbbb,anoo
REAL::ta(2)
REAL*8::hene(nx,ny),pne(ny),lne(nx),u2(nx,ny),h000,&
omega,zita(nx,ny),delta(nx,ny),v2(nx,ny)&
,uphi(ny),hcapit(ny),hprime(nx,ny)&
,uphiin(ky),hcapitin(ky),hprimein(kx,ky),u1ein(kx,ky),v1ein(kx,ky),hprimeein(kx,ky),hcapitein(kx,ky)&
,u1(nx,ny),v1(nx,ny),uphie(kx,ky)
real*8::hen1(mx,my),hen2(mx,my),hee1(mx,my),hee2(mx,my)&
,u2e(mx,my),v2e(mx,my),ln(mx),pn(my),le(mx),pe(my),fon(mx,my),foe(mx,my),u2n(mx,my),v2n(mx,my)&
,hcapitn(mx,my),hprimen(mx,my),u1n(mx,my),v1n(mx,my),hcapite(mx,my),hprimee(mx,my),u1e(mx,my),v1e(mx,my)&
,zitan(mx,my),deltan(mx,my),zitae(mx,my),deltae(mx,my)&
,dux(nx,ny),duy(nx,ny),dvx(nx,ny),dvy(nx,ny),dhx(nx,ny),dhy(nx,ny),hun1(mx,my),hvn1(mx,my)&
,dunx(mx,my),duny(mx,my),dvnx(mx,my),dvny(mx,my),dhnx(mx,my),dhny(mx,my)&
,duex(mx,my),duey(mx,my),dvex(mx,my),dvey(mx,my),dhex(mx,my),dhey(mx,my)&
,hhn1(mx,my),hhun1(mx,my),hhvn1(mx,my),hhhn1(mx,my),hue1(mx,my),hve1(mx,my),hhe1(mx,my),hhue1(mx,my)&
,hhve1(mx,my),hhhe1(mx,my),hun2(mx,my),hvn2(mx,my)&
,hhn2(mx,my),hhun2(mx,my),hhvn2(mx,my),hhhn2(mx,my),hue2(mx,my),hve2(mx,my),hhe2(mx,my),hhue2(mx,my)&
,hhve2(mx,my),hhhe2(mx,my),hun3(mx,my),hvn3(mx,my)&
,hhn3(mx,my),hhun3(mx,my),hhvn3(mx,my),hhhn3(mx,my),hue3(mx,my),hve3(mx,my),hhe3(mx,my),hhue3(mx,my)&
,hhve3(mx,my),hhhe3(mx,my),hun4(mx,my),hvn4(mx,my)&
,hhn4(mx,my),hhun4(mx,my),hhvn4(mx,my),hhhn4(mx,my),hue4(mx,my),hve4(mx,my),hhe4(mx,my),hhue4(mx,my)&
,hhve4(mx,my),hhhe4(mx,my),hhu(nx,ny),hhv(nx,ny),hhh(nx,ny)

character(5)::int_method
!int_method='bilin'
int_method='bicub'
open(10,file='yin_height.plt')
open(20,file='yang_height.plt')
open(30,file='yin_yang_height.plt')

open(50,file='comptime.dat')

open(100,file='yin_zita.plt')
open(110,file='yang_zita.plt')
open(120,file='yin_yang_zita.plt')
open(130,file='yin_delta.plt')
open(140,file='yang_delta.plt')
open(150,file='yin_yang_delta.plt')

open(250,file='yin_u.plt')
open(260,file='yang_u.plt')
open(270,file='yin_yang_u.plt')
open(280,file='yin_v.plt')
open(290,file='yang_v.plt')
open(300,file='yin_yang_v.plt')

knex=(nx-mx)/2
kney=(ny-my)/2
k_cof=kx/nx

print*,knex,kney,k_cof
print*,mx,my
print*,n_epsilon_x,n_epsilon_y
aa=6.37122d6
pi=2.0d0*dasin(1.0d0);twopi=2.0d0*pi
lx=twopi;ly=pi
dlanda=lx/dble(nx);dphi=ly/dble(ny-1)
t=dble(days)*86400.0d0;l=int(t/delt)
anoo=0.0d0
!anoo=1.0d5
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)
!********************************************************************************************************************
!Providing initial conditions
!********************************************************************************************************************
call initial_conditions(kx,ky,uphiin,hcapitin,hprimein,h000)

do j=1,ny
hcapit(j)=hcapitin(k_cof*(j-1)+1)
uphi(j)=uphiin(k_cof*(j-1)+1)

do i=1,nx
u1(i,j)=uphi(j)
v1(i,j)=0.0d0
hprime(i,j)=hprimein(k_cof*(i-1)+1,k_cof*(j-1)+1)
hene(i,j)=hcapit(j)+hprime(i,j)
enddo
enddo

do i=1,nx
ln(i)=-top_x_limit+dble(i-1)*dlanda
enddo
omega=7.292d-5
do j=1,my
pn(j)=dble(j-1)*dphi-top_y_limit
do i=1,nx
fon(i,j)=2.0d0*omega*dsin(pn(j))
enddo
enddo


do i=1,nx
le(i)=-top_x_limit+dble(i-1)*dlanda
enddo


do j=1,my
pe(j)=dble(j-1)*dphi-top_y_limit
do i=1,nx
a12=pe(j)
a13=le(i)
foe(i,j)=-2.0d0*omega*dcos(a12)*dsin(a13)
enddo
enddo
call initial_yang_total(kx,ky,uphie,hcapitein,hprimeein,h000,u1ein,v1ein)

do j=1,my

do i=1,nx
hcapite(i,j)=hcapitein(k_cof*(i-1)+1,k_cof*(j+kney-1)+1)
hprimee(i,j)=hprimeein(k_cof*(i-1)+1,k_cof*(j+kney-1)+1)
u1e(i,j)=u1ein(k_cof*(i-1)+1,k_cof*(j+kney-1)+1)
v1e(i,j)=v1ein(k_cof*(i-1)+1,k_cof*(j+kney-1)+1)
hee1(i,j)=hcapite(i,j)+hprimee(i,j)
hcapitn(i,j)=hcapit(j+kney)
hprimen(i,j)=hprime(i,j+kney)
u1n(i,j)=uphi(j+kney)
v1n(i,j)=0.0d0
hen1(i,j)=hcapitn(i,j)+hprimen(i,j)
enddo
enddo
do i=1,nx
lne(i)=-pi+(i-1)*dlanda
enddo
do j=1,ny
pne(j)=(j-1)*dphi-(pi/2.0d0)
enddo
!********************************************************************************************************************
!computing Vorticity and Divergence
!********************************************************************************************************************
call derivative3(u1,v1,u1e,v1e,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dunx,dvnx,duny&
,dvny,duex,dvex,duey,dvey,int_method,n_epsilon_x,n_epsilon_y)
 
call zita_delta(u1n,v1n,dunx,dvnx,duny,dvny,pn,mx,my,zitan,deltan)
call zita_delta(u1e,v1e,duex,dvex,duey,dvey,pe,mx,my,zitae,deltae)
call zita_delta(u1,v1,dux,dvx,duy,dvy,pne,nx,ny,zita,delta)

!********************************************************************************************************************
!merging zita and delta on the complete Yin-Yang grid
!********************************************************************************************************************
!call merging_Yin_Yang(zitan,zitae,zita,nx,ny,my,n_epsilon,dphi,dlanda,int_method)
!call merging_Yin_Yang(deltan,deltae,delta,nx,ny,my,n_epsilon,dphi,dlanda,int_method)

!********************************************************************************************************************
!Print data to the files
!********************************************************************************************************************


write(10,*)'zone            i=',nx,'   j=',my
write(20,*)'zone            i=',my,'   j=',nx
write(30,*)'zone            i=',nx,'   j=',ny

write(100,*)'zone            i=',nx,'   j=',my
write(110,*)'zone            i=',my,'   j=',nx
write(120,*)'zone            i=',nx,'   j=',ny
write(130,*)'zone            i=',nx,'   j=',my
write(140,*)'zone            i=',my,'   j=',nx
write(150,*)'zone            i=',nx,'   j=',ny

write(250,*)'zone            i=',nx,'   j=',my
write(260,*)'zone            i=',my,'   j=',nx
write(270,*)'zone            i=',nx,'   j=',ny
write(280,*)'zone            i=',nx,'   j=',my
write(290,*)'zone            i=',my,'   j=',nx
write(300,*)'zone            i=',nx,'   j=',ny


do j=1,my
do i=1,nx
write(10,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,hen1(i,j)
write(100,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,zitan(i,j)
write(130,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,deltan(i,j)
write(250,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,u1n(i,j)
write(280,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,v1n(i,j)
enddo
enddo
do i=1,nx
do j=1,my
write(20,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,hee1(i,j)
write(110,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,zitae(i,j)
write(140,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,deltae(i,j)
write(260,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,u1e(i,j)
write(290,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,v1e(i,j)
enddo
enddo
do j=1,ny
do i=1,nx
write(30,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,hene(i,j)
write(120,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,zita(i,j)
write(150,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,delta(i,j)
write(270,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,u1(i,j)
write(300,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,v1(i,j)
enddo
enddo



!********************************************************************************************************************
!Time marching
!********************************************************************************************************************
a2=0.5d0;a3=0.5d0;a4=1.0d0

do n=1,l
!********************************************************************************************************************
!First step
!********************************************************************************************************************
!Yin component

call derivative1(hene,u1,v1,hee1,u1e,v1e,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)

call RK_init_swe(u1n,v1n,hen1,dunx,duny,dvnx,dvny,dhnx,dhny,mx,my,a2,hun1,hvn1,hhn1,hhun1,hhvn1,hhhn1,pn&
,fon,delt,u1n,v1n,hen1)


!Yang component

call RK_init_swe(u1e,v1e,hee1,duex,duey,dvex,dvey,dhex,dhey,mx,my,a2,hue1,hve1,hhe1,hhue1,hhve1,hhhe1,pe&
,foe,delt,u1e,v1e,hee1)
call merging_vector_Yin_Yang(hhhn1,hhhe1,hhh,hhun1,hhue1,hhu,hhvn1,hhve1,hhv,nx,ny,mx,my,n_epsilon_x&
,n_epsilon_y,dphi,dlanda,int_method)
!call decompose_yang(hhh,hhu,hhv,hhhe1,hhue1,hhve1,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
!********************************************************************************************************************
!Second step
!********************************************************************************************************************

 !Yin component
call derivative2(hhh,hhu,hhv,hhhe1,hhue1,hhve1,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)

call RK_init_swe(u1n,v1n,hen1,dunx,duny,dvnx,dvny,dhnx,dhny,mx,my,a3,hun2,hvn2,hhn2,hhun2,hhvn2,hhhn2,pn&
,fon,delt,hhun1,hhvn1,hhhn1)

!Yang component

call RK_init_swe(u1e,v1e,hee1,duex,duey,dvex,dvey,dhex,dhey,mx,my,a3,hue2,hve2,hhe2,hhue2,hhve2,hhhe2,pe,foe&
,delt,hhue1,hhve1,hhhe1)
call merging_vector_Yin_Yang(hhhn2,hhhe2,hhh,hhun2,hhue2,hhu,hhvn2,hhve2,hhv,nx,ny,mx,my,n_epsilon_x&
,n_epsilon_y,dphi,dlanda,int_method)
!call decompose_yang(hhh,hhu,hhv,hhhe2,hhue2,hhve2,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
!********************************************************************************************************************
!Third step
!********************************************************************************************************************

 !Yin component
 call derivative1(hhh,hhu,hhv,hhhe2,hhue2,hhve2,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)

call RK_init_swe(u1n,v1n,hen1,dunx,duny,dvnx,dvny,dhnx,dhny,mx,my,a4,hun3,hvn3,hhn3,hhun3,hhvn3,hhhn3,pn&
,fon,delt,hhun2,hhvn2,hhhn2)

!Yang component

call RK_init_swe(u1e,v1e,hee1,duex,duey,dvex,dvey,dhex,dhey,mx,my,a4,hue3,hve3,hhe3,hhue3,hhve3,hhhe3&
,pe,foe,delt,hhue2,hhve2,hhhe2)
call merging_vector_Yin_Yang(hhhn3,hhhe3,hhh,hhun3,hhue3,hhu,hhvn3,hhve3,hhv,nx,ny,mx,my,n_epsilon_x&
,n_epsilon_y,dphi,dlanda,int_method)
!call decompose_yang(hhh,hhu,hhv,hhhe3,hhue3,hhve3,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
!********************************************************************************************************************
!Fourth step
!********************************************************************************************************************

 !Yin component
call derivative2(hhh,hhu,hhv,hhhe3,hhue3,hhve3,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)


call RK_init_swe(u1n,v1n,hen1,dunx,duny,dvnx,dvny,dhnx,dhny,mx,my,a4,hun4,hvn4,hhn4,hhun4,hhvn4,hhhn4,pn&
,fon,delt,hhun3,hhvn3,hhhn3)
call RK4_final(hen1,hen2,hhn1,hhn2,hhn3,hhn4,nx,my)
call RK4_final(u1n,u2n,hun1,hun2,hun3,hun4,nx,my)
call RK4_final(v1n,v2n,hvn1,hvn2,hvn3,hvn4,nx,my)
!Yang component

call RK_init_swe(u1e,v1e,hee1,duex,duey,dvex,dvey,dhex,dhey,mx,my,a4,hue4,hve4,hhe4,hhue4,hhve4,hhhe4,pe&
,foe,delt,hhue3,hhve3,hhhe3)
call RK4_final(hee1,hee2,hhe1,hhe2,hhe3,hhe4,nx,my)
call RK4_final(u1e,u2e,hue1,hue2,hue3,hue4,nx,my)
call RK4_final(v1e,v2e,hve1,hve2,hve3,hve4,nx,my)
!********************************************************************************************************************
!Filtering the noises that are coming from boundaries
!********************************************************************************************************************
call merging_vector_Yin_Yang(hen2,hee2,hene,u2n,u2e,u2,v2n,v2e,v2,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
!call Filter_boundaries(u2n,u2e,v2n,v2e,hen2,hee2,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
call decompose_yang(hene,u2,v2,hee2,u2e,v2e,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
!********************************************************************************************************************
!damping numerical instability
!********************************************************************************************************************
!Yin grid
!if(dabs(anoo)>=1.0d-1)then
!derivative4(hene,u1,v1,u1e,v1e,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
!,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)



!call provide_Yin_vector_to_derive(hen1,hee1,u1n,v1n,u1e,v1e,dlanda,dphi,nx,ny,my,n_epsilon&
!,h_meridion,u_meridion,ucosphi,v_meridion,vcosphi,int_method)

!call provide_cos_phi(dphi,ny,cosphi)

!call laplasian_longitude(hen1,nx,my,dlanda,dccx,pn)
!call laplasian_latitude(h_meridion,nx,2*(ny-1),dphi,dccy,cosphi)
!do j=1,my
!do i=1,nx
!hen2(i,j)=hen2(i,j)+anoo*(dccx(i,j)+dccy(i,j+kney))*delt
!enddo
!enddo
!call laplasian_longitude(u1n,nx,my,dlanda,dccx,pn)
!call laplasian_latitude(u_meridion,nx,2*(ny-1),dphi,dccy,cosphi)
!do j=1,my
!do i=1,nx
!u2n(i,j)=u2n(i,j)+anoo*(dccx(i,j)+dccy(i,j+kney))*delt
!enddo
!enddo
!call laplasian_longitude(v1n,nx,my,dlanda,dccx,pn)
!call laplasian_latitude(v_meridion,nx,2*(ny-1),dphi,dccy,cosphi)
!do j=1,my
!do i=1,nx
!v2n(i,j)=v2n(i,j)+anoo*(dccx(i,j)+dccy(i,j+kney))*delt
!enddo
!enddo
!***********************
!Yang grid
!call provide_Yang_vector_to_derive(hen1,hee1,u1n,v1n,u1e,v1e,dlanda,dphi,nx,ny,my,n_epsilon&
!,h_meridion,u_meridion,ucosphi,v_meridion,vcosphi,int_method)


!call laplasian_longitude(hee1,nx,my,dlanda,dccx,pe)
!call laplasian_latitude(h_meridion,nx,2*(ny-1),dphi,dccy,cosphi)
!do j=1,my
!do i=1,nx
!hee2(i,j)=hee2(i,j)+anoo*(dccx(i,j)+dccy(i,j+kney))*delt
!enddo
!enddo
!call laplasian_longitude(u1e,nx,my,dlanda,dccx,pe)
!!call laplasian_latitude(u_meridion,nx,2*(ny-1),dphi,dccy,cosphi)
!do j=1,my
!do i=1,nx
!u2e(i,j)=u2e(i,j)+anoo*(dccx(i,j)+dccy(i,j+kney))*delt
!enddo
!enddo
!call laplasian_longitude(v1e,nx,my,dlanda,dccx,pe)
!call laplasian_latitude(v_meridion,nx,2*(ny-1),dphi,dccy,cosphi)
!do j=1,my
!do i=1,nx
!v2e(i,j)=v2e(i,j)+anoo*(dccx(i,j)+dccy(i,j+kney))*delt
!enddo
!enddo
!endif
!********************************************************************************************************************
!Computing Vorticity and Divergence
!********************************************************************************************************************

if(dabs(dble(n)*delt/21600.0d0-dble(floor(dble(n)*delt/21600.0d0)))<=1.0d-8)then

call derivative3(u2,v2,u2e,v2e,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dunx,dvnx,duny&
,dvny,duex,dvex,duey,dvey,int_method,n_epsilon_x,n_epsilon_y)

call zita_delta(u2n,v2n,dunx,dvnx,duny,dvny,pn,mx,my,zitan,deltan)
call zita_delta(u2e,v2e,duex,dvex,duey,dvey,pe,mx,my,zitae,deltae)
call zita_delta(u2,v2,dux,dvx,duy,dvy,pne,nx,ny,zita,delta)



!********************************************************************************************************************
!Merging Data
!********************************************************************************************************************
!call merging_Yin_Yang(zitan,zitae,zita,nx,ny,my,n_epsilon,dphi,dlanda,int_method)
!call merging_Yin_Yang(deltan,deltae,delta,nx,ny,my,n_epsilon,dphi,dlanda,int_method)
!call merging_Yin_Yang(hen2,hee2,hene,nx,ny,my,n_epsilon,dphi,dlanda,int_method)
!call merging_vector_Yin_Yang(hee2,u2n,u2e,u2,v2n,v2e,v2,nx,ny,my,n_epsilon,dphi,dlanda,int_method)
!call merging_vector_Yin_Yang(hen2,hee2,hene,u2n,u2e,u2,v2n,v2e,v2,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
!********************************************************************************************************************
!printing Data
!********************************************************************************************************************
write(10,*)'zone            i=',nx,'   j=',my
write(20,*)'zone            i=',my,'   j=',nx
write(30,*)'zone            i=',nx,'   j=',ny

write(100,*)'zone            i=',nx,'   j=',my
write(110,*)'zone            i=',my,'   j=',nx
write(120,*)'zone            i=',nx,'   j=',ny
write(130,*)'zone            i=',nx,'   j=',my
write(140,*)'zone            i=',my,'   j=',nx
write(150,*)'zone            i=',nx,'   j=',ny

write(250,*)'zone            i=',nx,'   j=',my
write(260,*)'zone            i=',my,'   j=',nx
write(270,*)'zone            i=',nx,'   j=',ny
write(280,*)'zone            i=',nx,'   j=',my
write(290,*)'zone            i=',my,'   j=',nx
write(300,*)'zone            i=',nx,'   j=',ny

do j=1,my
do i=1,nx
write(10,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,hen2(i,j)
write(100,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,zitan(i,j)
write(130,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,deltan(i,j)
write(250,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,u2n(i,j)
write(280,*)ln(i)*180.0d0/pi,pn(j)*180.0d0/pi,v2n(i,j)
enddo
enddo
do i=1,nx
do j=1,my
write(20,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,hee2(i,j)
write(110,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,zitae(i,j)
write(140,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,deltae(i,j)
write(260,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,u2e(i,j)
write(290,*)pe(j)*180.0d0/pi,le(i)*180.0d0/pi,v2e(i,j)
enddo
enddo
do j=1,ny
do i=1,nx
write(30,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,hene(i,j)
write(120,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,zita(i,j)
write(150,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,delta(i,j)
write(270,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,u2(i,j)
write(300,*)lne(i)*180.0d0/pi,pne(j)*180.0d0/pi,v2(i,j)
enddo
enddo
endif



!********************************************************************************************************************
!Changing and Converting
!********************************************************************************************************************


print*,n*delt,hen2(1,1),u2n(2,1)

do i=1,mx

do j=1,my

hen1(i,j)=hen2(i,j)
hee1(i,j)=hee2(i,j)
u1n(i,j)=u2n(i,j)
u1e(i,j)=u2e(i,j)
v1n(i,j)=v2n(i,j)
v1e(i,j)=v2e(i,j)

enddo
enddo

do j=1,ny
do i=1,nx
u1(i,j)=u2(i,j)
v1(i,j)=v2(i,j)

enddo
enddo 
enddo

bbbb=dtime(ta)
write(50,*)'computational time=',bbbb
end program SWs_CMC42_Yin_Yang_GSP04_3rd
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dbx42p(e,dbx,m,k,dx)
integer::i,j,m,k
real*8::e(m,k),dbx(m,k),a,dx,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do j=1,k
sum=(e(1,j)-e(m,j))/dx
do i=m,2,-1
sum=sum+((e(i,j)-e(i-1,j))/dx)*((a/(a-1.0d0))**(m-i+1))
enddo
dbx(1,j)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(m-1))*a))
do i=2,m
dbx(i,j)=(((e(i,j)-e(i-1,j))/dx)-(a*dbx(i-1,j)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dfx42p(e,dfx,m,k,dx)
integer::i,j,m,k
real*8::e(m,k),dfx(m,k),a,dx,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do j=1,k
sum=(e(1,j)-e(m,j))/dx
do i=2,m
sum=sum+((e(i,j)-e(i-1,j))/dx)*((a/(a-1.0d0))**(i-1))
enddo
dfx(m,j)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(m-1))*a))
do i=m-1,1,-1
dfx(i,j)=(((e(i+1,j)-e(i,j))/dx)-(a*dfx(i+1,j)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dfy42p(e,dfy,m,k,dy)
integer::i,j,m,k
real*8::e(m,k),dfy(m,k),a,dy,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do i=1,m
sum=(e(i,1)-e(i,k))/dy
do j=2,k
sum=sum+((e(i,j)-e(i,j-1))/dy)*((a/(a-1.0d0))**(j-1))
enddo
dfy(i,k)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(k-1))*a))
do j=k-1,1,-1
dfy(i,j)=(((e(i,j+1)-e(i,j))/dy)-(a*dfy(i,j+1)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dby42p(e,dby,m,k,dy)
integer::i,j,m,k
real*8::e(m,k),dby(m,k),a,dy,sum
a=0.5d0-0.5d0/dsqrt(3.0d0)
do i=1,m
sum=(e(i,1)-e(i,k))/dy
do j=k,2,-1
sum=sum+((e(i,j)-e(i,j-1))/dy)*((a/(a-1.0d0))**(k-j+1))
enddo
dby(i,1)=sum/((1.0d0-a)+(((a/(a-1.0d0))**(k-1))*a))
do j=2,k
dby(i,j)=(((e(i,j)-e(i,j-1))/dy)-(a*dby(i,j-1)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yinvector_to_yangvector(phin,landan,phie,landae,vphin,vlandan,vphie,vlandae)
implicit none
real*8::phin,landan,phie,landae,pi,phin1,phie1,vphin,vlandan,vphie,vlandae,cossaay,sinsaay
pi=dasin(1.0d0)*2.0d0
phin1=phin+(pi/2.0d0)
phie1=phie+(pi/2.0d0)
cossaay=-dsin(landae)*dsin(landan)
if(dabs(dsin(phin1))<1.0d-8)then
sinsaay=dcos(landan)/dsin(phie1)
else
sinsaay=-dcos(landae)/dsin(phin1)
endif
vphie=(vphin*cossaay)-(vlandan*sinsaay)
vlandae=(vphin*sinsaay)+(vlandan*cossaay)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yangvector_to_yinvector(phin,landan,phie,landae,vphin,vlandan,vphie,vlandae)
implicit none
real*8::phin,landan,phie,landae,pi,phin1,phie1,vphin,vlandan,vphie,vlandae,cossaay,sinsaay
pi=dasin(1.0d0)*2.0d0
phin1=phin+(pi/2.0d0)
phie1=phie+(pi/2.0d0)
cossaay=-dsin(landae)*dsin(landan)
if(dabs(dsin(phin1))<1.0d-8)then
sinsaay=dcos(landan)/dsin(phie1)
else
sinsaay=-dcos(landae)/dsin(phin1)
endif
vphin=(vphie*cossaay)+(vlandae*sinsaay)
vlandan=-(vphie*sinsaay)+(vlandae*cossaay)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine num(phi3,landa3,k1,k2,dlanda,dphi,n_epsilon)
implicit none
integer::k1,k2,n_epsilon
real*8::pi,landa3,phi3,dphi,dlanda,top_x_limit,top_y_limit
pi=2.d0*dasin(1.0d0)
top_x_limit=pi
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)
k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine cubic1(h1,h2,h3,h4,x1,x2,x3,x4,xm,hm)
implicit none
real*8::h1,h2,h3,h4,x1,x2,x3,x4,xm,hm
real*8::a,b,c,d,c1,c2,c3,c4,c5
c1=(x1+x2+x3)
c2=(x2*x2)+(x1*x1)+(x2*x1)
c3=((h3-h1))/((x3-x1)*(x3-x2))
c4=((h2-h1))/((x2-x1)*(x3-x2))
c5=(h2-h1)/(x2-x1)
a=(h4-(x4*x4*(c3-c4))-(x4*c5)+(x4*(x1+x2)*(c3-c4))-h1+(x1*x1*(c3-c4))+(x1*c5)-(x1*(x1+x2)*(c3-c4)))/&
((x4*x4*x4)-(x4*x4*c1)-(x4*c2)+(x4*(x1+x2)*c1)-(x1*x1*x1)+(x1*x1*c1)+(x1*c2)-(x1*(x1+x2)*c1))
b=c3-c4-(c1*a)
c=c5-(c2*a)-(b*(x1+x2))
d=h1-(a*x1*x1*x1)-(b*x1*x1)-(c*x1)
hm=(a*xm*xm*xm)+(b*xm*xm)+(c*xm)+d
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bilinear_glob1(phi3,landa3,heb,nx,ny,hem,dlanda,dphi)
implicit none
integer::nx,ny
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hem,ddx,ddy,dddx,dddy,dd1,dd2,dd3,dd4,aa,phi,landa
real*8::heb(nx,ny)

aa=6.37122d6
pi=2.d0*dasin(1.0d0)


k1=int(((landa3+pi+1.0d-10)/dlanda))+1
k2=int(((phi3+(pi/2.0d0)+1.0d-10)/dphi))+1

landa=dble(k1-1)*dlanda-pi
phi=(dble(k2-1))*dphi-(pi/2.0d0)

call distance(phi3,phi3,landa,landa+dlanda,ddx)
call distance(phi,phi+dphi,landa3,landa3,ddy)
call distance(phi3,phi3,landa3,landa,dddx)
call distance(phi3,phi,landa3,landa3,dddy)
if(dabs(ddx)<=1.0d-8)then
hem=heb(k1,k2)
else
dd1=(ddx-dddx)*(ddy-dddy)/(ddx*ddy)
dd2=dddx*(ddy-dddy)/(ddx*ddy)
dd3=(ddx-dddx)*dddy/(ddx*ddy)
dd4=dddx*dddy/(ddx*ddy)
 

hem=(dd1*heb(k1,k2)+dd2*heb(k1+1,k2)+dd3*heb(k1,k2+1)+dd4*heb(k1+1,k2+1))
endif
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine distance(phi1,phi2,landa1,landa2,rr)
implicit none
real*8::phi1,phi2,landa1,landa2,rr,pi,aa,arg
aa=6.37122d6
pi=2.d0*dasin(1.0d0)
arg=dsin(phi1)*dsin(phi2)+dcos(phi1)*dcos(phi2)*dcos(landa1-landa2)
if(arg>=1.0d0)then
arg=1.0d0
elseif(arg<=-1.0d0)then
arg=-1.0d0
endif
rr=aa*dacos(arg)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dcxp(e,dcx,nx,ny,dx)
implicit none
integer::i,j,nx,ny
real*8::e(nx,ny),dcx(nx,ny),dx
do j=1,ny
do i=1,nx
if(i==1)then
dcx(i,j)=(e(i+1,j)-e(nx,j))/(2.0d0*dx)
elseif(i==nx)then
dcx(i,j)=(e(1,j)-e(i-1,j))/(2.0d0*dx)
else
dcx(i,j)=(e(i+1,j)-e(i-1,j))/(2.0d0*dx)
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dcyp(e,dcy,nx,ny,dy)
implicit none
integer::i,j,nx,ny
real*8::e(nx,ny),dcy(nx,ny),dy
do j=1,ny
do i=1,nx
if(j==1)then
dcy(i,j)=(e(i,j+1)-e(i,ny))/(2.0d0*dy)
elseif(j==ny)then
dcy(i,j)=(e(i,1)-e(i,j-1))/(2.0d0*dy)
else
dcy(i,j)=(e(i,j+1)-e(i,j-1))/(2.0d0*dy)
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  RK_init(he,dcx,dcy,nx,ny,cx,cy,a,h,hh,phi,delt)
implicit none
integer::i,j,nx,ny
real*8::he(nx,ny),dcx(nx,ny),dcy(nx,ny),cx(nx,ny),cy(nx,ny),h(nx,ny)&
,hh(nx,ny),phi(ny),a,delt

do i=1,nx
do j=1,ny


call Advection_sphere(h(i,j),dcx(i,j),dcy(i,j),cx(i,j),cy(i,j),phi(j),delt)




hh(i,j)=a*h(i,j)+he(i,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  RK4_final(he1,he2,h1,h2,h3,h4,nx,ny)
implicit none
integer::i,j,nx,ny
real*8::he1(nx,ny),he2(nx,ny),h1(nx,ny),h2(nx,ny),h3(nx,ny),h4(nx,ny)&
,b1,b2,b3,b4,aa
aa=6.37122d6
b1=1.0d0/6.0d0;b2=1.0d0/3.0d0;b3=1.0d0/3.0d0;b4=1.0d0/6.0d0
do i=1,nx
do j=1,ny
he2(i,j)=he1(i,j)+b1*h1(i,j)+b2*h2(i,j)+b3*h3(i,j)+b4*h4(i,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Advection_sphere(h,dcx,dcy,cx,cy,phi,delt)
implicit none
real*8::ux,uy,aa,dcx,dcy,cx,cy,phi,h,delt

aa=6.37122d6
if(dabs(dcos(phi))<=1.0d-3)then
ux=0.0d0
else
ux=dcx*cx/(aa*dcos(phi))
endif
uy=dcy*cy/aa
h=-delt*(ux+uy)

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  initial_conditions(nx,ny,uphiin,hcapitin,hprimein,h00)
implicit none
integer::nx,ny,i,j
real*8::aa,omega,pi,twopi,hhat,g,umax,phi0,phi1,en,phi2,alpha,beta,h00,h_mean(ny),&
dlandain,dphiin,landain(nx),phiin(ny),foin(ny),uphiin(ny),hcapitin(ny),hprimein(nx,ny),&
sum

aa=6.37122d6;omega=7.292d-5
pi=2.0d0*dasin(1.0d0);twopi=2.0d0*pi;hhat=120.0d0
g=9.80616d0;umax=80.0d0;phi0=pi/7.0d0
phi1=pi/2.0d0-phi0;en=dexp(-4/((phi1-phi0)*(phi1-phi0)))
phi2=pi/4.0d0;alpha=1.0d0/3.0d0;beta=1.0d0/15.0d0
dlandain=twopi/dble(nx);dphiin=pi/dble((ny-1))
do i=1,nx
landain(i)=dble(i-1)*dlandain-pi
enddo

do j=1,ny
phiin(j)=dble(j-1)*dphiin-(pi/2.0d0)
foin(j)=2.0d0*omega*dsin(phiin(j))
enddo
do j=1,ny
if(phiin(j)>phi0.and.phiin(j)<phi1)then
uphiin(j)=(umax/en)*dexp(1.0d0/((phiin(j)-phi0)*(phiin(j)-phi1)))
else
uphiin(j)=0.0d0
endif

enddo
do i=1,ny
sum=0.0d0
do j=1,i
sum=sum+uphiin(j)*(foin(j)+(dtan(phiin(j))*uphiin(j)/aa))*dphiin
enddo
h_mean(i)=aa*sum/g

enddo

sum=0.0d0
do i=1,ny
sum=sum+h_mean(i)*dcos(phiin(i))*dphiin
enddo
h00=(sum/2.0d0)+1.0d4
print*,'global h00= ',h00
do i=1,ny
hcapitin(i)=h00-h_mean(i)

enddo

do j=1,ny
do i=1,nx
hprimein(i,j)=hhat*dcos(phiin(j))*dexp(-(landain(i)/alpha)*(landain(i)/alpha))*&
dexp(-(((phi2-phiin(j))/beta)*((phi2-phiin(j))/beta)))

enddo
enddo

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  swe_sphere(hu,hv,hh,u,v,h,dux,duy,dvx,dvy,dhx,dhy,phi,f,delt)
implicit none
real*8::hu,hv,hh,u,v,h,dux,duy,dvx,dvy,dhx,dhy,aa,f,g,phi,delt

g=9.80616d0
aa=6.37122d6

hu=-delt*((u*dux/(aa*dcos(phi)))+(v*duy/(aa))-(f*v)-(u*v*dtan(phi)/aa)+(g*dhx/(aa*dcos(phi))))
hv=-delt*((u*dvx/(aa*dcos(phi)))+(v*dvy/(aa))+(f*u)+(u*u*dtan(phi)/aa)+(g*dhy/(aa)))
hh=-delt*((u*dhx/(aa*dcos(phi)))+(v*dhy/(aa))+(h*dux/(aa*dcos(phi)))+(h*dvy/aa)-(h*v*dsin(phi)/(aa*dcos(phi))))

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  RK_init_swe(u,v,h,dux,duy,dvx,dvy,dhx,dhy,nx,ny,a,hu,hv,hh,hhu,hhv,hhh,phi,f,delt,hhu0,hhv0,hhh0)
implicit none
integer::i,j,nx,ny
real*8::u(nx,ny),v(nx,ny),h(nx,ny),dux(nx,ny),duy(nx,ny),dvx(nx,ny),dvy(nx,ny)&
,dhx(nx,ny),dhy(nx,ny),hu(nx,ny),hv(nx,ny),hh(nx,ny),hhu(nx,ny),hhv(nx,ny),hhh(nx,ny)&
,phi(ny),a,delt,f(nx,ny),a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a12,a13,a14,a15,hhu0(nx,ny)&
,hhv0(nx,ny),hhh0(nx,ny),a16,a17,a18

do j=1,ny
do i=1,nx


a4=u(i,j)
a5=v(i,j)
a6=h(i,j)
a7=dux(i,j)
a8=duy(i,j)
a9=dvx(i,j)
a10=dvy(i,j)

a12=dhx(i,j)
a13=dhy(i,j)
a14=phi(j)
a15=f(i,j)
a16=hhu0(i,j)
a17=hhv0(i,j)
a18=hhh0(i,j)
call swe_sphere(a1,a2,a3,a16,a17,a18,a7,a8,a9,a10,a12,a13,a14,a15,delt)



hhu(i,j)=a*a1+a4
hhv(i,j)=a*a2+a5
hhh(i,j)=a*a3+a6
hu(i,j)=a1
hv(i,j)=a2
hh(i,j)=a3

enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  zita_delta(u,v,dux,dvx,duy,dvy,phi,nx,ny,zita,delta)
integer::i,j,nx,ny
real*8::dux(nx,ny),dvx(nx,ny),duy(nx,ny),dvy(nx,ny),zita(nx,ny),delta(nx,ny)&
,phi(ny),aa,u(nx,ny),v(nx,ny)
aa=6.37122d6
do i=1,nx
do j=1,ny
if(dabs(dcos(phi(j)))<=1.0d-3)then
zita(i,j)=0.0d0
delta(i,j)=0.0d0
else
zita(i,j)=(dvx(i,j)-duy(i,j)*dcos(phi(j))+u(i,j)*dsin(phi(j)))/(aa*dcos(phi(j)))
delta(i,j)=(dux(i,j)+dvy(i,j)*dcos(phi(j))-v(i,j)*dsin(phi(j)))/(aa*dcos(phi(j)))
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  initial_yang_total(nx,ny,uphie,hcapite,hprimee,h00,ue,ve)
implicit none
integer::nx,ny,i,j,k1,k2

real*8::aa,omega,pi,twopi,hhat,g,umax,phi0,phi1,en,phi2,alpha,beta,h00,&
dlandae,dphie,landae(nx),phie(ny),foe(nx,ny),uphie(nx,ny),hcapite(nx,ny),hprimee(nx,ny)&
,sum1,sum2,sum,ue(nx,ny),ve(nx,ny),a10,a11,a12,a13,a14,a15,a16

aa=6.37122d6;omega=7.292d-5
pi=2.0d0*dasin(1.0d0);twopi=2.0d0*pi;hhat=120.0d0
g=9.80616d0;umax=80.0d0;phi0=pi/7.0d0
phi1=pi/2.0d0-phi0;en=dexp(-4/((phi1-phi0)*(phi1-phi0)))
phi2=pi/4.0d0;alpha=1.0d0/3.0d0;beta=1.0d0/15.0d0

dlandae=twopi/dble(nx);dphie=pi/dble((ny-1))

do i=1,nx
landae(i)=-pi+dble(i-1)*dlandae
enddo

do j=1,ny
phie(j)=dble(j-1)*dphie-(pi/2.0d0)

do i=1,nx
a12=phie(j)
a13=landae(i)
foe(i,j)=-2.0d0*omega*dcos(a12)*dsin(a13)
enddo
enddo
print*,'Yang h00= ',h00



do j=1,ny
a12=phie(j)
do i=1,nx
a13=landae(i)
call yang_to_yin(a10,a11,a12,a13)
if(a10>phi0.and.a10<phi1)then
uphie(i,j)=(umax/en)*dexp(1.0d0/((a10-phi0)*(a10-phi1)))
else
uphie(i,j)=0.0d0
endif
a14=uphie(i,j)
call yinvector_to_yangvector(a10,a11,a12,a13,0.0d0,a14,a15,a16)
ve(i,j)=a15
ue(i,j)=a16

enddo
enddo
do j=1,ny
do i=1,nx
sum1=0.0d0
do k1=1,j
a12=phie(k1)
a13=landae(i)
call yang_to_yin(a10,a11,a12,a13)
if(dabs(dsin(a10+(pi/2.0d0)))<=1.0d-6)then
sum1=sum1
else
sum1=sum1+uphie(i,k1)*(foe(i,k1)+(dtan(a10)*uphie(i,k1)/aa))*(dcos(a12+(pi/2.0d0))*dsin(a13)/&
dsin(a10+(pi/2.0d0)))*dphie
endif
enddo
sum2=0.0d0
do k2=1,i

a12=phie(j)
a13=landae(k2)
call yang_to_yin(a10,a11,a12,a13)


if(dabs(dsin(a10+(pi/2.0d0)))<=1.0d-6)then
sum2=sum2
else
sum2=sum2+uphie(k2,j)*(foe(k2,j)+(dtan(a10)*uphie(k2,j)/aa))*(dsin(a12+(pi/2.0d0))*dcos(a13)/&
dsin(a10+(pi/2.0d0)))*dlandae
endif
enddo

sum=(-sum1-sum2)/2.0d0
hcapite(i,j)=h00-(aa*sum/g)

enddo
enddo
do j=1,ny
do i=1,nx
a12=phie(j)
a13=landae(i)
call yang_to_yin(a10,a11,a12,a13)
hprimee(i,j)=hhat*dcos(a10)*dexp(-(a11/alpha)*(a11/alpha))*&
dexp(-(((phi2-a10)/beta)*((phi2-a10)/beta)))

enddo
enddo

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dec2x(func,dcx,nx,ny,dlanda)
implicit none
integer::nx,ny
real*8::func(nx,ny),dfx(nx,ny),dlanda,dcx(nx,ny)
call dfx42p(func,dfx,nx,ny,dlanda)
call dbx42p(dfx,dcx,nx,ny,dlanda)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dec2y(func,dcy,nx,ny,dphi)
implicit none
integer::nx,ny
real*8::func(nx,ny),dfy(nx,ny),dphi,dcy(nx,ny)
call dfy42p(func,dfy,nx,ny,dphi)
call dby42p(dfy,dcy,nx,ny,dphi)
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine laplasian(func,nx,ny,dphi,dlanda,laplas,phi)
implicit none
integer::nx,ny,i,j
real*8::func(nx,ny),dfx(nx,ny),dphi,dcx(nx,ny),laplas(nx,ny),aa,phi(ny),dby(nx,ny)&
,dcy(nx,ny),funcy(nx,ny),dlanda
aa=6.37122d6
call dfx42p(func,dfx,nx,ny,dlanda)
call dbx42p(dfx,dcx,nx,ny,dlanda)
call dby42p(func,dby,nx,ny,dphi)
do j=1,ny
do i=1,nx
funcy(i,j)=dby(i,j)*dcos(phi(j))
enddo
enddo
call dfy42p(funcy,dcy,nx,ny,dphi)
do j=1,ny
do i=1,nx
if(dabs(dcos(phi(j)))<=1.0d-3)then
laplas(i,j)=0.0d0
else
laplas(i,j)=(dcx(i,j)/(aa*aa*dcos(phi(j))*dcos(phi(j))))+(dcy(i,j)/(aa*aa*dcos(phi(j))))
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine laplasian_longitude(func,nx,ny,dlanda,laplas_x,phi)
implicit none
integer::nx,ny,i,j
real*8::func(nx,ny),dfx(nx,ny),dlanda,dcx(nx,ny),laplas_x(nx,ny),aa,phi(ny)
aa=6.37122d6
call dfx42p(func,dfx,nx,ny,dlanda)
call dbx42p(dfx,dcx,nx,ny,dlanda)
do j=1,ny
do i=1,nx
if(dabs(dcos(phi(j)))<=1.0d-3)then
laplas_x(i,j)=0.0d0
else
laplas_x(i,j)=dcx(i,j)/(aa*aa*dcos(phi(j))*dcos(phi(j)))
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine laplasian_latitude(func,nx,ny,dphi,laplas_y,phi)
implicit none
integer::nx,ny,i,j
real*8::func(nx,ny),dphi,laplas_y(nx,ny),aa,phi(ny),dby(nx,ny),dcy(nx,ny),funcy(nx,ny)
aa=6.37122d6

call dby42p(func,dby,nx,ny,dphi)
do j=1,ny
do i=1,nx
funcy(i,j)=dby(i,j)*dcos(phi(j))
enddo
enddo
call dfy42p(funcy,dcy,nx,ny,dphi)
do j=1,ny
do i=1,nx
if(dabs(dcos(phi(j)))<=1.0d-3)then
laplas_y(i,j)=0.0d0
else
laplas_y(i,j)=dcy(i,j)/(aa*aa*dcos(phi(j)))
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  decx(func,dcx,nx,ny,dlanda)
implicit none
integer::nx,ny,i,j
real*8::func(nx,ny),dfx(nx,ny),dbx(nx,ny),dlanda,dcx(nx,ny)
call dfx42p(func,dfx,nx,ny,dlanda)
call dbx42p(func,dbx,nx,ny,dlanda)
do j=1,ny
do i=1,nx
dcx(i,j)=(dfx(i,j)+dbx(i,j))/2.0d0
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  decy(func,dcy,nx,ny,dphi)
implicit none
integer::nx,ny,i,j
real*8::func(nx,ny),dfy(nx,ny),dby(nx,ny),dphi,dcy(nx,ny)
call dfy42p(func,dfy,nx,ny,dphi)
call dby42p(func,dby,nx,ny,dphi)
do j=1,ny
do i=1,nx
dcy(i,j)=(dfy(i,j)+dby(i,j))/2.0d0
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yin_to_yang(phin,landan,phie,landae)
real*8::phin,landan,phie,landae,pi
pi=dasin(1.0d0)*2.0d0
if(dabs(phin+(pi/2.0d0))<=1.0d-8)then
phie=0.0d0
landae=pi/2.0d0
elseif(dabs(phin-(pi/2.0d0))<=1.0d-8)then
phie=0.0d0
landae=-pi/2.0d0
elseif(dabs(phin)<=1.0d-8.and.dabs(landan-(pi/2.0d0))<=1.0d-8)then
phie=-pi/2.0d0
landae=0.0d0
elseif(dabs(phin)<=1.0d-8.and.dabs(landan+(pi/2.0d0))<=1.0d-8)then
phie=pi/2.0d0
landae=0.0d0
elseif(dabs(phin)<=1.0d-8.and.dabs(landan)<=1.0d-8)then
phie=0.0d0
landae=-pi
elseif(dabs(phin)<=1.0d-8.and.dabs(landan+pi)<=1.0d-8)then
phie=0.0d0
landae=0.0d0
elseif(dabs(phin)<=1.0d-8.and.landan>-pi.and.landan<-(pi/2.0d0))then
phie=landan+pi
landae=0.0d0
elseif(dabs(phin)<=1.0d-8.and.landan>-(pi/2.0d0).and.landan<0.0d0)then
phie=-landan
landae=-pi
elseif(dabs(phin)<=1.0d-8.and.landan>0.0d0.and.landan<(pi/2.0d0))then
phie=-landan
landae=-pi
elseif(dabs(phin)<=1.0d-8.and.landan>(pi/2.0d0))then
phie=landan-pi
landae=0.0d0

elseif(dabs(landan)<=1.0d-8.and.phin>0.0d0)then
phie=0.0d0
landae=phin-pi
elseif(dabs(landan)<=1.0d-8.and.phin<0.0d0)then
phie=0.0d0
landae=phin+pi
elseif(dabs(landan+pi)<=1.0d-8.and.dabs(phin)<=1.0d-8)then
phie=0.0d0
landae=0.0d0
elseif(dabs(landan+pi)<=1.0d-8.and.phin>0.0d0)then
phie=0.0d0
landae=-phin
elseif(dabs(landan+pi)<=1.0d-8.and.phin<0.0d0)then
phie=0.0d0
landae=-phin
elseif(phin>0.0d0.and.dabs(landan+pi/2.0d0)<=1.0d-8)then
phie=dasin(dcos(phin))
landae=-pi/2.0d0
elseif(phin>0.0d0.and.dabs(landan-pi/2.0d0)<=1.0d-8)then
phie=-dasin(dcos(phin))
landae=-pi/2.0d0
elseif(phin<0.0d0.and.dabs(landan+pi/2.0d0)<=1.0d-8)then
phie=dasin(dcos(phin))
landae=pi/2.0d0
elseif(phin<0.0d0.and.dabs(landan-pi/2.0d0)<=1.0d-8)then
phie=-dasin(dcos(phin))
landae=pi/2.0d0
elseif(phin>0.0d0.and.landan>-pi.and.landan<-(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))
elseif(phin>0.0d0.and.landan>-(pi/2.0d0).and.landan<0.0d0)then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))-pi
elseif(phin>0.0d0.and.landan>0.0d0.and.landan<(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))-pi
elseif(phin>0.0d0.and.landan>(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))

elseif(phin<0.0d0.and.landan>-pi.and.landan<-(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))
elseif(phin<0.0d0.and.landan>-(pi/2.0d0).and.landan<0.0d0)then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))+pi
elseif(phin<0.0d0.and.landan>0.0d0.and.landan<(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))+pi
elseif(phin<0.0d0.and.landan>(pi/2.0d0))then
phie=-dasin(dcos(phin)*dsin(landan))
landae=datan(dtan(phin)/dcos(landan))
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  yang_to_yin(phin,landan,phie,landae)
real*8::phin,landan,phie,landae,pi
pi=dasin(1.0d0)*2.0d0
if(dabs(phie+(pi/2.0d0))<=1.0d-8)then
phin=0.0d0
landan=pi/2.0d0
elseif(dabs(phie-(pi/2.0d0))<=1.0d-8)then
phin=0.0d0
landan=-pi/2.0d0
elseif(dabs(phie)<=1.0d-8.and.dabs(landae-(pi/2.0d0))<=1.0d-8)then
phin=-pi/2.0d0
landan=0.0d0
elseif(dabs(phie)<=1.0d-8.and.dabs(landae+(pi/2.0d0))<=1.0d-8)then
phin=pi/2.0d0
landan=0.0d0
elseif(dabs(phie)<=1.0d-8.and.dabs(landae)<=1.0d-8)then
phin=0.0d0
landan=-pi
elseif(dabs(phie)<=1.0d-8.and.dabs(landae+pi)<=1.0d-8)then
phin=0.0d0
landan=0.0d0
elseif(dabs(phie)<=1.0d-8.and.landae>-pi.and.landae<-(pi/2.0d0))then
phin=landae+pi
landan=0.0d0
elseif(dabs(phie)<=1.0d-8.and.landae>-(pi/2.0d0).and.landae<0.0d0)then
phin=-landae
landan=-pi
elseif(dabs(phie)<=1.0d-8.and.landae>0.0d0.and.landae<(pi/2.0d0))then
phin=-landae
landan=-pi
elseif(dabs(phie)<=1.0d-8.and.landae>(pi/2.0d0))then
phin=landae-pi
landan=0.0d0
elseif(dabs(landae)<=1.0d-8.and.phie>0.0d0)then
phin=0.0d0
landan=phie-pi
elseif(dabs(landae)<=1.0d-8.and.phie<0.0d0)then
phin=0.0d0
landan=phie+pi
elseif(dabs(landae+pi)<=1.0d-8.and.dabs(phie)<=1.0d-8)then
phin=0.0d0
landan=0.0d0
elseif(dabs(landae+pi)<=1.0d-8.and.phie>0.0d0)then
phin=0.0d0
landan=-phie
elseif(dabs(landae+pi)<=1.0d-8.and.phie<0.0d0)then
phin=0.0d0
landan=-phie
elseif(phie>0.0d0.and.dabs(landae+pi/2.0d0)<=1.0d-8)then
phin=dasin(dcos(phie))
landan=-pi/2.0d0
elseif(phie>0.0d0.and.dabs(landae-pi/2.0d0)<=1.0d-8)then
phin=-dasin(dcos(phie))
landan=-pi/2.0d0
elseif(phie<0.0d0.and.dabs(landae+pi/2.0d0)<=1.0d-8)then
phin=dasin(dcos(phie))
landan=pi/2.0d0
elseif(phie<0.0d0.and.dabs(landae-pi/2.0d0)<=1.0d-8)then
phin=-dasin(dcos(phie))
landan=pi/2.0d0

elseif(phie>0.0d0.and.landae>-pi.and.landae<-(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))

elseif(phie>0.0d0.and.landae>-(pi/2.0d0).and.landae<0.0d0)then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))-pi

elseif(phie>0.0d0.and.landae>0.0d0.and.landae<(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))-pi

elseif(phie>0.0d0.and.landae>(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))

elseif(phie<0.0d0.and.landae>-pi.and.landae<-(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))

elseif(phie<0.0d0.and.landae>-(pi/2.0d0).and.landae<0.0d0)then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))+pi

elseif(phie<0.0d0.and.landae>0.0d0.and.landae<(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))+pi

elseif(phie<0.0d0.and.landae>(pi/2.0d0))then
phin=-dasin(dcos(phie)*dsin(landae))
landan=datan(dtan(phie)/dcos(landae))
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  merging_Yin_Yang(Fn,Fe,Fne,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
implicit none
integer::nx,ny,mx,my,knex,kney,i,j,n_epsilon_x,n_epsilon_y
real*8::pi,top_x_limit,top_y_limit,dphi,dlanda,Fn(mx,my),Fe(mx,my),Fne(nx,ny)&
,a10,a11,a12,a13,hem,phi(ny),landa(nx)
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+dble(n_epsilon_x)*dlanda
top_y_limit=(pi/4.0d0)+dble(n_epsilon_y)*dphi

do j=1,ny
phi(j)=-(pi/2.0d0)+dble(j-1)*dphi
do i=1,nx
landa(i)=-pi+dble(i-1)*dlanda
if(phi(j)>=-top_y_limit.and.phi(j)<=top_y_limit.and.landa(i)>=-top_x_limit.and.landa(i)<=top_x_limit)then
Fne(i,j)=Fn(i-knex,j-kney)
else
a10=phi(j)
a11=landa(i)
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,Fe,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
Fne(i,j)=hem
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Filter_boundaries(Un,Ue,Vn,Ve,Hen,Hee,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
implicit none
integer::my,mx,i,j,n_epsilon_x,n_epsilon_y
real*8::pi,top_x_limit,top_y_limit,dphi,dlanda,Un(mx,my),Ue(mx,my),Vn(mx,my)&
,Ve(mx,my),Hen(mx,my),Hee(mx,my),a10,a11,a12,a13,hem,um,vm,ux,vy
character(5)::int_method
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)
do i=1,mx
do j=1,3
a10=-top_y_limit+dble(j-1)*dphi
a11=-top_x_limit+dble(i-1)*dlanda
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
call select_int_method(a12,a13,Ue,mx,my,um,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
call select_int_method(a12,a13,Ve,mx,my,vm,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
Hen(i,j)=hem
call yangvector_to_yinvector(a10,a11,a12,a13,vy,ux,vm,um)
Un(i,j)=ux
Vn(i,j)=vy
enddo
do j=my-2,my
a10=-top_y_limit+dble(j-1)*dphi
a11=-top_x_limit+dble(i-1)*dlanda
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
call select_int_method(a12,a13,Ue,mx,my,um,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
call select_int_method(a12,a13,Ve,mx,my,vm,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
Hen(i,j)=hem
call yangvector_to_yinvector(a10,a11,a12,a13,vy,ux,vm,um)
Un(i,j)=ux
Vn(i,j)=vy
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  merging_vector_Yin_Yang(hen,hee,hene,un,ue,une,vn,ve,vne,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
implicit none
integer::nx,ny,mx,my,knex,kney,i,j,n_epsilon_x,n_epsilon_y
real*8::pi,top_x_limit,top_y_limit,dphi,dlanda,un(mx,my),ue(mx,my),une(nx,ny)&
,vn(mx,my),ve(mx,my),vne(nx,ny),a10,a11,a12,a13,hm,um,vm,ux,vy,phi(ny),landa(nx)&
,hen(mx,my),hee(mx,my),hene(nx,ny)
character(5)::int_method
knex=(nx-mx)/2
kney=(ny-my)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+dble(n_epsilon_x)*dlanda
top_y_limit=(pi/4.0d0)+dble(n_epsilon_y)*dphi

do j=1,ny
phi(j)=-(pi/2.0d0)+dble(j-1)*dphi
do i=1,nx
landa(i)=-pi+dble(i-1)*dlanda
if(phi(j)>=-top_y_limit.and.phi(j)<=top_y_limit.and.landa(i)>=-top_x_limit.and.landa(i)<=top_x_limit)then
hene(i,j)=hen(i-knex,j-kney)
une(i,j)=un(i-knex,j-kney)
vne(i,j)=vn(i-knex,j-kney)
else
a10=phi(j)
a11=landa(i)
call yin_to_yang(a10,a11,a12,a13)
call select_int_method(a12,a13,hee,mx,my,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
hene(i,j)=hm
call select_int_method(a12,a13,ue,mx,my,um,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
call select_int_method(a12,a13,ve,mx,my,vm,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
call yangvector_to_yinvector(a10,a11,a12,a13,vy,ux,vm,um)
une(i,j)=ux
vne(i,j)=vy
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  provide_Yin_scalar_to_derive(hen,hee,dlanda,dphi,nx,ny,mx,my,n_epsilon,he_zone,he_meridion,int_method)
implicit none
integer::i,j,nx,ny,mx,my,n_epsilon,k1,k2,knex,kney
real*8::hen(mx,my),hee(mx,my),he_zone(nx,my),he_meridion(mx,2*(ny-1)),a10,a11,a12,a13,hem,pi,dlanda,dphi&
,top_x_limit,top_y_limit
character(5)::int_method
kney=int((real(ny-my)+0.01)/2.0)
knex=int((real(nx-mx+1)+0.01)/2.0)
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)
do j=1,my
a10=-top_y_limit+dble(j-1)*dphi
do i=1,nx
a11=-pi+dble(i-1)*dlanda
if(a11>=-top_x_limit.and.a11<=top_x_limit)then
he_zone(i,j)=hen(i-knex,j)
else

call yin_to_yang(a10,a11,a12,a13)
call interpolation_scalar(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon,int_method)
he_zone(i,j)=hem
endif
enddo
enddo
do i=1,mx
a11=-top_x_limit+dble(i-1)*dlanda
do j=1,ny
a10=-(pi/2.0d0)+dble(j-1)*dphi
if(a10>=-top_y_limit.and.a10<=top_y_limit)then
he_meridion(i,j)=hen(i,j-kney)
else

call yin_to_yang(a10,a11,a12,a13)
call interpolation_scalar(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon,int_method)
he_meridion(i,j)=hem
endif
enddo
do j=ny+1,2*(ny-1)
a10=(pi/2.0d0)-dble(j-ny)*dphi
if(a11<0.0d0)then
a11=a11+pi
elseif(a11>=0.0d0)then
a11=a11-pi
endif
if(a10>=-top_y_limit.and.a10<=top_y_limit.and.a11>=-top_x_limit.and.a11<=top_x_limit)then
call num(a10,a11,k1,k2,dlanda,dphi,n_epsilon)
he_meridion(i,j)=hen(k1,k2)
else
call yin_to_yang(a10,a11,a12,a13)
call interpolation_scalar(a12,a13,hee,mx,my,hem,dlanda,dphi,n_epsilon,int_method)
he_meridion(i,j)=hem
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  provide_Yang_scalar_to_derive(hen,hee,dlanda,dphi,nx,ny,mx,my,n_epsilon,he_zone,he_meridion,int_method)
implicit none
integer::i,j,nx,ny,mx,my,n_epsilon,k1,k2,knex,kney
real*8::hen(mx,my),hee(mx,my),he_zone(nx,my),he_meridion(mx,2*(ny-1)),a10,a11,a12,a13,hem,pi,dlanda,dphi&
,top_x_limit,top_y_limit
character(5)::int_method
kney=int((real(ny-my)+0.01)/2.0)
knex=int((real(nx-mx+1)+0.01)/2.0)
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)
do j=1,my
a12=-top_y_limit+dble(j-1)*dphi
do i=1,nx
a13=-pi+dble(i-1)*dlanda
if(a13>=-top_x_limit.and.a13<=top_x_limit)then
he_zone(i,j)=hee(i-knex,j)
else

call yang_to_yin(a10,a11,a12,a13)
call interpolation_scalar(a10,a11,hen,mx,my,hem,dlanda,dphi,n_epsilon,int_method)

he_zone(i,j)=hem
endif
enddo
enddo


do i=1,mx
a13=-top_x_limit+dble(i-1)*dlanda
do j=1,ny
a12=-(pi/2.0d0)+dble(j-1)*dphi
if(a12>=-top_y_limit.and.a12<=top_y_limit)then
he_meridion(i,j)=hee(i,j-kney)
else

call yang_to_yin(a10,a11,a12,a13)
call interpolation_scalar(a10,a11,hen,mx,my,hem,dlanda,dphi,n_epsilon,int_method)

he_meridion(i,j)=hem
endif
enddo
do j=ny+1,2*(ny-1)
a12=(pi/2.0d0)-dble(j-ny)*dphi
if(a13<0.0d0)then
a13=a13+pi
elseif(a13>=0.0d0)then
a13=a13-pi
endif
if(a12>=-top_y_limit.and.a12<=top_y_limit.and.a13>=-top_x_limit.and.a13<=top_x_limit)then
call num(a12,a13,k1,k2,dlanda,dphi,n_epsilon)
he_meridion(i,j)=hee(k1,k2)
else
call yang_to_yin(a10,a11,a12,a13)
call interpolation_scalar(a10,a11,hen,mx,my,hem,dlanda,dphi,n_epsilon,int_method)
he_meridion(i,j)=hem
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  provide_Yin_vector_to_derive(hen,hee,vxn,vyn,vxe,vye,dlanda,dphi,nx,ny,my,n_epsilon&
,he_meridion,vx_meridion,vxcosy_meridion,vy_meridion,vycosy_meridion,int_method)
implicit none
integer::i,j,nx,ny,my,n_epsilon,k1,k2,kney
real*8::hen(nx,my),hee(nx,my),vxn(nx,my),vyn(nx,my),vxe(nx,my),vye(nx,my),vx_meridion(nx,2*(ny-1))&
,vy_meridion(nx,2*(ny-1)),a10,a11,a12,a13,hem,vxm,vym,vx,vy,pi,dlanda,dphi&
,top_x_limit,top_y_limit,vxcosy_meridion(nx,2*(ny-1)),vycosy_meridion(nx,2*(ny-1))&
,he_meridion(nx,2*(ny-1))
character(5)::int_method

kney=int((real(ny-my)+0.01)/2.0)

pi=2.0d0*dasin(1.0d0)
top_x_limit=pi
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)

do i=1,nx
a11=-top_x_limit+dble(i-1)*dlanda
do j=1,ny
a10=-(pi/2.0d0)+dble(j-1)*dphi
if(a10>=-top_y_limit.and.a10<=top_y_limit)then
he_meridion(i,j)=hen(i,j-kney)
vx_meridion(i,j)=vxn(i,j-kney)
vy_meridion(i,j)=vyn(i,j-kney)
vxcosy_meridion(i,j)=vxn(i,j-kney)*dcos(a10)
vycosy_meridion(i,j)=vyn(i,j-kney)*dcos(a10)
else

call yin_to_yang(a10,a11,a12,a13)
call new_interpolation(a12,a13,hee,vxe,vye,nx,my,hem,vxm,vym,dlanda,dphi,n_epsilon,int_method)
call yangvector_to_yinvector(a10,a11,a12,a13,vy,vx,vym,vxm)
he_meridion(i,j)=hem
vx_meridion(i,j)=vx
vy_meridion(i,j)=vy
vxcosy_meridion(i,j)=vx*dcos(a10)
vycosy_meridion(i,j)=vy*dcos(a10)
endif
enddo
do j=ny+1,2*(ny-1)
a10=(pi/2.0d0)-dble(j-ny)*dphi
if(a11<0.0d0)then
a11=a11+pi
elseif(a11>=0.0d0)then
a11=a11-pi
endif
if(a10>=-top_y_limit.and.a10<=top_y_limit)then
call num(a10,a11,k1,k2,dlanda,dphi,n_epsilon)
he_meridion(i,j)=hen(k1,k2)
vx_meridion(i,j)=-vxn(k1,k2)
vy_meridion(i,j)=-vyn(k1,k2)
vxcosy_meridion(i,j)=-vxn(k1,k2)*dcos(a10)
vycosy_meridion(i,j)=-vyn(k1,k2)*dcos(a10)
else
call yin_to_yang(a10,a11,a12,a13)
call new_interpolation(a12,a13,hee,vxe,vye,nx,my,hem,vxm,vym,dlanda,dphi,n_epsilon,int_method)
call yangvector_to_yinvector(a10,a11,a12,a13,vy,vx,vym,vxm)
he_meridion(i,j)=hem
vx_meridion(i,j)=-vx
vy_meridion(i,j)=-vy
vxcosy_meridion(i,j)=-vx*dcos(a10)
vycosy_meridion(i,j)=-vy*dcos(a10)
endif
enddo
enddo
!do i=1,nx
!vx_meridion(i,1)=(vx_meridion(i,2)+vx_meridion(i,2*(ny-1)))/2.0d0
!vx_meridion(i,ny)=(vx_meridion(i,ny-1)+vx_meridion(i,ny+1))/2.0d0
!vy_meridion(i,1)=(vy_meridion(i,2)+vy_meridion(i,2*(ny-1)))/2.0d0
!vy_meridion(i,ny)=(vy_meridion(i,ny-1)+vy_meridion(i,ny+1))/2.0d0
!vxcosy_meridion(i,1)=0.0d0
!vxcosy_meridion(i,ny)=0.0d0
!vycosy_meridion(i,1)=0.0d0
!vycosy_meridion(i,ny)=0.0d0
!enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  provide_Yang_vector_to_derive(hen,hee,vxn,vyn,vxe,vye,dlanda,dphi,nx,ny,my,n_epsilon&
,he_meridion,vx_meridion,vxcosy_meridion,vy_meridion,vycosy_meridion,int_method)
implicit none
integer::i,j,nx,ny,my,n_epsilon,k1,k2,kney
real*8::hen(nx,my),hee(nx,my),vxn(nx,my),vyn(nx,my),vxe(nx,my),vye(nx,my),vx_meridion(nx,2*(ny-1))&
,vy_meridion(nx,2*(ny-1)),a10,a11,a12,a13,hem,vxm,vym,vx,vy,pi,dlanda,dphi&
,top_x_limit,top_y_limit,vxcosy_meridion(nx,2*(ny-1)),vycosy_meridion(nx,2*(ny-1))&
,he_meridion(nx,2*(ny-1))
character(5)::int_method



kney=int((real(ny-my)+0.01)/2.0)

pi=2.0d0*dasin(1.0d0)
top_x_limit=pi
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)

do i=1,nx
a13=-top_x_limit+dble(i-1)*dlanda
do j=1,ny
a12=-(pi/2.0d0)+dble(j-1)*dphi
if(a12>=-top_y_limit.and.a12<=top_y_limit)then
he_meridion(i,j)=hee(i,j-kney)
vx_meridion(i,j)=vxe(i,j-kney)
vy_meridion(i,j)=vye(i,j-kney)
vxcosy_meridion(i,j)=vxe(i,j-kney)*dcos(a12)
vycosy_meridion(i,j)=vye(i,j-kney)*dcos(a12)
else

call yang_to_yin(a10,a11,a12,a13)
call new_interpolation(a10,a11,hen,vxn,vyn,nx,my,hem,vxm,vym,dlanda,dphi,n_epsilon,int_method)
call yinvector_to_yangvector(a10,a11,a12,a13,vym,vxm,vy,vx)
he_meridion(i,j)=hem
vx_meridion(i,j)=vx
vy_meridion(i,j)=vy
vxcosy_meridion(i,j)=vx*dcos(a12)
vycosy_meridion(i,j)=vy*dcos(a12)
endif
enddo
do j=ny+1,2*(ny-1)
a12=(pi/2.0d0)-dble(j-ny)*dphi
if(a13<0.0d0)then
a13=a13+pi
elseif(a13>=0.0d0)then
a13=a13-pi
endif
if(a12>=-top_y_limit.and.a12<=top_y_limit)then
call num(a12,a13,k1,k2,dlanda,dphi,n_epsilon)
he_meridion(i,j)=hee(k1,k2)
vx_meridion(i,j)=-vxe(k1,k2)
vy_meridion(i,j)=-vye(k1,k2)
vxcosy_meridion(i,j)=-vxe(k1,k2)*dcos(a12)
vycosy_meridion(i,j)=-vye(k1,k2)*dcos(a12)
else
call yang_to_yin(a10,a11,a12,a13)
call new_interpolation(a10,a11,hen,vxn,vyn,nx,my,hem,vxm,vym,dlanda,dphi,n_epsilon,int_method)
call yinvector_to_yangvector(a10,a11,a12,a13,vym,vxm,vy,vx)
he_meridion(i,j)=hem
vx_meridion(i,j)=-vx
vy_meridion(i,j)=-vy
vxcosy_meridion(i,j)=-vx*dcos(a12)
vycosy_meridion(i,j)=-vy*dcos(a12)
endif
enddo
enddo
!do i=1,nx
!vx_meridion(i,1)=(vx_meridion(i,2)+vx_meridion(i,2*(ny-1)))/2.0d0
!vx_meridion(i,ny)=(vx_meridion(i,ny-1)+vx_meridion(i,ny+1))/2.0d0
!vy_meridion(i,1)=(vy_meridion(i,2)+vy_meridion(i,2*(ny-1)))/2.0d0
!vy_meridion(i,ny)=(vy_meridion(i,ny-1)+vy_meridion(i,ny+1))/2.0d0
!vxcosy_meridion(i,1)=0.0d0
!vxcosy_meridion(i,ny)=0.0d0
!vycosy_meridion(i,1)=0.0d0
!vycosy_meridion(i,ny)=0.0d0
!enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine provide_cos_phi(dphi,ny,cosphi)
implicit none
integer::ny,j
real*8::pi,dphi,cosphi(2*(ny-1)),a10,a12
pi=2.d0*dasin(1.0d0)


do j=1,2*(ny-1)
a10=-(pi/2.0d0)+dble(j-1)*dphi
if(a10>pi/2.0d0.and.a10<3.0d0*pi/2.0d0)then
a12=pi-a10
else
a12=a10
endif
a10=a12
cosphi(j)=dcos(a10)
enddo
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine new_interpolation(phi3,landa3,heb,vxb,vyb,nx,ny,hem,vxm,vym,dlanda,dphi,n_epsilon,int_method)
implicit none
integer::nx,ny,n_epsilon,k1,k2
real*8::landa3,phi3,dphi,dlanda,heb(nx,ny),vxb(nx,ny),vyb(nx,ny),phi,landa,dd,ddd&
,top_x_limit,top_y_limit,pi,h1,h2,h3,h4,p1,p2,p3,p4,hem1,hem2,hem3,hem4,pxm1,pxm2,pxm3&
,pxm4,pym1,pym2,pym3,pym4,hem,pxm,pym,vxm,vym
character(5)::int_method
pi=2.d0*dasin(1.0d0)
top_x_limit=pi
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)


k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1

landa=dble(k1-1)*dlanda-top_x_limit
phi=(dble(k2-1))*dphi-top_y_limit


if(int_method=='bilin')then
call distance(phi,phi,landa,landa+dlanda,dd)
call distance(phi,phi,landa3,landa,ddd)
h1=heb(k1,k2)
h2=heb(k1+1,k2)
call linear(h1,h2,dd,ddd,hem1)
p1=heb(k1,k2)*vxb(k1,k2)
p2=heb(k1+1,k2)*vxb(k1+1,k2)
call linear(p1,p2,dd,ddd,pxm1)
p1=heb(k1,k2)*vyb(k1,k2)
p2=heb(k1+1,k2)*vyb(k1+1,k2)
call linear(p1,p2,dd,ddd,pym1)
call distance(phi+dphi,phi+dphi,landa,landa+dlanda,dd)
call distance(phi+dphi,phi+dphi,landa3,landa,ddd)
h1=heb(k1,k2+1)
h2=heb(k1+1,k2+1)
call linear(h1,h2,dd,ddd,hem2)
p1=heb(k1,k2+1)*vxb(k1,k2+1)
p2=heb(k1+1,k2+1)*vxb(k1+1,k2+1)
call linear(p1,p2,dd,ddd,pxm2)
p1=heb(k1,k2+1)*vyb(k1,k2+1)
p2=heb(k1+1,k2+1)*vyb(k1+1,k2+1)
call linear(p1,p2,dd,ddd,pym2)
call distance(phi,phi+dphi,landa3,landa3,dd)
call distance(phi,phi3,landa3,landa3,ddd)
call linear(hem1,hem2,dd,ddd,hem)
call linear(pxm1,pxm2,dd,ddd,pxm)
vxm=pxm/hem
call linear(pym1,pym2,dd,ddd,pym)
vym=pym/hem

elseif(int_method=='bicub')then
call distance(phi-dphi,phi-dphi,landa,landa+dlanda,dd)
call distance(phi-dphi,phi-dphi,landa3,landa,ddd)
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem1)
p1=heb(k1-1,k2-1)*vxb(k1-1,k2-1)
p2=heb(k1,k2-1)*vxb(k1,k2-1)
p3=heb(k1+1,k2-1)*vxb(k1+1,k2-1)
p4=heb(k1+2,k2-1)*vxb(k1+2,k2-1)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pxm1)
p1=heb(k1-1,k2-1)*vyb(k1-1,k2-1)
p2=heb(k1,k2-1)*vyb(k1,k2-1)
p3=heb(k1+1,k2-1)*vyb(k1+1,k2-1)
p4=heb(k1+2,k2-1)*vyb(k1+2,k2-1)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pym1)


call distance(phi,phi,landa,landa+dlanda,dd)
call distance(phi,phi,landa3,landa,ddd)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem2)
p1=heb(k1-1,k2)*vxb(k1-1,k2)
p2=heb(k1,k2)*vxb(k1,k2)
p3=heb(k1+1,k2)*vxb(k1+1,k2)
p4=heb(k1+2,k2)*vxb(k1+2,k2)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pxm2)
p1=heb(k1-1,k2)*vyb(k1-1,k2)
p2=heb(k1,k2)*vyb(k1,k2)
p3=heb(k1+1,k2)*vyb(k1+1,k2)
p4=heb(k1+2,k2)*vyb(k1+2,k2)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pym2)

call distance(phi+dphi,phi+dphi,landa,landa+dlanda,dd)
call distance(phi+dphi,phi+dphi,landa3,landa,ddd)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem3)
p1=heb(k1-1,k2+1)*vxb(k1-1,k2+1)
p2=heb(k1,k2+1)*vxb(k1,k2+1)
p3=heb(k1+1,k2+1)*vxb(k1+1,k2+1)
p4=heb(k1+2,k2+1)*vxb(k1+2,k2+1)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pxm3)
p1=heb(k1-1,k2+1)*vyb(k1-1,k2+1)
p2=heb(k1,k2+1)*vyb(k1,k2+1)
p3=heb(k1+1,k2+1)*vyb(k1+1,k2+1)
p4=heb(k1+2,k2+1)*vyb(k1+2,k2+1)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pym3)


call distance(phi+(2.0d0*dphi),phi+(2.0d0*dphi),landa,landa+dlanda,dd)
call distance(phi+(2.0d0*dphi),phi+(2.0d0*dphi),landa3,landa,ddd)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem4)
p1=heb(k1-1,k2+2)*vxb(k1-1,k2+2)
p2=heb(k1,k2+2)*vxb(k1,k2+2)
p3=heb(k1+1,k2+2)*vxb(k1+1,k2+2)
p4=heb(k1+2,k2+2)*vxb(k1+2,k2+2)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pxm4)
p1=heb(k1-1,k2+2)*vyb(k1-1,k2+2)
p2=heb(k1,k2+2)*vyb(k1,k2+2)
p3=heb(k1+1,k2+2)*vyb(k1+1,k2+2)
p4=heb(k1+2,k2+2)*vyb(k1+2,k2+2)
call cubic(p1,p2,p3,p4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pym4)

call distance(phi,phi+dphi,landa3,landa3,dd)
call distance(phi,phi3,landa3,landa3,ddd)
call cubic(hem1,hem2,hem3,hem4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem)
call cubic(pxm1,pxm2,pxm3,pxm4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pxm)
vxm=pxm/hem
call cubic(pym1,pym2,pym3,pym4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,pym)
vym=pym/hem
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine cubic(h1,h2,h3,h4,x1,x2,x3,x4,xm,hm)
implicit none
real*8::h1,h2,h3,h4,x1,x2,x3,x4,xm,hm
real*8::a,b,c,d,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11
c1=(h4-h1)/(x4-x1)
c2=(h2-h1)/(x2-x1)
c3=((h3-h1)*(x4-x2))/((x3-x1)*(x3-x2))
c4=((h2-h1)*(x4-x2))/((x2-x1)*(x3-x2))
c5=(x4*x4)-(x2*x4)-(x3*x4)+(x2*x3)
a=(c1/c5)-(c2/c5)-(c3/c5)+(c4/c5)
c6=(h3-h1)/((x3-x1)*(x3-x2))
c7=(h2-h1)/((x2-x1)*(x3-x2))
c8=x1+x2+x3
b=c6-c7-(c8*a)
c9=(h2-h1)/(x2-x1)
c10=(x2*x2)+(x1*x1)+(x2*x1)
c11=x1+x2
c=c9-(c10*a)-(b*c11)
d=h1-(a*x1*x1*x1)-(b*x1*x1)-(c*x1)
hm=(a*xm*xm*xm)+(b*xm*xm)+(c*xm)+d
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine interpolation_scalar(phi3,landa3,heb,nx,ny,hem,dlanda,dphi,n_epsilon,int_method)
implicit none
integer::nx,ny,n_epsilon,k1,k2
real*8::landa3,phi3,dphi,dlanda,heb(nx,ny),phi,landa,dd,ddd&
,top_x_limit,top_y_limit,pi,h1,h2,h3,h4,hem1,hem2,hem3,hem4,hem
character(5)::int_method
pi=2.d0*dasin(1.0d0)
top_x_limit=pi
top_y_limit=(pi/4.0d0)+(dble(n_epsilon)*dphi)


k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1

landa=dble(k1-1)*dlanda-top_x_limit
phi=(dble(k2-1))*dphi-top_y_limit


if(int_method=='bilin')then
call distance(phi,phi,landa,landa+dlanda,dd)
call distance(phi,phi,landa3,landa,ddd)
h1=heb(k1,k2)
h2=heb(k1+1,k2)
call linear(h1,h2,dd,ddd,hem1)
call distance(phi+dphi,phi+dphi,landa,landa+dlanda,dd)
call distance(phi+dphi,phi+dphi,landa3,landa,ddd)
h1=heb(k1,k2+1)
h2=heb(k1+1,k2+1)


call distance(phi,phi+dphi,landa3,landa3,dd)
call distance(phi,phi3,landa3,landa3,ddd)
call linear(hem1,hem2,dd,ddd,hem)


elseif(int_method=='bicub')then
call distance(phi-dphi,phi-dphi,landa,landa+dlanda,dd)
call distance(phi-dphi,phi-dphi,landa3,landa,ddd)
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem1)



call distance(phi,phi,landa,landa+dlanda,dd)
call distance(phi,phi,landa3,landa,ddd)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem2)


call distance(phi+dphi,phi+dphi,landa,landa+dlanda,dd)
call distance(phi+dphi,phi+dphi,landa3,landa,ddd)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem3)


call distance(phi+(2.0d0*dphi),phi+(2.0d0*dphi),landa,landa+dlanda,dd)
call distance(phi+(2.0d0*dphi),phi+(2.0d0*dphi),landa3,landa,ddd)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem4)


call distance(phi,phi+dphi,landa3,landa3,dd)
call distance(phi,phi3,landa3,landa3,ddd)
call cubic(hem1,hem2,hem3,hem4,0.0d0,dd,2.0d0*dd,3.0d0*dd,dd+ddd,hem)

endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bicubic(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y)
implicit none
integer::nx,ny,n_epsilon_x,n_epsilon_y
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hm,hm1,hm2,phi,landa
real*8::heb(nx,ny),hm3,hm4,h1,h2,h3,h4

pi=2.d0*dasin(1.0d0)

k1=int(((landa3+(3.0d0*pi/4.0d0)+1.0d-10)/dlanda))+1+n_epsilon_x
k2=int(((phi3+(pi/4.0d0)+1.0d-10)/dphi))+1+n_epsilon_y
phi=(dble(k2-1-n_epsilon_y))*dphi-(pi/4.0d0)
landa=dble(k1-1-n_epsilon_x)*dlanda-3.0d0*pi/4.0d0
if(k1==1)then
h1=heb(nx,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(nx,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(nx,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(nx,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
elseif(k1==nx-1)then
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(1,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(1,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(1,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(1,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
elseif(k1==nx)then
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(1,k2-1)
h4=heb(2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(1,k2)
h4=heb(2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(1,k2+1)
h4=heb(2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(1,k2+2)
h4=heb(2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
else
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
endif
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Trans_Deriv_Yin_to_Yang(phin,landan,phie,landae,dphin,dlandan,dphie,dlandae)
implicit none
real*8::phin,landan,phie,landae,dlandan,dphin,dlandae,dphie
if(dabs(dcos(phin))<=1.0d-8)then
dlandae=dlandan
dphie=dphin
else
dlandae=-((dsin(phin)*dsin(phie)/(dcos(phin)*dcos(phin)))*dlandan)+(dcos(landan)*dphin)
dphie=((dcos(landae)/(dcos(phin)*dcos(phin)))*dlandan)+((dsin(landae)*dsin(phie)/dcos(phin))*dphin)
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Trans_Deriv_Yang_to_Yin(phin,landan,phie,landae,dphin,dlandan,dphie,dlandae)
implicit none
real*8::phin,landan,phie,landae,dlandan,dphin,dlandae,dphie
if(dabs(dcos(phie))<=1.0d-8)then
dlandan=dlandae
dphin=dphie
else
dlandan=-((dsin(phie)*dsin(phin)/(dcos(phie)*dcos(phie)))*dlandae)+(dcos(landae)*dphie)
dphin=((dcos(landan)/(dcos(phie)*dcos(phie)))*dlandae)+((dsin(landan)*dsin(phin)/dcos(phie))*dphie)
endif
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Trans_vector_Deriv_Yin_to_Yang(phin,landan,phie,landae,un,vn,dunx,dvnx,duny,dvny,duex,dvex,duey,dvey)
implicit none
real*8::phin,landan,phie,landae,dunx,dvnx,duny,dvny,duex,dvex,duey,dvey,cossaay,sinsaay,&
dcosx,dcosy,dsinx,dsiny,dxm,dym,un,vn

!call yangvector_to_yinvector(phin,landan,phie,landae,vn,un,ve,ue)
call Trans_Deriv_Yin_to_Yang(phin,landan,phie,landae,duny,dunx,dym,dxm)
dunx=dxm
duny=dym
call Trans_Deriv_Yin_to_Yang(phin,landan,phie,landae,dvny,dvnx,dym,dxm)
dvnx=dxm
dvny=dym
cossaay=-dsin(landae)*dsin(landan)
if(dabs(dcos(phie))>=dabs(dcos(landan)))then
sinsaay=dcos(landan)/dcos(phie)
else
sinsaay=-dcos(landae)/dcos(phin)
endif
if(dabs(dcos(phin))<1.0d-8)then
dcosx=(-dcos(landae)*dsin(landan))
dcosy=0.0d0
dsinx=0.0d0
dsiny=0.0d0
else
dcosx=-dcos(landae)*dsin(landan)+((dsin(landae)*dcos(landan)*dsin(phin)*dsin(phie))/(dcos(phin)*dcos(phin)))
dcosy=-((dsin(landae)*dcos(landan)*dcos(landae))/(dcos(phin)*dcos(phin)))
dsinx=(dsin(landan)*dsin(landan)*dsin(landae))/(dcos(phin))
dsiny=(dsin(landae)*dsin(landae)*dsin(landan)*dcos(landan))/(dcos(phin))
endif



duex=un*dcosx+cossaay*dunx+vn*dsinx+sinsaay*dvnx
dvex=-un*dsinx-sinsaay*dunx+vn*dcosx+cossaay*dvnx
duey=un*dcosy+cossaay*duny+vn*dsiny+sinsaay*dvny
dvey=-un*dsiny-sinsaay*duny+vn*dcosy+cossaay*dvny
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  Trans_vector_Deriv_Yang_to_Yin(phin,landan,phie,landae,ue,ve,dunx,dvnx,duny,dvny,duex,dvex,duey,dvey)
implicit none
real*8::phin,landan,phie,landae,dunx,dvnx,duny,dvny,duex,dvex,duey,dvey,cossaay,sinsaay,&
dcosx,dcosy,dsinx,dsiny,dxm,dym,ue,ve
!call yinvector_to_yangvector(phin,landan,phie,landae,vn,un,ve,ue)
call Trans_Deriv_Yang_to_Yin(phin,landan,phie,landae,dym,dxm,duey,duex)
duex=dxm
duey=dym
call Trans_Deriv_Yang_to_Yin(phin,landan,phie,landae,dym,dxm,dvey,dvex)
dvex=dxm
dvey=dym
cossaay=-dsin(landae)*dsin(landan)
if(dabs(dcos(phin))>=dabs(dcos(landae)))then
sinsaay=-dcos(landae)/dcos(phin)
else
sinsaay=dcos(landan)/dcos(phie)
endif
if(dabs(dcos(phin))<1.0d-8)then
dcosx=(-dcos(landae)*dsin(landan))
dcosy=0.0d0
dsinx=0.0d0
dsiny=0.0d0
else
dcosx=-dcos(landan)*dsin(landae)+((dsin(landan)*dcos(landae)*dsin(phin)*dsin(phie))/(dcos(phie)*dcos(phie)))
dcosy=-((dsin(landan)*dcos(landan)*dcos(landae))/(dcos(phie)*dcos(phie)))
dsinx=-(dsin(landae)*dsin(landae)*dsin(landan))/(dcos(phie))
dsiny=-(dsin(landan)*dsin(landan)*dsin(landae)*dcos(landae))/(dcos(phie))
endif
dunx=ue*dcosx+cossaay*duex-ve*dsinx-sinsaay*dvex
dvnx=ue*dsinx+sinsaay*duex+ve*dcosx+cossaay*dvex
duny=ue*dcosy+cossaay*duey-ve*dsiny-sinsaay*dvey
dvny=ue*dsiny+sinsaay*duey+ve*dcosy+cossaay*dvey
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bilinear(phi3,landa3,heb,nx,ny,hem,dlanda,dphi,n_epsilon_x,n_epsilon_y)
implicit none
integer::nx,ny,n_epsilon_x,n_epsilon_y
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hem,dl,dp,phi,landa,hem1,hem2
real*8::heb(nx,ny),top_x_limit,top_y_limit


pi=2.d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)


k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1

landa=dble(k1-1)*dlanda-top_x_limit
phi=(dble(k2-1))*dphi-top_y_limit
dl=landa3-landa
dp=phi3-phi
if(k1==nx)then
call linear(heb(k1,k2),heb(1,k2),dlanda,dl,hem1)
call linear(heb(k1,k2+1),heb(1,k2+1),dlanda,dl,hem2)
call linear(hem1,hem2,dphi,dp,hem)
else
call linear(heb(k1,k2),heb(k1+1,k2),dlanda,dl,hem1)
call linear(heb(k1,k2+1),heb(k1+1,k2+1),dlanda,dl,hem2)
call linear(hem1,hem2,dphi,dp,hem)
endif

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine linear(h1,h2,d,dd,hem)
implicit none
real*8::hem,h1,h2,d,dd
hem=((h1*(d-dd))+(h2*dd))/d
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bicubic_glob(phi3,landa3,heb,nx,ny,hm,dlanda,dphi)
implicit none
integer::nx,ny
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hm,hm1,hm2,phi,landa
real*8::heb(nx,ny),hm3,hm4,h1,h2,h3,h4,dl,dp

pi=2.d0*dasin(1.0d0)

k1=int(((landa3+pi+1.0d-10)/dlanda))+1
k2=int(((phi3+(pi/2.0d0)+1.0d-10)/dphi))+1
phi=(dble(k2-1))*dphi-(pi/2.0d0)
landa=dble(k1-1)*dlanda-pi
if(k2==1.and.k1==nx)then
dl=landa3-landa
dp=phi3-phi
call linear(heb(k1,k2),heb(1,k2),dlanda,dl,hm1)
call linear(heb(k1,k2+1),heb(1,k2+1),dlanda,dl,hm2)
call linear(hm1,hm2,dphi,dp,hm)
elseif(k2==1)then
dl=landa3-landa
dp=phi3-phi
call linear(heb(k1,k2),heb(k1+1,k2),dlanda,dl,hm1)
call linear(heb(k1,k2+1),heb(k1+1,k2+1),dlanda,dl,hm2)
call linear(hm1,hm2,dphi,dp,hm)
elseif(k2==ny.and.k1==1)then
h1=heb(nx,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm)
elseif(k2==ny.and.k1==nx-1)then
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(1,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm)
elseif(k2==ny.and.k1==nx)then
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(1,k2)
h4=heb(2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm)
elseif(k2==ny)then
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm)
elseif(k2==ny-1.and.k1==nx)then
dl=landa3-landa
dp=phi3-phi
call linear(heb(k1,k2),heb(1,k2),dlanda,dl,hm1)
call linear(heb(k1,k2+1),heb(1,k2+1),dlanda,dl,hm2)
call linear(hm1,hm2,dphi,dp,hm)
elseif(k2==ny-1)then
dl=landa3-landa
dp=phi3-phi
call linear(heb(k1,k2),heb(k1+1,k2),dlanda,dl,hm1)
call linear(heb(k1,k2+1),heb(k1+1,k2+1),dlanda,dl,hm2)
call linear(hm1,hm2,dphi,dp,hm)
elseif(k1==1)then
h1=heb(nx,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(nx,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(nx,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(nx,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
elseif(k1==nx)then
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(1,k2-1)
h4=heb(2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(1,k2)
h4=heb(2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(1,k2+1)
h4=heb(2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(1,k2+2)
h4=heb(2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
elseif(k1==nx-1)then
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(1,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(1,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(1,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(1,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
else
h1=heb(k1-1,k2-1)
h2=heb(k1,k2-1)
h3=heb(k1+1,k2-1)
h4=heb(k1+2,k2-1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm1)
h1=heb(k1-1,k2)
h2=heb(k1,k2)
h3=heb(k1+1,k2)
h4=heb(k1+2,k2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm2)
h1=heb(k1-1,k2+1)
h2=heb(k1,k2+1)
h3=heb(k1+1,k2+1)
h4=heb(k1+2,k2+1)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm3)
h1=heb(k1-1,k2+2)
h2=heb(k1,k2+2)
h3=heb(k1+1,k2+2)
h4=heb(k1+2,k2+2)
call cubic(h1,h2,h3,h4,landa-dlanda,landa,landa+dlanda,landa+2.0d0*dlanda,landa3,hm4)
call cubic(hm1,hm2,hm3,hm4,phi-dphi,phi,phi+dphi,phi+2.0d0*dphi,phi3,hm)
endif
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  tu(u,he,nx,ny,kkx,kky)
implicit none
integer::i,j,nx,ny,kkx,kky
real*8::u(kkx,kky),he(nx,ny)
do j=1,ny
do i=1,kkx

u(i,j)=he(i,j)
u(i,kky-j+1)=he(i+kkx,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  td(d2,d1,nx,ny,kkx,kky)
implicit none
integer::i,j,nx,ny,kkx,kky
real*8::d1(kkx,kky),d2(nx,ny)
do j=1,ny
do i=1,kkx

d2(i,j)=d1(i,j)
d2(i+kkx,j)=-d1(i,kky-j+1)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  tuv(u,v,nx,ny,kkx,kky,phi)
implicit none
integer::i,j,nx,ny,kkx,kky
real*8::u(kkx,kky),v(nx,ny),phi(ny)
do i=1,kkx
do j=1,ny
u(i,j)=v(i,j)*dcos(phi(j))
u(i,kky-j+1)=-v(i+kkx,j)*dcos(phi(j))
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  tdv(d2,d1,nx,ny,kkx,kky)
implicit none
integer::i,j,nx,ny,kkx,kky
real*8::d1(kkx,kky),d2(nx,ny)
do j=1,ny
do i=1,kkx

d2(i,j)=d1(i,j)
d2(i+kkx,j)=d1(i,kky-j+1)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  tuu(u,v,nx,ny,kkx,kky)
implicit none
integer::i,j,nx,ny,kkx,kky
real*8::u(kkx,kky),v(nx,ny)
do j=1,ny
do i=1,kkx

u(i,j)=v(i,j)
u(i,kky-j+1)=-v(i+kkx,j)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  tdu(d2,d1,nx,ny,kkx,kky)
implicit none
integer::i,j,nx,ny,kkx,kky
real*8::d1(kkx,kky),d2(nx,ny)
do i=1,kkx
do j=1,ny
d2(i,j)=d1(i,j)
d2(i+kkx,j)=d1(i,kky-j+1)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  derivative1(hene,u,v,hee,ue,ve,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)
implicit none
integer::i,j,nx,ny,mx,my,kkx,kky,knex,kney,n_epsilon_x,n_epsilon_y,l_pole_variate
real*8::hene(nx,ny),u(nx,ny),v(nx,ny),dux(nx,ny),dvx(nx,ny),duy(nx,ny),dvy(nx,ny),dhx(nx,ny),dhy(nx,ny)&
,dunx(mx,my),dvnx(mx,my),duny(mx,my),dvny(mx,my),dhnx(mx,my),dhny(mx,my)&
,duex(mx,my),dvex(mx,my),duey(mx,my),dvey(mx,my),dhex(mx,my),dhey(mx,my)&
,a10,a11,a12,a13,hmnx,hmny,hmex,hmey,umnx,umny,vmnx,vmny,pi,dlanda,dphi,dccy(kkx,kky)&
,top_x_limit,top_y_limit,u1(kkx,kky),umex,umey,vmex,vmey,ue(mx,my),ve(mx,my),vn,un,hee(mx,my),a40
!,dux_b(my),dux_f(my),dvx_b(my),dvx_f(my),dhx_b(my),dhx_f(my),duy_b(mx),duy_f(mx),dvy_b(mx),dvy_f(mx)&
!,dhy_b(mx),dhy_f(mx)
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)


call tuu(u1,u,nx,ny,kkx,kky)
call dby42p(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call tdu(duy,dccy,nx,ny,kkx,kky)
call dfx42p(u,dux,nx,ny,dlanda)

call tuu(u1,v,nx,ny,kkx,kky)
call dby42p(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call tdu(dvy,dccy,nx,ny,kkx,kky)
call dfx42p(v,dvx,nx,ny,dlanda)

call tu(u1,hene,nx,ny,kkx,kky)
call dby42p(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call td(dhy,dccy,nx,ny,kkx,kky)
call dfx42p(hene,dhx,nx,ny,dlanda)
do i=1,nx
dux(i,1)=0.0d0
dvx(i,1)=0.0d0
dhx(i,1)=0.0d0
dux(i,ny)=0.0d0
dvx(i,ny)=0.0d0
dhx(i,ny)=0.0d0
enddo
do j=1,my
a12=-top_y_limit+dble(j-1)*dphi
do i=1,mx
dunx(i,j)=dux(i+knex,j+kney)
dvnx(i,j)=dvx(i+knex,j+kney)
duny(i,j)=duy(i+knex,j+kney)
dvny(i,j)=dvy(i+knex,j+kney)
dhnx(i,j)=dhx(i+knex,j+kney)
dhny(i,j)=dhy(i+knex,j+kney)

a13=-top_x_limit+dble(i-1)*dlanda
call yang_to_yin(a10,a11,a12,a13)
if(dabs(dcos(a10))<=1.0d-8)then
duex(i,j)=0.0d0
dvex(i,j)=0.0d0
duey(i,j)=0.0d0
dvey(i,j)=0.0d0
dhex(i,j)=0.0d0
dhey(i,j)=0.0d0
else
call select_int_method_glob(a10,a11,dhx,nx,ny,hmnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dhy,nx,ny,hmny,dlanda,dphi,int_method)
call Trans_Deriv_Yin_to_Yang(a10,a11,a12,a13,hmny,hmnx,hmey,hmex)
dhex(i,j)=hmex
dhey(i,j)=hmey
call select_int_method_glob(a10,a11,dux,nx,ny,umnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,duy,nx,ny,umny,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dvx,nx,ny,vmnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dvy,nx,ny,vmny,dlanda,dphi,int_method)
call yangvector_to_yinvector(a10,a11,a12,a13,vn,un,ve(i,j),ue(i,j))
call Trans_vector_Deriv_Yin_to_Yang(a10,a11,a12,a13,un,vn,umnx,vmnx,umny,vmny,umex,vmex,umey,vmey)
duex(i,j)=umex
duey(i,j)=umey
dvex(i,j)=vmex
dvey(i,j)=vmey
endif
enddo
enddo
!***************************************************************************
!pole problem solution 1
!***************************************************************************
a40=0.5d0-0.5d0/dsqrt(3.0d0)
if(nx==256)then
l_pole_variate=3
elseif(nx==512)then
l_pole_variate=4
elseif(nx==128)then
l_pole_variate=3
elseif(nx==64)then
l_pole_variate=1
elseif(nx==1024)then
l_pole_variate=5
elseif(nx==32)then
l_pole_variate=1
else
l_pole_variate=3
endif
do j=(my+1)/2-l_pole_variate,(my+1)/2+l_pole_variate
do i=(nx/4)+1-knex+l_pole_variate,(nx/4)+1-knex-l_pole_variate,-1
duex(i,j)=(((ue(i+1,j)-ue(i,j))/dlanda)-(a40*duex(i+1,j)))/(1.0d0-a40)
dvex(i,j)=(((ve(i+1,j)-ve(i,j))/dlanda)-(a40*dvex(i+1,j)))/(1.0d0-a40)
dhex(i,j)=(((hee(i+1,j)-hee(i,j))/dlanda)-(a40*dhex(i+1,j)))/(1.0d0-a40)
enddo
enddo
do i=(nx/4)+1-knex-l_pole_variate,(nx/4)+1-knex+l_pole_variate
do j=(my+1)/2-l_pole_variate,(my+1)/2+l_pole_variate
duey(i,j)=(((ue(i,j)-ue(i,j-1))/dphi)-(a40*duey(i,j-1)))/(1.0d0-a40)
dvey(i,j)=(((ve(i,j)-ve(i,j-1))/dphi)-(a40*dvey(i,j-1)))/(1.0d0-a40)
dhey(i,j)=(((hee(i,j)-hee(i,j-1))/dphi)-(a40*dhey(i,j-1)))/(1.0d0-a40)
enddo
enddo
do j=(my+1)/2-l_pole_variate,(my+1)/2+l_pole_variate
do i=(3*nx/4)+1-knex+l_pole_variate,(3*nx/4)+1-knex-l_pole_variate,-1
duex(i,j)=(((ue(i+1,j)-ue(i,j))/dlanda)-(a40*duex(i+1,j)))/(1.0d0-a40)
dvex(i,j)=(((ve(i+1,j)-ve(i,j))/dlanda)-(a40*dvex(i+1,j)))/(1.0d0-a40)
dhex(i,j)=(((hee(i+1,j)-hee(i,j))/dlanda)-(a40*dhex(i+1,j)))/(1.0d0-a40)
enddo
enddo
do i=(3*nx/4)+1-knex-l_pole_variate,(3*nx/4)+1-knex+l_pole_variate
do j=(my+1)/2-l_pole_variate,(my+1)/2+l_pole_variate
duey(i,j)=(((ue(i,j)-ue(i,j-1))/dphi)-(a40*duey(i,j-1)))/(1.0d0-a40)
dvey(i,j)=(((ve(i,j)-ve(i,j-1))/dphi)-(a40*dvey(i,j-1)))/(1.0d0-a40)
dhey(i,j)=(((hee(i,j)-hee(i,j-1))/dphi)-(a40*dhey(i,j-1)))/(1.0d0-a40)
enddo
enddo
!***************************************************************************
!pole problem solution 2
!***************************************************************************
!do j=1,my
!dux_b(j)=duex(1,j)
!dux_f(j)=duex(mx,j)
!dvx_b(j)=dvex(1,j)
!dvx_f(j)=dvex(mx,j)
!dhx_b(j)=dhex(1,j)
!dhx_f(j)=dhex(mx,j)
!enddo
!do i=1,mx

!duy_b(i)=duey(i,1)
!duy_f(i)=duey(i,my)
!dvy_b(i)=dvey(i,1)
!dvy_f(i)=dvey(i,my)
!dhy_b(i)=dhey(i,1)
!dhy_f(i)=dhey(i,my)

!enddo

!call dfx42(ue,duex,mx,my,dlanda,dux_f)
!call dfx42(ve,dvex,mx,my,dlanda,dvx_f)
!call dfx42(hee,dhex,mx,my,dlanda,dhx_f)

!call dby42(ue,duey,mx,my,dphi,duy_b)
!call dby42(ve,dvey,mx,my,dphi,dvy_b)

!call dby42(hee,dhey,mx,my,dphi,dhy_b)
!***************************************************************************
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  derivative2(hene,u,v,hee,ue,ve,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dhx,dhy,dunx,dvnx,duny&
,dvny,dhnx,dhny,duex,dvex,duey,dvey,dhex,dhey,int_method,n_epsilon_x,n_epsilon_y)
implicit none
integer::i,j,nx,ny,mx,my,kkx,kky,knex,kney,n_epsilon_x,n_epsilon_y,l_pole_variate
real*8::hene(nx,ny),u(nx,ny),v(nx,ny),dux(nx,ny),dvx(nx,ny),duy(nx,ny),dvy(nx,ny),dhx(nx,ny),dhy(nx,ny)&
,dunx(mx,my),dvnx(mx,my),duny(mx,my),dvny(mx,my),dhnx(mx,my),dhny(mx,my)&
,duex(mx,my),dvex(mx,my),duey(mx,my),dvey(mx,my),dhex(mx,my),dhey(mx,my)&
,a10,a11,a12,a13,hmnx,hmny,hmex,hmey,umnx,umny,vmnx,vmny,pi,dlanda,dphi,dccy(kkx,kky)&
,top_x_limit,top_y_limit,u1(kkx,kky),umex,umey,vmex,vmey,ue(mx,my),ve(mx,my),un,vn,hee(mx,my),a40
!,dux_b(my),dux_f(my),dvx_b(my),dvx_f(my),dhx_b(my),dhx_f(my),duy_b(mx),duy_f(mx),dvy_b(mx),dvy_f(mx)&
!,dhy_b(mx),dhy_f(mx)
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)


call tuu(u1,u,nx,ny,kkx,kky)
call dfy42p(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call tdu(duy,dccy,nx,ny,kkx,kky)
call dbx42p(u,dux,nx,ny,dlanda)

call tuu(u1,v,nx,ny,kkx,kky)
call dfy42p(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call tdu(dvy,dccy,nx,ny,kkx,kky)
call dbx42p(v,dvx,nx,ny,dlanda)

call tu(u1,hene,nx,ny,kkx,kky)
call dfy42p(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call td(dhy,dccy,nx,ny,kkx,kky)
call dbx42p(hene,dhx,nx,ny,dlanda)
do i=1,nx
dux(i,1)=0.0d0
dvx(i,1)=0.0d0
dhx(i,1)=0.0d0
dux(i,ny)=0.0d0
dvx(i,ny)=0.0d0
dhx(i,ny)=0.0d0
enddo
do j=1,my
a12=-top_y_limit+dble(j-1)*dphi
do i=1,mx
dunx(i,j)=dux(i+knex,j+kney)
dvnx(i,j)=dvx(i+knex,j+kney)
duny(i,j)=duy(i+knex,j+kney)
dvny(i,j)=dvy(i+knex,j+kney)
dhnx(i,j)=dhx(i+knex,j+kney)
dhny(i,j)=dhy(i+knex,j+kney)
a13=-top_x_limit+dble(i-1)*dlanda
call yang_to_yin(a10,a11,a12,a13)
if(dabs(dcos(a10))<=1.0d-8)then
duex(i,j)=0.0d0
dvex(i,j)=0.0d0
duey(i,j)=0.0d0
dvey(i,j)=0.0d0
dhex(i,j)=0.0d0
dhey(i,j)=0.0d0
else
call select_int_method_glob(a10,a11,dhx,nx,ny,hmnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dhy,nx,ny,hmny,dlanda,dphi,int_method)
call Trans_Deriv_Yin_to_Yang(a10,a11,a12,a13,hmny,hmnx,hmey,hmex)
dhex(i,j)=hmex
dhey(i,j)=hmey
call select_int_method_glob(a10,a11,dux,nx,ny,umnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,duy,nx,ny,umny,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dvx,nx,ny,vmnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dvy,nx,ny,vmny,dlanda,dphi,int_method)
call yangvector_to_yinvector(a10,a11,a12,a13,vn,un,ve(i,j),ue(i,j))
call Trans_vector_Deriv_Yin_to_Yang(a10,a11,a12,a13,un,vn,umnx,vmnx,umny,vmny,umex,vmex,umey,vmey)
duex(i,j)=umex
duey(i,j)=umey
dvex(i,j)=vmex
dvey(i,j)=vmey
endif
enddo
enddo


!***************************************************************************
!pole problem solution 1
!***************************************************************************
a40=0.5d0-0.5d0/dsqrt(3.0d0)
if(nx==256)then
l_pole_variate=3
elseif(nx==512)then
l_pole_variate=4
elseif(nx==128)then
l_pole_variate=3
elseif(nx==64)then
l_pole_variate=1
elseif(nx==1024)then
l_pole_variate=5
elseif(nx==32)then
l_pole_variate=1
else
l_pole_variate=3
endif
do j=(my+1)/2-l_pole_variate,(my+1)/2+l_pole_variate
do i=(nx/4)+1-knex-l_pole_variate,(nx/4)+1-knex+l_pole_variate
duex(i,j)=(((ue(i,j)-ue(i-1,j))/dlanda)-(a40*duex(i-1,j)))/(1.0d0-a40)
dvex(i,j)=(((ve(i,j)-ve(i-1,j))/dlanda)-(a40*dvex(i-1,j)))/(1.0d0-a40)
dhex(i,j)=(((hee(i,j)-hee(i-1,j))/dlanda)-(a40*dhex(i-1,j)))/(1.0d0-a40)
enddo
enddo
do i=(nx/4)+1-knex-l_pole_variate,(nx/4)+1-knex+l_pole_variate
do j=(my+1)/2+l_pole_variate,(my+1)/2-l_pole_variate,-1
duey(i,j)=(((ue(i,j+1)-ue(i,j))/dphi)-(a40*duey(i,j+1)))/(1.0d0-a40)
dvey(i,j)=(((ve(i,j+1)-ve(i,j))/dphi)-(a40*dvey(i,j+1)))/(1.0d0-a40)
dhey(i,j)=(((hee(i,j+1)-hee(i,j))/dphi)-(a40*dhey(i,j+1)))/(1.0d0-a40)
enddo
enddo
do j=(my+1)/2-l_pole_variate,(my+1)/2+l_pole_variate
do i=(3*nx/4)+1-knex-l_pole_variate,(3*nx/4)+1-knex+l_pole_variate
duex(i,j)=(((ue(i,j)-ue(i-1,j))/dlanda)-(a40*duex(i-1,j)))/(1.0d0-a40)
dvex(i,j)=(((ve(i,j)-ve(i-1,j))/dlanda)-(a40*dvex(i-1,j)))/(1.0d0-a40)
dhex(i,j)=(((hee(i,j)-hee(i-1,j))/dlanda)-(a40*dhex(i-1,j)))/(1.0d0-a40)
enddo
enddo
do i=(3*nx/4)+1-knex-l_pole_variate,(3*nx/4)+1-knex+l_pole_variate
do j=(my+1)/2+l_pole_variate,(my+1)/2-l_pole_variate,-1
duey(i,j)=(((ue(i,j+1)-ue(i,j))/dphi)-(a40*duey(i,j+1)))/(1.0d0-a40)
dvey(i,j)=(((ve(i,j+1)-ve(i,j))/dphi)-(a40*dvey(i,j+1)))/(1.0d0-a40)
dhey(i,j)=(((hee(i,j+1)-hee(i,j))/dphi)-(a40*dhey(i,j+1)))/(1.0d0-a40)
enddo
enddo
!***************************************************************************
!pole problem solution 2
!***************************************************************************
!do j=1,my
!dux_b(j)=duex(1,j)
!dux_f(j)=duex(mx,j)
!dvx_b(j)=dvex(1,j)
!dvx_f(j)=dvex(mx,j)
!dhx_b(j)=dhex(1,j)
!dhx_f(j)=dhex(mx,j)
!enddo
!do i=1,mx
!duy_b(i)=duey(i,1)
!duy_f(i)=duey(i,my)
!dvy_b(i)=dvey(i,1)
!dvy_f(i)=dvey(i,my)
!dhy_b(i)=dhey(i,1)
!dhy_f(i)=dhey(i,my)
!enddo


!call dbx42(ue,duex,mx,my,dlanda,dux_b)
!call dbx42(ve,dvex,mx,my,dlanda,dvx_b)
!call dbx42(hee,dhex,mx,my,dlanda,dhx_b)

!call dfy42(ue,duey,mx,my,dphi,duy_f)
!call dfy42(ve,dvey,mx,my,dphi,dvy_f)
!call dfy42(hee,dhey,mx,my,dphi,dhy_f)
!***************************************************************************
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  derivative3(u,v,ue,ve,dlanda,dphi,nx,ny,mx,my,kkx,kky,dux,dvx,duy,dvy,dunx,dvnx,duny&
,dvny,duex,dvex,duey,dvey,int_method,n_epsilon_x,n_epsilon_y)
implicit none
integer::i,j,nx,ny,mx,my,kkx,kky,knex,kney,n_epsilon_x,n_epsilon_y
real*8::u(nx,ny),v(nx,ny),dux(nx,ny),dvx(nx,ny),duy(nx,ny),dvy(nx,ny)&
,dunx(mx,my),dvnx(mx,my),duny(mx,my),dvny(mx,my)&
,duex(mx,my),dvex(mx,my),duey(mx,my),dvey(mx,my)&
,a10,a11,a12,a13,umnx,umny,vmnx,vmny,pi,dlanda,dphi,dccy(kkx,kky)&
,top_x_limit,top_y_limit,u1(kkx,kky),umex,umey,vmex,vmey,ue(mx,my),ve(mx,my),un,vn&
,dux_b(my),dux_f(my),dvx_b(my),dvx_f(my),duy_b(mx),duy_f(mx),dvy_b(mx),dvy_f(mx)&
,dc1(mx,my),dc2(mx,my)
character(5)::int_method
kney=(ny-my)/2
knex=(nx-mx)/2
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)


call tuu(u1,u,nx,ny,kkx,kky)

call decy(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call tdu(duy,dccy,nx,ny,kkx,kky)

call decx(u,dux,nx,ny,dlanda)

call tuu(u1,v,nx,ny,kkx,kky)

call decy(u1,dccy,kkx,kky,dphi)
do i=1,kkx
dccy(i,ny)=(dccy(i,ny-1)+dccy(i,ny+1))/2.0d0
dccy(i,1)=(dccy(i,kky)+dccy(i,2))/2.0d0
enddo
call tdu(dvy,dccy,nx,ny,kkx,kky)
call decx(v,dvx,nx,ny,dlanda)
do i=1,nx
dux(i,1)=0.0d0
dvx(i,1)=0.0d0

dux(i,ny)=0.0d0
dvx(i,ny)=0.0d0

enddo

do j=1,my
a12=-top_y_limit+dble(j-1)*dphi
do i=1,mx
dunx(i,j)=dux(i+knex,j+kney)
dvnx(i,j)=dvx(i+knex,j+kney)
duny(i,j)=duy(i+knex,j+kney)
dvny(i,j)=dvy(i+knex,j+kney)

a13=-top_x_limit+dble(i-1)*dlanda
call yang_to_yin(a10,a11,a12,a13)
!if(dabs(dcos(a10))<=1.0d-8)then
!duex(i,j)=0.0d0
!dvex(i,j)=0.0d0
!duey(i,j)=0.0d0
!dvey(i,j)=0.0d0
!else
call select_int_method_glob(a10,a11,dux,nx,ny,umnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,duy,nx,ny,umny,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dvx,nx,ny,vmnx,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,dvy,nx,ny,vmny,dlanda,dphi,int_method)
call yangvector_to_yinvector(a10,a11,a12,a13,vn,un,ve(i,j),ue(i,j))
call Trans_vector_Deriv_Yin_to_Yang(a10,a11,a12,a13,un,vn,umnx,vmnx,umny,vmny,umex,vmex,umey,vmey)
duex(i,j)=umex
duey(i,j)=umey
dvex(i,j)=vmex
dvey(i,j)=vmey
!endif
enddo
enddo
do j=1,my
!a12=-top_y_limit+dble(j-1)*dphi
dux_b(j)=duex(1,j)
dux_f(j)=duex(mx,j)
dvx_b(j)=dvex(1,j)
dvx_f(j)=dvex(mx,j)

enddo
do i=1,mx

duy_b(i)=duey(i,1)
duy_f(i)=duey(i,my)
dvy_b(i)=dvey(i,1)
dvy_f(i)=dvey(i,my)
enddo
!call yang_to_yin(a10,a11,a12,a13)
!if(dabs(dcos(a10))<=1.0d-10)then
!duex(i,j)=(ue(i+1,j)-ue(i-1,j))/(2.0d0*dlanda)
!dvex(i,j)=(ve(i+1,j)-ve(i-1,j))/(2.0d0*dlanda)
!duey(i,j)=(ue(i,j+1)-ue(i,j-1))/(2.0d0*dphi)
!dvey(i,j)=(ve(i,j+1)-ve(i,j-1))/(2.0d0*dphi)
!endif
call dfx42(ue,dc1,mx,my,dlanda,dux_f)
call dbx42(ue,dc2,mx,my,dlanda,dux_b)
do j=1,my
do i=1,mx
duex(i,j)=(dc1(i,j)+dc2(i,j))/2.0d0
enddo
enddo
call dfx42(ve,dc1,mx,my,dlanda,dvx_f)
call dbx42(ve,dc2,mx,my,dlanda,dvx_b)
do j=1,my
do i=1,mx
dvex(i,j)=(dc1(i,j)+dc2(i,j))/2.0d0
enddo
enddo
call dby42(ue,dc1,mx,my,dphi,duy_b)
call dfy42(ue,dc2,mx,my,dphi,duy_f)
do j=1,my
do i=1,mx
duey(i,j)=(dc1(i,j)+dc2(i,j))/2.0d0
enddo
enddo
call dby42(ve,dc1,mx,my,dphi,dvy_b)
call dfy42(ve,dc2,mx,my,dphi,dvy_f)
do j=1,my
do i=1,mx
dvey(i,j)=(dc1(i,j)+dc2(i,j))/2.0d0
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine select_int_method_glob(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,int_method)
implicit none
integer::nx,ny

real*8::landa3,phi3,dphi,dlanda,hm
real*8::heb(nx,ny)
character(5)::int_method
if(int_method=='bilin')then
call bilinear_glob(phi3,landa3,heb,nx,ny,hm,dlanda,dphi)
elseif(int_method=='bicub')then
call bicubic_glob(phi3,landa3,heb,nx,ny,hm,dlanda,dphi)
endif
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine select_int_method(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y,int_method)
implicit none
integer::nx,ny,n_epsilon_x,n_epsilon_y

real*8::landa3,phi3,dphi,dlanda,hm
real*8::heb(nx,ny)
character(5)::int_method
if(int_method=='bilin')then
call bilinear(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y)
elseif(int_method=='bicub')then
call bicubic(phi3,landa3,heb,nx,ny,hm,dlanda,dphi,n_epsilon_x,n_epsilon_y)
endif
return
end subroutine 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine bilinear_glob(phi3,landa3,heb,nx,ny,hem,dlanda,dphi)
implicit none
integer::nx,ny
integer::k1,k2
real*8::pi,landa3,phi3,dphi,dlanda,hem,dl,dp,phi,landa,hem1,hem2
real*8::heb(nx,ny),top_x_limit,top_y_limit


pi=2.d0*dasin(1.0d0)
top_x_limit=pi
top_y_limit=(pi/2.0d0)


k1=int(((landa3+top_x_limit+1.0d-10)/dlanda))+1
k2=int(((phi3+top_y_limit+1.0d-10)/dphi))+1

landa=dble(k1-1)*dlanda-top_x_limit
phi=(dble(k2-1))*dphi-top_y_limit
dl=landa3-landa
dp=phi3-phi
if(k2==ny.and.k1==nx)then
call linear(heb(k1,k2),heb(1,k2),dlanda,dl,hem)
elseif(k2==ny)then
call linear(heb(k1,k2),heb(k1+1,k2),dlanda,dl,hem)
elseif(k1==nx)then
call linear(heb(k1,k2),heb(1,k2),dlanda,dl,hem1)
call linear(heb(k1,k2+1),heb(1,k2+1),dlanda,dl,hem2)
call linear(hem1,hem2,dphi,dp,hem)
else
call linear(heb(k1,k2),heb(k1+1,k2),dlanda,dl,hem1)
call linear(heb(k1,k2+1),heb(k1+1,k2+1),dlanda,dl,hem2)
call linear(hem1,hem2,dphi,dp,hem)
endif

return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  dbx42(e,dbx,m,k,dx,db)
integer::i,j,m,k
real*8::e(m,k),dbx(m,k),a,dx,db(k)
a=0.5d0-0.5d0/dsqrt(3.0d0)
do j=1,k
dbx(1,j)=db(j)
do i=2,m
dbx(i,j)=(((e(i,j)-e(i-1,j))/dx)-(a*dbx(i-1,j)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dfx42(e,dfx,m,k,dx,df)
integer::i,j,m,k
real*8::e(m,k),dfx(m,k),a,dx,df(k)
a=0.5d0-0.5d0/dsqrt(3.0d0)
do j=1,k
dfx(m,j)=df(j)
do i=m-1,1,-1
dfx(i,j)=(((e(i+1,j)-e(i,j))/dx)-(a*dfx(i+1,j)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dfy42(e,dfy,m,k,dy,df)
integer::i,j,m,k
real*8::e(m,k),dfy(m,k),a,dy,df(m)
a=0.5d0-0.5d0/dsqrt(3.0d0)
do i=1,m
dfy(i,k)=df(i)
do j=k-1,1,-1
dfy(i,j)=(((e(i,j+1)-e(i,j))/dy)-(a*dfy(i,j+1)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
subroutine  dby42(e,dby,m,k,dy,db)
integer::i,j,m,k
real*8::e(m,k),dby(m,k),a,dy,db(m)
a=0.5d0-0.5d0/dsqrt(3.0d0)
do i=1,m
dby(i,1)=db(i)
do j=2,k
dby(i,j)=(((e(i,j)-e(i,j-1))/dy)-(a*dby(i,j-1)))/(1.0d0-a)
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
subroutine  decompose_yang(hene,u,v,hee,ue,ve,nx,ny,mx,my,n_epsilon_x,n_epsilon_y,dphi,dlanda,int_method)
integer::i,j,nx,ny,mx,my,n_epsilon_x,n_epsilon_y
real*8::hene(nx,ny),u(nx,ny),v(nx,ny),hee(mx,my),ue(mx,my),ve(mx,my),dlanda,dphi,top_x_limit,top_y_limit,pi&
,a10,a11,a12,a13,hm,um,vm,ux,vy
character(5)::int_method
pi=2.0d0*dasin(1.0d0)
top_x_limit=(3.0d0*pi/4.0d0)+(dble(n_epsilon_x)*dlanda)
top_y_limit=(pi/4.0d0)+(dble(n_epsilon_y)*dphi)
do j=1,my
a12=-top_y_limit+dble(j-1)*dphi
do i=1,mx
a13=-top_x_limit+dble(i-1)*dlanda
call yang_to_yin(a10,a11,a12,a13)
if(a10>=-top_y_limit.and.a10<=top_y_limit)then
call select_int_method_glob(a10,a11,hene,nx,ny,hm,dlanda,dphi,int_method)
hee(i,j)=hm
call select_int_method_glob(a10,a11,u,nx,ny,um,dlanda,dphi,int_method)
call select_int_method_glob(a10,a11,v,nx,ny,vm,dlanda,dphi,int_method)
call yinvector_to_yangvector(a10,a11,a12,a13,vm,um,vy,ux)
ue(i,j)=ux
ve(i,j)=vy
else
hee(i,j)=hee(i,j)
ue(i,j)=ue(i,j)
ve(i,j)=ve(i,j)
endif
enddo
enddo
return
end subroutine
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
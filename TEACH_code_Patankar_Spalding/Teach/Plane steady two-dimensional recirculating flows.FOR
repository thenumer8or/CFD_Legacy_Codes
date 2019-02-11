c**********************************************************************
c                             teach - t
c**********************************************************************
c     A computer program for the calculation of plane or axisymmetric
c            steady two-dimensional recirculating flows                        
c**********************************************************************
c
c
c     subroutine main
      include'common.inc'
c
chapter  0  0  0  0  0  0  0  0  preliminaries  0  0  0  0  0  0  0  0
c
      write(*,*)'select method:(1=hybrid,2=power)'
      read(*,*)ihp
      if(ihp.eq.1) then
          open(1,file='hybrid.dat')
        else if(ihp.eq.2) then
          open(2,file='power.dat')
        else
          stop
      endif
      great=1.e30
      niter=0
      nswpu=3 
      nswpv=3
      nswpp=5
      nswpk=3
      nswpd=3
c
chapter  1  1  1  1  parameters and control indices  1  1  1  1  1  1
c
c-----grid 
      nim1=ni-1
      njm1=nj-1
      indcos=1
      ista=20
      istb=30
      jst=25
      istap1=ista+1
      istam1=ista-1
      istbp1=istb+1
      istbm1=istb-1
      jstp1=jst+1
      jstm1=jst-1
      wl=30
      w=1.5
      wr=30
      h=3
      h0=15
      s=1.2
      sumx=0.5+s**(ista-1)/2+(s-s**(ista-1))/(1-s)
      dxl=wl/sumx
      dx=dxl
      x(ista)=wl-dx/2
      do i=istam1,1,-1
         dx=dx*s
         x(i)=x(i+1)-dx
      end do
      dxm=w/float(istb-istap1)
      x(istap1)=x(ista)+dxl/2+dxm/2
      do i=istap1+1,istbm1
         x(i)=x(i-1)+dxm
      end do
      sumx=0.5+s**(ni-istb)/2+(s-s**(ni-istb))/(1-s)
      dxr=wr/sumx
      dx=dxr
      x(istb)=x(istbm1)+dxm/2+dxr/2
      do i=istbp1,ni
         dx=dx*s
         x(i)=x(i-1)+dx
      end do
      dyb=h/float(jst-2)
      y(1)=-0.5*dyb
      do j=2,jstm1
         y(j)=y(j-1)+dyb
      end do
      sumy=0.5+s**(nj-jst)/2+(s-s**(nj-jst))/(1-s)
      dyt=(h0-h)/sumy
      dy=dyt
      y(jst)=y(jstm1)+dyb/2+dyt/2
      do j=jstp1,nj
         dy=dy*s
         y(j)=y(j-1)+dy
      end do
c-----fiuid   properties
      densit=1.16
      viscos=15e-6
      prandt=0.70
c-----turbulence  constants
      cmu=0.09
      cd=1.00
      c1=1.44
      c2=1.92
      cappa=0.4187
      elog=9.793
      pred=1.3
c     pred=cappa*cappa/(c2-c1)/(cmu**0.5)
      prte=1.0
c-----boundary  values
      write(*,*)'enter velocity:'
      read(*,*)uin
      tein=0.0
      edin=0.0
c      factor=0.00001
c      tein=factor*uin**2
c      edin=factor*uin**3/h0
c-----pressure   calculation
      ipref=2
      jpref=2
c-----program   control and monitor
      write(*,*)'enter no. of iteration:'
      read(*,*) maxit
      sormax=0.001
c
chater  2  2  2  2  2  2  initial   operations  2  2  2  2  2  2  2  2
c
c-----calculate geometrical quantities and set variables to zero
      call init
c-----initialise variable  fields
      flowin=0.0
      do j=2,nj
          if (y(j).le.6) then
               u(2,j)=uin*(y(j)/6)**0.28
             else
                u(2,j)=uin
          end if
      end do
      do j=2,njm1
         te(1,j)=tein
         ed(1,j)=edin
         arden=0.5*(den(1,j)+den(2,j))*r(j)*sns(j)
         flowin=flowin+arden*u(2,j)
      end do
      do i=3,ista
         do j=2,nj
            u(i,j)=uin  
         end do
      end do
      do i=istap1,istb
         do j=jst,nj
            u(i,j)=uin*h0/(h0-h)   
         end do
      end do
      do i=istbp1,ni
         do j=2,nj
            u(i,j)=uin   
         end do
      end do
      do i=1,ni
         do j=2,nj
            if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 10
            te(i,j)=tein
            ed(i,j)=edin
10       end do
      end do
      do  i=2,nim1
          ypluss(i)=11.0
      end do
      do j=2,jst
         xplusw(j)=11.0
         xpluse(j)=11.0
      end do   
      call props
c-----calculate  source-normalisation factors
      snormm=flowin
      snormu=flowin*uin
c
chapter  3  3  3  3  3  3  iteration loop  3  3  3  3  3  3  3  3  3  3
c
  300 niter=niter+1 
      urfu=0.15
      urfv=0.15
      urfp=0.5
      urfe=0.95
      urfk=0.95
      urfvis=0.9
      write(*,*)niter
c-----update  main  dependent  variables
      call calcu     
      call calcv
      call calcp
      call calcte
      call calced
c-----update fluid properties
      call props
c-----normalise  residual-source  terms
      resoru=resoru/snormu
      resorv=resorv/snormu
      resorm=resorm/snormm
c-----termination  tests
      source=dmax1(resoru,resorv,resorm,resork)
      if (niter.ge.maxit) go to 302
      if (niter.lt.20.or.source.gt.sormax) go to 300
  302 continue
c
chapter  4  4  4  4  4  final operations and output  4  4  4  4  4  4
c
       write(*,*) source
       if(ihp.eq.1) then
            write(1,*)'variables="x","y","u","v","p","k","d","mu"'
            write(1,*)'zone  i=50  j=40 '
         else if(ihp.eq.2) then
            write(2,*)'variables="x","y","u","v","p","k","d","mu"'
            write(2,*)'zone  i=50  j=40 '
       end if     
c.....................................      
      do j=2,njm1
         do i=2,nim1
            u(i,j)=(u(i,j)+u(i+1,j))/2
            v(i,j)=(v(i,j)+v(i,j+1))/2
            vis(i,j)=den(i,j)*te(i,j)**2*cmu/ed(i,j)
         end do
      end do
      do j=2,njm1
         do  i=2,nim1
              if(ihp.eq.1) then
                  write(1,7) x(i),y(j),u(i,j),v(i,j),p(i,j),
     &                        te(i,j),ed(i,j),vis(i,j)
                else if(ihp.eq.2) then
                  write(2,7) x(i),y(j),u(i,j),v(i,j),p(i,j),       
     &                       te(i,j),ed(i,j),vis(i,j)
              end if                 
         end do
      end do   
7     format(8(5x,1pe11.4))
      stop
      end
c***********************************************************************
      subroutine init
      include'common.inc'
c
chapter  1  1  1  1  1  geometrical quantities  1  1  1  1  1  1  1  1
c
      do   100  j=1,nj
      r(j)=y(j)
  100 if (indcos.eq.1) r(j)=1.0
      dxpw(1)=0.0
      dxep(ni)=0.0
      do  101  i=1,nim1
      dxep(i)=x(i+1)-x(i)
 101  dxpw(i+1)=dxep(i)
      dyps(1)=0.0
      dynp(nj)=0.0
      do  102  j=1,njm1
      dynp(j)=y(j+1)-y(j)
 102  dyps(j+1)=dynp(j)
      sew(1)=0.0
      sew(ni)=0.0
      do  103  i=2,nim1
 103  sew(i)=0.5*(dxep(i)+dxpw(i))
      sns(1)=0.0
      sns(nj)=0.0
      do  104  j=2,njm1
 104  sns(j)=0.5*(dynp(j)+dyps(j))
      xu(1)=0.0
      do  105 i=2,ni
  105 xu(i)=0.5*(x(i)+x(i-1))
      dxpwu(1)=0.0
      dxpwu(2)=0.0
      dxepu(1)=0.0
      dxepu(ni)=0.0
      do  106  i=2,nim1
      dxepu(i)=xu(i+1)-xu(i)
 106  dxpwu(i+1)=dxepu(i)
      sewu(1)=0.0
      sewu(2)=0.0
      do  107  i=3,nim1
 107  sewu(i)=0.5*(dxepu(i)+dxpwu(i))
      yv(1)=0.0
      rv(1)=0.0
      do  108  j=2,nj
      rv(j)=0.5*(r(j)+r(j-1))
      rcv(j)=0.5*(rv(j)+rv(j-1))
 108  yv(j)=0.5*(y(j)+y(j-1))
      dypsv(1)=0.0
      dypsv(2)=0.0
      dynpv(nj)=0.0
      do  109  j=2,njm1
      dynpv(j)=yv(j+1)-yv(j)
 109  dypsv(j+1)=dynpv(j)
      snsv(1)=0.0
      snsv(2)=0.0
      snsv(nj)=0.0
      do  110  j=3,njm1
  110 snsv(j)=0.5*(dynpv(j)+dypsv(j))
c
chapter  2  2  2  2  2  2  set variable to zero  2  2  2  2  2  2  2
c
      do  200  i=1,ni
      do  200  j=1,nj
         u(i,j)=0.0
         v(i,j)=0.0
         p(i,j)=0.0
         pp(i,j)=0.0
         te(i,j)=1.0e-10
         ed(i,j)=0.0001
         den(i,j)=densit
         vis(i,j)=viscos
         du(i,j)=0.0
         dv(i,j)=0.0
         su(i,j)=0.0
         sp(i,j)=0.0
  200 continue
      return
      end
c*********************************************************************
      subroutine props
      include'common.inc'
c
chapter  1  1  1  1  1   viscosity  1  1  1  1  1  1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 101 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 101
      visold=vis(i,j)
      if (ed(i,j).eq.0.) go to 102
      vis(i,j)=den(i,j)*te(i,j)**2*cmu/ed(i,j)+viscos
      go to 103
  102 vis(i,j)=viscos
c-----under-relax viscosity
  103 vis(i,j)=urfvis*vis(i,j)+(1.-urfvis)*visold
  101 continue
  100 continue
      return
      end
c**********************************************************************
      subroutine calcu
      include'common.inc'
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      do 100 i=3,nim1
      do 101 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.le.istb)) goto 101
c-----compute areas and volume
      arean=rv(j+1)*sewu(i)
      areas=rv(j)*sewu(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sewu(i)*sns(j)
c-----calculate convection coefficients
      gn=0.5*(den(i,j+1)+den(i,j))*v(i,j+1)
      gnw=0.5*(den(i-1,j)+den(i-1,j+1))*v(i-1,j+1)
      gs=0.5*(den(i,j-1)+den(i,j))*v(i,j)
      gsw=0.5*(den(i-1,j)+den(i-1,j-1))*v(i-1,j)
      ge=0.5*(den(i+1,j)+den(i,j))*u(i+1,j)
      gp=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      gw=0.5*(den(i-1,j)+den(i-2,j))*u(i-1,j)
      cn=0.5*(gn+gnw)*arean
      cs=0.5*(gs+gsw)*areas
      ce=0.5*(ge+gp)*areaew
      cw=0.5*(gp+gw)*areaew
c-----calculate diffusion coefficients
      visn=0.25*(vis(i,j)+vis(i,j+1)+vis(i-1,j)+vis(i-1,j+1))
      viss=0.25*(vis(i,j)+vis(i,j-1)+vis(i-1,j)+vis(i-1,j-1))
      dn=visn*arean/dynp(j)
      ds=viss*areas/dyps(j)
      de=vis(i,j)*areaew/dxepu(i)
      dw=vis(i-1,j)*areaew/dxpwu(i)
c-----calculate coefficients of source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0,smp)
      cpo=cp
c-----assemble main coefficients
      if(ihp.eq.1) then
          an(i,j)=dn*dmax1(0,1-0.5*abs(cn/dn))+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,1-0.5*abs(cs/ds))+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,1-0.5*abs(ce/de))+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,1-0.5*abs(cw/dw))+dmax1(cw,0)
        else if(ihp.eq.2) then
          an(i,j)=dn*dmax1(0,(1-0.1*abs(cn/dn))**5)+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,(1-0.1*abs(cs/ds))**5)+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,(1-0.1*abs(ce/de))**5)+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,(1-0.1*abs(cw/dw))**5)+dmax1(cw,0)
      end if
      du(i,j)=areaew
      su(i,j)=cpo*u(i,j)+du(i,j)*(p(i-1,j)-p(i,j))
      sp(i,j)=-cp
  101 continue
  100 continue   

c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modu
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resoru=0.0
      do 300 i=3,nim1
      do 301 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.le.istb)) goto 301
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j) 
      du(i,j)=du(i,j)/ap(i,j)
      resor=an(i,j)*u(i,j+1)+as(i,j)*u(i,j-1)+ae(i,j)*u(i+1,j)
     1    +aw(i,j)*u(i-1,j)-ap(i,j)*u(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resoru=resoru+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfu
      su(i,j)=su(i,j)+(1.-urfu)*ap(i,j)*u(i,j)
      du(i,j)=du(i,j)*urfu
 301  continue
 300  continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
      index=1
      do  400  n=1,nswpu
 400  call lisolv(index,3,2,u)
      return     
      end
c**********************************************************************
      subroutine calcv
      include'common.inc'
c
chapter  1  1  1  1  1  assembly of coefficients  1  1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 101 j=3,njm1
      if(j.le.jst.and.(i.gt.ista.and.i.lt.istb)) goto 101
c-----compute areas and volume
      arean=rcv(j+1)*sew(i)
      areas=rcv(j)*sew(i)
      areaew=rv(j)*snsv(j)
      vol=rv(j)*sew(i)*snsv(j)
c-----calculate convection coefficients
      gn=0.5*(den(i,j+1)+den(i,j))*v(i,j+1)
      gp=0.5*(den(i,j)+den(i,j-1))*v(i,j)
      gs=0.5*(den(i,j-1)+den(i,j-2))*v(i,j-1)
      ge=0.5*(den(i+1,j)+den(i,j))*u(i+1,j)
      gse=0.5*(den(i,j-1)+den(i+1,j-1))*u(i+1,j-1)
      gw=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      gsw=0.5*(den(i,j-1)+den(i-1,j-1))*u(i,j-1)
      cn=0.5*(gn+gp)*arean
      cs=0.5*(gs+gp)*areas
      ce=0.5*(ge+gse)*areaew
      cw=0.5*(gw+gsw)*areaew
c-----calculate diffusion coefficients
      vise=0.25*(vis(i,j)+vis(i+1,j)+vis(i,j-1)+vis(i+1,j-1))
      visw=0.25*(vis(i,j)+vis(i-1,j)+vis(i,j-1)+vis(i-1,j-1))
      dn=vis(i,j)*arean/dynpv(j)
      ds=vis(i,j-1)*areas/dypsv(j)
      de=vise*areaew/dxep(i)
      dw=visw*areaew/dxpw(i)
c-----calculate coefficients of source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0,smp)
      cpo=cp
c-----assemble main coefficients
      if(ihp.eq.1) then
          an(i,j)=dn*dmax1(0,1-0.5*abs(cn/dn))+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,1-0.5*abs(cs/ds))+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,1-0.5*abs(ce/de))+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,1-0.5*abs(cw/dw))+dmax1(cw,0)
        else if(ihp.eq.2) then
          an(i,j)=dn*dmax1(0,(1-0.1*abs(cn/dn))**5)+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,(1-0.1*abs(cs/ds))**5)+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,(1-0.1*abs(ce/de))**5)+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,(1-0.1*abs(cw/dw))**5)+dmax1(cw,0)
      end if
      dv(i,j)=0.5*(arean+areas)
      su(i,j)=cpo*v(i,j)+dv(i,j)*(p(i,j-1)-p(i,j))
      sp(i,j)=-cp
      if (indcos.eq.2) sp(i,j)=sp(i,j)-2*vis(i,j)*vol/rv(j)**2
 101  continue
 100  continue
      arean=rcv(3)*sew(2)
      areas=rcv(2)*sew(2)
      dv(2,2)=0.5*(arean+areas)
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modv
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resorv=0.0
      do 300 i=2,nim1
      do 301 j=3,njm1
      if(j.le.jst.and.(i.gt.ista.and.i.lt.istb)) goto 301
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      dv(i,j)=dv(i,j)/ap(i,j)
      resor=an(i,j)*v(i,j+1)+as(i,j)*v(i,j-1)+ae(i,j)*v(i+1,j)
     1    +aw(i,j)*v(i-1,j)-ap(i,j)*v(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resorv=resorv+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfv
      su(i,j)=su(i,j)+(1.-urfv)*ap(i,j)*v(i,j)
      dv(i,j)=dv(i,j)*urfv
  301 continue
  300 continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
      index=2
      do 400 n=1,nswpv
 400  call lisolv(index,2,3,v)
      return
      end
c***********************************************************************
        subroutine calcp
        include'common.inc'
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c     resorm=0.0
      do 100 i=2,nim1
      do 101 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 101
c-----compute areas and volume
      arean=rv(j+1)*sew(i)
      areas=rv(j)*sew(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sns(j)*sew(i)
c-----calculate coefficients
      denn=0.5*(den(i,j)+den(i,j+1))
      dens=0.5*(den(i,j)+den(i,j-1))
      dene=0.5*(den(i,j)+den(i+1,j))
      denw=0.5*(den(i,j)+den(i-1,j))
      an(i,j)=denn*arean*dv(i,j+1)
      as(i,j)=dens*areas*dv(i,j)
      ae(i,j)=dene*areaew*du(i+1,j)
      aw(i,j)=denw*areaew*du(i,j)
c-----calculate source term
      cn=denn*v(i,j+1)*arean
      cs=dens*v(i,j)*areas
      ce=dene*u(i+1,j)*areaew
      cw=denw*u(i,j)*areaew
      smp=cn-cs+ce-cw
      sp(i,j)=0.0
      su(i,j)=-smp
c-----compute sum of absolute mass sources
      resorm=resorm+abs(smp)
 101  continue
 100  continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modp
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      do 300 i=2,nim1
      do 301 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 301
  301 ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
  300 continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
chapter  5  5  5  5  correct velocities and pressure  5  5  5  5  5  5
c      
      index=3
      do 400 n=1,nswpp
  400 call lisolv(index,2,2,pp)
c
c-----velocities
      do 500 i=2,nim1
      do 501 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 501
      if(i.eq.2) then
          du(i,j)=0.0
        else if(i.eq.istb.and.j.lt.jst) then
          du(i,j)=0.0
      end if
      u(i,j)=u(i,j)+du(i,j)*(pp(i-1,j)-pp(i,j))
      if(j.eq.2) then
           dv(i,j)=0.0
        else if(j.eq.jst.and.(i.gt.ista.and.i.lt.istb)) then
           dv(i,j)=0.0
      end if
         v(i,j)=v(i,j)+dv(i,j)*(pp(i,j-1)-pp(i,j))
  501 continue
  500 continue
c-----pressures (with provision for under-relaxation)
      ppref=pp(ipref,jpref)
      do 502 i=2,nim1
      do 503 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 503
      p(i,j)=p(i,j)+urfp*(pp(i,j)-ppref)
      pp(i,j)=0.0
  503 continue
  502 continue
      return
      end
c***********************************************************************
      subroutine calcte
      include'common.inc'
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      prte=1.0
      do 100 i=2,nim1
      do 101 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 101
c-----compute areas and volume
      arean=rv(j+1)*sew(i)
      areas=rv(j)*sew(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sns(j)*sew(i)
c-----calculate convection coefficients
      gn=0.5*(den(i,j)+den(i,j+1))*v(i,j+1)
      gs=0.5*(den(i,j)+den(i,j-1))*v(i,j)
      ge=0.5*(den(i,j)+den(i+1,j))*u(i+1,j)
      gw=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      cn=gn*arean
      cs=gs*areas
      ce=ge*areaew
      cw=gw*areaew
c-----calculate diffusion coefficients
      gamn=0.5*(vis(i,j)+vis(i,j+1))/prte
      gams=0.5*(vis(i,j)+vis(i,j-1))/prte
      game=0.5*(vis(i,j)+vis(i+1,j))/prte
      gamw=0.5*(vis(i,j)+vis(i-1,j))/prte
      dn=gamn*arean/dynp(j)
      ds=gams*areas/dyps(j)
      de=game*areaew/dxep(i)
      dw=gamw*areaew/dxpw(i)
c-----source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0,smp)
      cpo=cp
      dudx=(u(i+1,j)-u(i,j))/sew(i)
      dvdy=(v(i,j+1)-v(i,j))/sns(j)
      dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+u(i+1,j)+
     1      u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
        dvdx=((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4.-(v(i,j)+v(i,j+1)+
     1      v(i-1,j)+v(i-1,j+1))/4.)/sew(i)
        gen(i,j)=(2.*(dudx**2+dvdy**2)+(dudy+dvdx)**2)*vis(i,j)
c-----assemble main coefficients
      if(ihp.eq.1) then
          an(i,j)=dn*dmax1(0,1-0.5*abs(cn/dn))+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,1-0.5*abs(cs/ds))+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,1-0.5*abs(ce/de))+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,1-0.5*abs(cw/dw))+dmax1(cw,0)
        else if(ihp.eq.2) then
          an(i,j)=dn*dmax1(0,(1-0.1*abs(cn/dn))**5)+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,(1-0.1*abs(cs/ds))**5)+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,(1-0.1*abs(ce/de))**5)+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,(1-0.1*abs(cw/dw))**5)+dmax1(cw,0)
      end if
      su(i,j)=cpo*te(i,j)
      sukd(i,j)=su(i,j)
      su(i,j)=su(i,j)+gen(i,j)*vol
      sp(i,j)=-cp
      spkd(i,j)=sp(i,j)
      sp(i,j)=sp(i,j)-cd*cmu*den(i,j)**2*te(i,j)*vol/vis(i,j)
  101 continue
  100 continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call modte
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resork=0.0
      do 300 i=2,nim1
      do 301 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 301
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      resor=an(i,j)*te(i,j+1)+as(i,j)*te(i,j-1)+ae(i,j)*te(i+1,j)
     1    +aw(i,j)*te(i-1,j)-ap(i,j)*te(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resork=resork+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfk
      su(i,j)=su(i,j)+(1.-urfk)*ap(i,j)*te(i,j)
 301  continue
 300  continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c      
      index=3
      do 400 n=1,nswpk
  400 call lisolv(index,2,2,te)
      return
      end
c***********************************************************************
      subroutine calced
      include'common.inc'
c
chapter  1  1  1  1  1  assembly of coefficients   1  1  1  1  1  1  1
c
      do 100 i=2,nim1
      do 101 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 101
c-----compute areas and volume
      arean=rv(j+1)*sew(i)
      areas=rv(j)*sew(i)
      areaew=r(j)*sns(j)
      vol=r(j)*sns(j)*sew(i)
c-----calculate convection coefficients
      gn=0.5*(den(i,j)+den(i,j+1))*v(i,j+1)
      gs=0.5*(den(i,j)+den(i,j-1))*v(i,j)
      ge=0.5*(den(i,j)+den(i+1,j))*u(i+1,j)
      gw=0.5*(den(i,j)+den(i-1,j))*u(i,j)
      cn=gn*arean
      cs=gs*areas
      ce=ge*areaew
      cw=gw*areaew
c-----calculate diffusion coefficients
      gamn=0.5*(vis(i,j)+vis(i,j+1))/pred
      gams=0.5*(vis(i,j)+vis(i,j-1))/pred
      game=0.5*(vis(i,j)+vis(i+1,j))/pred
      gamw=0.5*(vis(i,j)+vis(i-1,j))/pred
      dn=gamn*arean/dynp(j)
      ds=gams*areas/dyps(j)
      de=game*areaew/dxep(i)
      dw=gamw*areaew/dxpw(i)
c-----source terms
      smp=cn-cs+ce-cw
      cp=dmax1(0.0,smp)
      cpo=cp
c-----assemble main coefficients
      if(ihp.eq.1) then
          an(i,j)=dn*dmax1(0,1-0.5*abs(cn/dn))+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,1-0.5*abs(cs/ds))+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,1-0.5*abs(ce/de))+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,1-0.5*abs(cw/dw))+dmax1(cw,0)
        else if(ihp.eq.2) then
          an(i,j)=dn*dmax1(0,(1-0.1*abs(cn/dn))**5)+dmax1(-cn,0)
          as(i,j)=ds*dmax1(0,(1-0.1*abs(cs/ds))**5)+dmax1(cs,0)
          ae(i,j)=de*dmax1(0,(1-0.1*abs(ce/de))**5)+dmax1(-ce,0)
          aw(i,j)=dw*dmax1(0,(1-0.1*abs(cw/dw))**5)+dmax1(cw,0)
      end if
      su(i,j)=cpo*ed(i,j)
      sukd(i,j)=su(i,j)
      su(i,j)=su(i,j)+c1*cmu*gen(i,j)*vol*den(i,j)*te(i,j)/vis(i,j)
      sp(i,j)=-cp
      spkd(i,j)=sp(i,j)
      sp(i,j)=sp(i,j)-c2*den(i,j)*ed(i,j)*vol/te(i,j)
 101  continue
 100  continue
c
chapter  2  2  2  2  2  2  2  problem modifications  2  2  2  2  2  2  2
c
      call moded
c
chapter  3  final coeff. assembly and residual source calculations  3  3
c
      resore=0.0
      do 300 i=2,nim1
      do 301 j=2,njm1
      if(j.lt.jst.and.(i.gt.ista.and.i.lt.istb)) goto 301
      ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
      resor=an(i,j)*ed(i,j+1)+as(i,j)*ed(i,j-1)+ae(i,j)*ed(i+1,j)
     1    +aw(i,j)*ed(i-1,j)-ap(i,j)*ed(i,j)+su(i,j)
      vol=r(j)*sew(i)*sns(j)
      sorvol=great*vol
      if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
      resore=resore+abs(resor)
c-----under-relaxation
      ap(i,j)=ap(i,j)/urfe
      su(i,j)=su(i,j)+(1.-urfe)*ap(i,j)*ed(i,j)
 301  continue
 300  continue
c
chapter  4  4  4  4  solution of difference equations  4  4  4  4  4  4
c
      index=3
      do 400 n=1,nswpd
 400  call lisolv(index,2,2,ed)
      return
      end
c***********************************************************************
        subroutine lisolv(index,istart,jstart,phi)
        include'common.inc'
      dimension phi(ni,nj),a(nj),b(nj),c(nj),d(nj)
      jm1=jstart-1
      a(jm1)=0
c-----commence w-e sweep
      do 100 i=istart,ista
      c(jm1)=phi(i,jm1)
c-----commence s-n traverse
      do 101 j=jstart,njm1
c-----assemble toma coefficients
      a(j)=an(i,j)
      b(j)=as(i,j)
      c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
      d(j)=ap(i,j)
c-----calculate coefficients of recorrence formula
      term=1./(d(j)-b(j)*a(j-1))
      a(j)=a(j)*term
 101  c(j)=(c(j)+b(j)*c(j-1))*term
c-----obtain new phi's
      do 102 jj=jstart,njm1
      j=nj+jm1-jj
 102  phi(i,j)=a(j)*phi(i,j+1)+c(j)
 100  continue
      if(index.eq.1) then
         iend=istb
         jbgn=jst
       else if(index.eq.2) then
         iend=istbm1
         jbgn=jstp1
       else
         iend=istbm1
         jbgn=jst
      end if
      jbgnm1=jbgn-1
      a(jbgnm1)=0
c-----commence w-e sweep
      do 200 i=istap1,iend
      c(jbgnm1)=phi(i,jbgnm1)
c-----commence s-n traverse
      do 201 j=jbgn,njm1
c-----assemble toma coefficients
      a(j)=an(i,j)
      b(j)=as(i,j)
      c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
      d(j)=ap(i,j)
c-----calculate coefficients of recorrence formula
      term=1./(d(j)-b(j)*a(j-1))
      a(j)=a(j)*term
 201  c(j)=(c(j)+b(j)*c(j-1))*term
c-----obtain new phi's
      do 202 jj=jbgn,njm1
      j=nj+jbgnm1-jj
 202  phi(i,j)=a(j)*phi(i,j+1)+c(j)
 200  continue
      a(jm1)=0.0
c-----commence w-e sweep
      do 300 i=iend+1,nim1
      c(jm1)=phi(i,jm1)
c-----commence s-n traverse
      do 301 j=jstart,njm1
c-----assemble toma coefficients
      a(j)=an(i,j)
      b(j)=as(i,j)
      c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
      d(j)=ap(i,j)
c-----calculate coefficients of recorrence formula
      term=1./(d(j)-b(j)*a(j-1))
      a(j)=a(j)*term
 301  c(j)=(c(j)+b(j)*c(j-1))*term
c-----obtain new phi's
      do 302 jj=jstart,njm1
      j=nj+jm1-jj
 302  phi(i,j)=a(j)*phi(i,j+1)+c(j)
 300  continue
      return
      end
c*************************************************************************
      subroutine promod
      include'common.inc'
c
chapter  1  1  1  1  1  geometrical quantities  1  1  1  1  1  1  1  1
c
      entry modpro
c-----no modifications for this problem
      return
c
chpter  2  2  2  2  2  2  2 u  momentum  2  2  2  2  2  2  2  2  2  2  2
c
      entry modu
c-----top boundry
      do i=3,nim1 
         an(i,njm1)=0.0
      end do
c-----east wall
      do j=2,jstm1
         ae(ista,j)=0.0
      end do
c-----west wall
      do j=2,jstm1
         aw(istbp1,j)=0.0
      end do
c-----south wall left
      cdterm=cmu**0.25
      yp=y(2)
      j=2
      do i=3,ista
         sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
         denu=0.5*(den(i,j)+den(i-1,j))
         yplusa=0.5*(ypluss(i)+ypluss(i-1))
         if (yplusa.le.11.63)  go to 311
         tmult=denu*cdterm*sqrtk*cappa/alog(elog*yplusa)
         go to 312
 311     tmult=viscos/yp
 312     ueff=u(i,j)
         if(u(i,j) .lt.0.) ueff=-ueff
         taus(i)=-tmult*ueff
         sp(i,j)=sp(i,j)-tmult*sewu(i)
         as(i,j)=0.0
      end do
      taus(2)=taus(3)
c-----south wall middle
      cdterm=cmu**0.25
      yp=y(jst)-yv(jst)
      j=jst
      do i=istap1,istb
         sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
         denu=0.5*(den(i,j)+den(i-1,j))
         if(i.eq.istap1) then
              yplusa=ypluss(i)
           else if(i.eq.istb) then
              yplusa=ypluss(i-1)
           else
              yplusa=0.5*(ypluss(i)+ypluss(i-1))
         end if
         if (yplusa.le.11.63)  go to 411
         tmult=denu*cdterm*sqrtk*cappa/alog(elog*yplusa)
         go to 412
 411     tmult=viscos/yp
 412     ueff=u(i,j)
         if(u(i,j) .lt.0.) ueff=-ueff
         taus(i)=-tmult*ueff
         sp(i,j)=sp(i,j)-tmult*sewu(i)
         as(i,j)=0.0
      end do
c-----south wall right
      cdterm=cmu**0.25
      yp=y(2)
      j=2
      do i=istbp1,nim1
         sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
         denu=0.5*(den(i,j)+den(i-1,j))
         yplusa=0.5*(ypluss(i)+ypluss(i-1))
         if (yplusa.le.11.63)  go to 511
         tmult=denu*cdterm*sqrtk*cappa/alog(elog*yplusa)
         go to 512
 511     tmult=viscos/yp
 512     ueff=u(i,j)
         if(u(i,j) .lt.0.) ueff=-ueff
         taus(i)=-tmult*ueff
         sp(i,j)=sp(i,j)-tmult*sewu(i)
         as(i,j)=0.0
      end do
      taus(istb)=taus(istbp1)
c-----outlet
      ardent=0.0
      flow=0.0
      do 204 j=2,njm1
      arden=0.5*(den(nim1,j)+den(nim1-1,j))*r(j)*sns(j)
      ardent=ardent+arden
 204  flow=flow+arden*u(nim1,j)
      uinc=(flowin-flow)/ardent
      do 205 j=2,njm1
 205  u(ni,j)=u(nim1,j)+uinc
      return
c
chapter  3  3  3  3  3  3  3  3  v  momentum  3  3  3  3  3  3  3  3  3
c
      entry modv
c-----top boundry
      do i=2,nim1
         an(i,njm1)=0.0
      end do
c-----south wall left
      do i=2,ista
         as(i,3)=0.0
      end do
c-----south wall middle
      do i=istap1,istbm1
         as(i,jstp1)=0.0
      end do
c-----south wall right
      do i=istb,nim1
         as(i,3)=0.0
      end do
c-----west wall
      cdterm=cmu**0.25
      xp=x(istb)-xu(istb)
      i=istb
      do j=3,jst
         sqrtk=sqrt(0.5*(te(i,j)+te(i,j-1)))
         denv=0.5*(den(i,j)+den(i,j-1))
         if(j.eq.jst) then
             xplusa=xplusw(j-1)  
           else
             xplusa=0.5*(xplusw(j)+xplusw(j-1))
         end if
         if (xplusa.le.11.63)  go to 611
         tmult=denv*cdterm*sqrtk*cappa/alog(elog*xplusa)
         go to 612
 611     tmult=viscos/xp
 612     veff=v(i,j)
         if(v(i,j).lt.0.) veff=-veff
         tauw(j)=-tmult*veff
         sp(i,j)=sp(i,j)-tmult*snsv(j)*rv(j)
         aw(i,j)=0.0
      end do
      tauw(2)=tauw(3)
c.....east wall
      cdterm=cmu**0.25
      xp=xu(istap1)-x(ista)
      i=ista
      do j=3,jst
         sqrtk=sqrt(0.5*(te(i,j)+te(i,j-1)))
         denv=0.5*(den(i,j)+den(i,j-1))
         if(j.eq.jst) then
             xplusa=xpluse(j-1)  
           else
             xplusa=0.5*(xpluse(j)+xpluse(j-1))
         end if
         if (xplusa.le.11.63)  go to 711
         tmult=denv*cdterm*sqrtk*cappa/alog(elog*xplusa)
         go to 712
 711     tmult=viscos/xp
 712     veff=v(i,j)
         if(v(i,j).lt.0.) veff=-veff
         taue(j)=-tmult*veff
         sp(i,j)=sp(i,j)-tmult*snsv(j)*rv(j)
         ae(i,j)=0.0
      end do
      taue(2)=taue(3)
      return
c
chapter  4  4  4  4  4  4  4  pressure correction  4  4  4  4  4  4  4
c
      entry modp
      return
c
chapter  5  5  5  5  5  5  5  5  thermal energy  5  5  5  5  5  5  5  5
c
      entry modt
c-----no modifications for this problem
      return
c
chapter  6  6  6  6  6  6  6  turbulent kinetic energy   6  6  6  6  6
c
      entry modte
c-----top boundry
      j=njm1
      do i=2,nim1
         dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4-
     &         (u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4)/sns(j)
         vol=r(j)*sns(j)*sew(i)
         gen(i,j)=gen(i,j)-vis(i,j)*dudy**2
         su(i,j)=sukd(i,j)+gen(i,j)*vol
         an(i,njm1)=0.0
      end do
c-----south wall left
      cdterm=cmu**0.25
      yp=y(2)
      j=2
      do i=2,ista
         uavg=abs((u(i,j)+u(i+1,j))/2.)
         ueff=uavg
         if((u(i,j)+u(i+1,j))/2. .lt.0.) ueff=-ueff
         denu=den(i,j)
         sqrtk=sqrt(te(i,j))
         vol=r(j)*sns(j)*sew(i)
         if(i.eq.ista) then
             tau=taus(i)
           else
             tau=(taus(i)+taus(i+1))/2
         end if
         gencou=abs(tau*ueff)/yp
         ypluss(i)=denu*sqrtk*cdterm*yp/viscos
         dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+
     1         u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
         genres=gen(i,j)-vis(i,j)*dudy**2
         gen(i,j)=genres+gencou
         if (ypluss(i).le.11.63)  go to 121
         diterm=den(i,j)*(cmu**0.75)*sqrtk*alog(elog*ypluss(i))/
     1          (cappa*yp)
         go to 122
 121     continue
         diterm=den(i,j)*(cmu**0.75)*sqrtk*ypluss(i)/yp
 122     continue
         su(i,j)=gen(i,j)*vol+sukd(i,j)
         sp(i,j)=-diterm*vol+spkd(i,j)
         as(i,j)=0.0
      end do
c-----south wall middle
      cdterm=cmu**0.25
      yp=y(jst)-yv(jst)
      j=jst
      do i=istap1,istbm1
         uavg=abs((u(i,j)+u(i+1,j))/2.)
         ueff=uavg
         if((u(i,j)+u(i+1,j))/2. .lt.0.) ueff=-ueff
         denu=den(i,j)
         sqrtk=sqrt(te(i,j))
         vol=r(j)*sns(j)*sew(i)
         if(i.eq.istap1) then
             tau=taus(i+1)
           else if(i.eq.istbm1) then
             tau=taus(i)
           else
             tau=(taus(i)+taus(i+1))/2
         end if
         gencou=abs(tau*ueff)/yp
         ypluss(i)=denu*sqrtk*cdterm*yp/viscos
         dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+
     1         u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
         genres=gen(i,j)-vis(i,j)*(dudy**2)
         gen(i,j)=genres+gencou
         if (ypluss(i).le.11.63)  go to 221
         diterm=den(i,j)*(cmu**0.75)*sqrtk*alog(elog*ypluss(i))/
     1          (cappa*yp)
         go to 222
 221     continue
         diterm=den(i,j)*(cmu**0.75)*sqrtk*ypluss(i)/yp
 222     continue
         su(i,j)=gen(i,j)*vol+sukd(i,j)
         sp(i,j)=-diterm*vol+spkd(i,j)
         as(i,j)=0.0
      end do
c-----south wall right
      cdterm=cmu**0.25
      yp=y(2)
      j=2
      do i=istb,nim1
         uavg=abs((u(i,j)+u(i+1,j))/2.)
         ueff=uavg
         if((u(i,j)+u(i+1,j))/2. .lt.0.) ueff=-ueff
         denu=den(i,j)
         sqrtk=sqrt(te(i,j))
         vol=r(j)*sns(j)*sew(i)
         if(i.eq.istb) then
             tau=taus(i+1)
           else if(i.eq.nim1) then
             tau=taus(i)
             tau=(taus(i)+taus(i+1))/2
         end if
         gencou=abs(tau*ueff)/yp
         ypluss(i)=denu*sqrtk*cdterm*yp/viscos
         dudy=((u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))/4.-(u(i,j)+
     1         u(i+1,j)+u(i,j-1)+u(i+1,j-1))/4.)/sns(j)
         genres=gen(i,j)-vis(i,j)*(dudy**2)
         gen(i,j)=genres+gencou
         if (ypluss(i).le.11.63)  go to 321
         diterm=den(i,j)*(cmu**0.75)*sqrtk*alog(elog*ypluss(i))/
     1          (cappa*yp)
         go to 322
 321     continue
         diterm=den(i,j)*(cmu**0.75)*sqrtk*ypluss(i)/yp
 322     continue
         su(i,j)=gen(i,j)*vol+sukd(i,j)
         sp(i,j)=-diterm*vol+spkd(i,j)
         as(i,j)=0.0
      end do
c-----west wall
      xp=x(istb)-xu(istb)
      i=istb
      do j=2,jstm1
         vavg=abs((v(i,j)+v(i,j+1)))/2.
         veff=vavg
         if((v(i,j)+v(i,j+1))/2 .lt.0.) veff=-veff
         denv=den(i,j)
         sqrtk=sqrt(te(i,j))
         vol=r(j)*sns(j)*sew(i)
         xplusw(j)=denv*sqrtk*cdterm*xp/viscos
         if(j.eq.jstam1) then
             tau=tauw(j)
           else
             tau=(tauw(j)+tauw(j+1))/2
         end if
         gencou=abs(tau*veff)/xp
         dvdx=((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4.-
     1         (v(i,j)+v(i,j+1)+v(i-1,j)+v(i-1,j+1))/4.)/sew(i)
         genres=gen(i,j)-vis(i,j)*(dvdx**2)
         gen(i,j)=genres+gencou
         if (xplusw(j).le.11.63) go to 131
         diterm=den(i,j)*(cmu**0.75)*sqrtk*alog(elog*xplusw(j))/
     1    (cappa*xp)
         go to 132
 131     continue
         diterm=den(i,j)*(cmu**0.75)*sqrtk*xplusw(j)/xp
 132     continue
         su(i,j)=sukd(i,j)+gen(i,j)*vol
         sp(i,j)=spkd(i,j)-diterm*vol
         aw(i,j)=0.0
      end do
c-----east wall
      xp=xu(istap1)-x(ista)
      i=ista
      do j=2,jstm1
         vavg=abs((v(i,j)+v(i,j+1)))/2.
         veff=vavg
         if((v(i,j)+v(i,j+1))/2 .lt.0.) veff=-veff
         denv=den(i,j)
         sqrtk=sqrt(te(i,j))
         vol=r(j)*sns(j)*sew(i)
         xpluse(j)=denv*sqrtk*cdterm*xp/viscos
         if(j.eq.jstm1) then
             tau=taue(j)
           else
             tau=(taue(j)+taue(j+1))/2
         end if
         gencou=abs(tau*veff)/xp
         dvdx=((v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1))/4.-
     1         (v(i,j)+v(i,j+1)+v(i-1,j)+v(i-1,j+1))/4.)/sew(i)
         genres=gen(i,j)-vis(i,j)*(dvdx**2)
         gen(i,j)=genres+gencou
         if (xpluse(j).le.11.63) go to 171
         diterm=den(i,j)*(cmu**0.75)*sqrtk*alog(elog*xpluse(j))/
     1    (cappa*xp)
         go to 172
 171     continue
         diterm=den(i,j)*(cmu**0.75)*sqrtk*xpluse(j)/xp
 172     continue
         su(i,j)=sukd(i,j)+gen(i,j)*vol
         sp(i,j)=spkd(i,j)-diterm*vol
         ae(i,j)=0.0
      end do
      return
chapter  7  7  7  7  7  7  7  7  dissipation  7  7  7  7  7  7  7  7  7
c
      entry moded
c-----top boundry
      do i=2,nim1
         an(i,njm1)=0.0
      end do
c-----south wall left
      yp=y(2)
      j=2
      term=(cmu**0.75)/(cappa*yp)
      do i=2,ista
        if(te(i,j).lt.0.) te(i,j)=abs(te(i,j))
        su(i,j)=great*term*te(i,j)**1.5
        sp(i,j)=-great
      end do
c-----south wall middle
      yp=y(jst)-yv(jst)
      j=jst
      term=(cmu**0.75)/(cappa*yp)
      do i=istap1,istbm1
        if(te(i,j).lt.0.) te(i,j)=abs(te(i,j))
        su(i,j)=great*term*te(i,j)**1.5
        sp(i,j)=-great
      end do
c-----south wall right
      yp=y(2)
      j=2
      term=(cmu**0.75)/(cappa*yp)
      do i=istb,nim1
        if(te(i,j).lt.0.) te(i,j)=abs(te(i,j))
        su(i,j)=great*term*te(i,j)**1.5
        sp(i,j)=-great
      end do
c-----west wall
      xp=x(istb)-xu(istb)
      i=istb
      term=(cmu**.75)/(cappa*xp)
      do j=2,jstm1
         if(te(i,j).lt.0.) te(i,j)=abs(te(i,j))
         su(i,j)=great*term*te(i,j)**1.5
         sp(i,j)=-great
      end do
c-----east wall
      xp=xu(istap1)-x(ista)
      i=ista
      term=(cmu**.75)/(cappa*xp)
      do j=2,jstm1
         if(te(i,j).lt.0.) te(i,j)=abs(te(i,j))
         su(i,j)=great*term*te(i,j)**1.5
         sp(i,j)=-great
      end do
      return
      end
c....................... end of program ..............................




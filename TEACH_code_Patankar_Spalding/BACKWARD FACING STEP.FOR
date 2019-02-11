c     ---------------------------------------------------------------
c     |                   BACKWARD FACING STEP                      |
c     ---------------------------------------------------------------     
c     |                                                             |
c     |               -----------------------------                 |
c     |                 --->    ^ y    ---->                        |
c     |               ----------|      ---->                        |
c     |                         |------------------> x              |
c     |                                                             |
c     ---------------------------------------------------------------
      program back step
      include 'wall.inc'
c     .................

      call input
      call constant
      call grid
      call geometry
      call initial
c-----main  loop
      urfu=0.65
      urfv=0.65
      urfp=0.1
      urfd=0.65
      urfk=0.65
      urfvis=0.7
      do iter=1,50000
	call props
	call calcu     
	call calcv
	call calcp
	call calcte
	call calced
c-----normalised  residuals
	resoru=resoru/(flowin*uin)
	resorv=resorv/(flowin*uin)
	resorm=resorm/(flowin)
c-----termination  tests
	source=dmax1(resoru,resorv,resorm)
	sour=dmax1(resork,resord)
	write(*,*)iter,source,sour
        write(3,*)iter,source
	if (iter.ge.maxit.or.source.lt.sormax) go to 10
      enddo
10    call output
      stop
      end
c     ---------------------------------------------------------------
c                          subroutine input
c     ---------------------------------------------------------------
      subroutine input
      include 'wall.inc'
c     .................
      open(1,file='c:\tecplot\wall.dat')
      open(2,file='yplus.dat')
      open(3,file='source.dat')
      open(5,file='ent.dat')
c-----
      write(*,*)'select method for u(1=upwind,2=smart)'
      read(*,*)indu
      if(indu.eq.2) then
        write(*,*)'select boundary order for u(1=first,2=second)'
        read(*,*)indbu
      endif
      write(*,*)'select method for v(1=upwind,2=smart)'
      read(*,*)indv
      if(indv.eq.2) then 
        write(*,*)'select boundary order for v(1=first,2=second)'
        read(*,*)indbv
      endif
      write(*,*)'select method for k(1=upwind,2=smart)'
      read(*,*)indk
      if(indk.eq.2) then
        write(*,*)'select boundary order for k(1=first,2=second)'
        read(*,*)indbk
      endif
      write(*,*)'select method for eps(1=upwind,2=smart)'
      read(*,*)indd
      if(indd.eq.2) then
        write(*,*)'select boundary order for eps(1=first,2=second)'
        read(*,*)indbd
      endif
      write(*,*)'enter step location ist:'
      read(*,*)ist
      write(*,*)'enter step location jst(an even number):'
      read(*,*)jst
      write(*,*)'enter no. of iteration:'
      read(*,*) maxit
      return
      end
c     ---------------------------------------------------------------
c                          subroutine constant
c     ---------------------------------------------------------------
      subroutine constant
      include 'wall.inc'
c     .................

c-----step declearation
      istp1=ist+1
      istm1=ist-1
      jstp1=jst+1
      jstm1=jst-1
c-----entrance  values
      uin=12
      tein=0.0
      edin=0.0
c-----a constant great number
      great=1.e30                           
c-----max error       
      sormax=1e-5
c-----no of sweeps      
      nswpu=2 
      nswpv=2
      nswpp=10
      nswpk=2
      nswpd=2
c-----fluid properties
      densit=1.23
      viscos=1.8e-5
      prandt=0.71
c-----geometry info
      h=5.08e-2
      w1=7.62e-2
      w2=12.7e-2
      dr=30*h
      dl=3*h
c-----turbulence  constants
      cmu=0.09
      cd=1.00
      c1=1.44
      c2=1.92
      cappa=0.41
      elog=7.768
      pred=1.3
      prte=1.0
c-----pressure   calculation
      ipref=2
      jpref=jst
      return      
      end
c     ---------------------------------------------------------------
c                            subroutine grid
c     ---------------------------------------------------------------
      subroutine grid
      include 'wall.inc'
c     .................      

c-----step right side
      s=1.05
      sumx=1+(s**(ni-ist)-s)/(s-1)
      dxr=dr/sumx
      xu(ist)=0
      do i=istp1,ni
	 xu(i)=xu(i-1)+dxr
	 dxr=dxr*s
      end do
c-----step left side
      s=1.05
      sumx=1+(s**(ist-2)-s)/(s-1)
      dxl=dl/sumx
      xu(ist)=0
      do i=istm1,1,-1
	 xu(i)=xu(i+1)-dxl
	 dxl=dxl*s 
      enddo
c-----down wall
      s=1.1
      sumy=1+(s**(jst/2-1)-s)/(s-1)
      dyd=h/2/sumy
      yv(1)=-dyd
      yv(2)=0
      yv(jst)=h
      do j=3,jst/2+1
	 yv(j)=yv(j-1)+dyd
	 yv(jst-j+2)=yv(jst-j+3)-dyd
	 dyd=dyd*s
      enddo
c-----top wall
      s=1.1
      sumy=1+(s**((nj-jst)/2)-s)/(s-1)
      dyt=w1/2/sumy
      yv(nj)=h+w1
      do j=jstp1,(nj+jst)/2
	 yv(j)=yv(j-1)+dyt
	 yv(nj-j+jst)=yv(nj-j+jst+1)-dyt
	 dyt=dyt*s
      enddo
c      write(1,*)'variables=x,y'
c      write(1,*)'zone i=60,j=40'
c      do j=1,nj
c       do i=1,ni
c         write(1,*)xu(i),yv(j)
c       enddo
c      enddo
c      stop
      return
      end

c     ---------------------------------------------------------------
c                          subroutine geometry
c     ---------------------------------------------------------------
      subroutine geometry
      include 'wall.inc'
c     .................
      do i=1,nim1
	x(i)=(xu(i)+xu(i+1))/2
      enddo
      x(ni)=xu(ni)+(xu(ni)-x(nim1))
c-----
      do j=1,njm1
	y(j)=(yv(j)+yv(j+1))/2
      enddo
      y(nj)=yv(nj)+(yv(nj)-y(njm1))
c-----
      dxpw(1)=0.0
      dxep(ni)=0.0
      do i=1,nim1
	dxep(i)=x(i+1)-x(i)
	dxpw(i+1)=dxep(i)
      enddo
c-----
      dyps(1)=0.0
      dynp(nj)=0.0
      do j=1,njm1
	dynp(j)=y(j+1)-y(j)
	dyps(j+1)=dynp(j)
      enddo
c-----
      sew(1)=0.0
      sew(ni)=0.0
      do i=2,nim1
	sew(i)=xu(i+1)-xu(i)
      enddo
c-----
      sns(1)=0.0
      sns(nj)=0.0
      do j=2,njm1
	sns(j)=yv(j+1)-yv(j)
      enddo
c-----
      dxpwu(1)=0.0
      dxepu(ni)=0.0
      do i=1,nim1
	dxepu(i)=xu(i+1)-xu(i)
	dxpwu(i+1)=dxepu(i)
      enddo
c-----
      sewu(1)=0.0
      sewu(ni)=0.0
      do i=2,nim1
	sewu(i)=x(i)-x(i-1)
      enddo
c-----
      dypsv(1)=0.0
      dynpv(nj)=0.0
      do j=1,njm1
	dynpv(j)=yv(j+1)-yv(j)
	dypsv(j+1)=dynpv(j)
      enddo
c-----
      snsv(1)=0.0
      snsv(nj)=0.0
      do j=2,njm1
	snsv(j)=y(j)-y(j-1)
      enddo
c-----set variable to zero  
      do  i=1,ni
	do  j=1,nj
	   u(i,j)=0.0
	   v(i,j)=0.0
	   p(i,j)=0.0
	   pp(i,j)=0.0
	   te(i,j)=1.0e-10
	   ed(i,j)=0.0001
           vis(i,j)=viscos
	   du(i,j)=0.0
	   dv(i,j)=0.0
	   su(i,j)=0.0
	   sp(i,j)=0.0
	enddo
      enddo
      return
      end
c     ---------------------------------------------------------------
c                            subroutine initial
c     ---------------------------------------------------------------
      subroutine initial
      include 'wall.inc'
      dimension ye(69),f1(69),f2(69)
c     .................
      do i=1,69
	read(5,*)ye(i),f1(i),f2(i)
	ye(i)=ye(i)/100+h
      enddo
      do j=jst,njm1
	 do i=1,69
	   if(ye(i).gt.y(j)) then
	     ii=i
	     goto 10
	   endif
	 enddo
10       i=ii
	 fup=(y(j)-ye(i-1))/(ye(i)-ye(i-1))
	 fdown=1-fup
	 u(2,j)=f1(i)*fup+f1(i-1)*fdown
	 te(1,j)=f2(i)*fup+f2(i-1)*fdown
	 ed(1,j)=te(1,j)**1.5/(0.005*w1)
      enddo

c-----flow rate 
      flowin=0.0
      do j=jst,njm1
         arden=densit*sns(j)
	 flowin=flowin+arden*u(2,j)
      enddo
      do j=jst,njm1
	u(1,j)=u(2,j)
      enddo
      do i=3,ist
	 do j=jst,njm1
	    u(i,j)=uin  
	 enddo
      enddo
      do i=istp1,ni
	 do j=2,njm1
	    u(i,j)=uin*w1/w2   
	 enddo
      enddo
      do  i=2,nim1
	 yplusn(i)=11.0
	 ypluss(i)=11.0
      enddo
      do j=2,jstm1
	 xplusw(j)=11.0
      enddo   
      return      
      end
c     ---------------------------------------------------------------
c                           subroutine props
c     ---------------------------------------------------------------
      subroutine props
      include 'wall.inc'
c     .................

      do i=1,ni
	do j=2,njm1
	  if(i.ge.ist.or.j.ge.jst) then
	    visold=vis(i,j)
	    if (ed(i,j).eq.0.) then
	      vis(i,j)=viscos
	    else  
              vis(i,j)=densit*te(i,j)**2*cmu/ed(i,j)+viscos
	    endif
	    vis(i,j)=urfvis*vis(i,j)+(1-urfvis)*visold
	  endif
	enddo
      enddo
      return
      end
c     ---------------------------------------------------------------
c                           subroutine calcu
c     ---------------------------------------------------------------
      subroutine calcu
      include 'wall.inc'
c     .................

c-----assembling  coefficients
      do i=3,nim1
	do j=2,njm1
	  if(i.gt.ist.or.j.ge.jst) then
c-----computing areas and volume
	    arean=sewu(i)
	    areas=sewu(i)
	    areaew=sns(j)
	    vol=sewu(i)*sns(j)
c-----calculating convection coefficients
            gn=densit*v(i,j+1)
            gnw=densit*v(i-1,j+1)
            gs=densit*v(i,j)
            gsw=densit*v(i-1,j)
            ge=densit*u(i+1,j)
            gp=densit*u(i,j)
            gw=densit*u(i-1,j)
	    cn=(gn*dxpwu(i)+gnw*dxepu(i))/(dxpwu(i)+dxepu(i))*arean
	    cs=(gs*dxpwu(i)+gsw*dxepu(i))/(dxpwu(i)+dxepu(i))*areas
	    ce=0.5*(ge+gp)*areaew
	    cw=0.5*(gp+gw)*areaew
c-----calculating diffusion coefficients
            visn=((vis(i,j)*sew(i-1)+vis(i-1,j)*sew(i))
     &             /(sew(i)+sew(i-1))*sns(j+1)
     &             +(vis(i,j+1)*sew(i-1)+vis(i-1,j+1)*sew(i))
     &             /(sew(i)+sew(i-1))*sns(j))
     &             /(sns(j)+sns(j+1))

            viss=((vis(i,j)*sew(i-1)+vis(i-1,j)*sew(i))
     &             /(sew(i)+sew(i-1))*sns(j-1)
     &             +(vis(i,j-1)*sew(i-1)+vis(i-1,j-1)*sew(i))
     &             /(sew(i)+sew(i-1))*sns(j))
     &             /(sns(j)+sns(j-1))

            dn=visn*arean/dynp(j)
	    ds=viss*areas/dyps(j)
	    de=vis(i,j)*areaew/dxepu(i)
	    dw=vis(i-1,j)*areaew/dxpwu(i)
c-----calculating coefficients of source terms
	    smp=cn-cs+ce-cw
	    cp=dmax1(0.0,smp)
	    cpo=cp
c-----assembling main coefficients
c           an(i,j)=dn*dmax1(0,1-0.5*abs(cn/dn))+dmax1(-cn,0)
c           as(i,j)=ds*dmax1(0,1-0.5*abs(cs/ds))+dmax1(cs,0)
c           ae(i,j)=de*dmax1(0,1-0.5*abs(ce/de))+dmax1(-ce,0)
c           aw(i,j)=dw*dmax1(0,1-0.5*abs(cw/dw))+dmax1(cw,0)

	    an(i,j)=dmax1(0,-cn)+dn
	    as(i,j)=dmax1(0,cs)+ds
	    ae(i,j)=dmax1(0,-ce)+de
	    aw(i,j)=dmax1(0,cw)+dw

c-----second order operation
c---------------------------
	    sum=0
	    if(indu.eq.1.or.iter.lt.200) goto 80
c-----east face
	    if(ce.ge.0) then
	      phic=u(i,j)
	      phid=u(i+1,j)
	      phiu=u(i-1,j)
	      rd=dxepu(i)
	      ru=dxepu(i-1)
	      rf=x(i)-xu(i)
	    else 
	      phic=u(i+1,j)
	      phid=u(i,j)
	      phiu=u(i+2,j)
	      rd=dxepu(i)
	      ru=dxepu(i+1)
	      rf=xu(i+1)-x(i)
	    endif
            call calc
            sum=sum-ce*(phiho-phic)
c-----west face
      if(cw.ge.0) then
          phic=u(i-1,j)
          phid=u(i,j)
          phiu=u(i-2,j)
          rd=dxpwu(i)
          ru=dxpwu(i-1)
          rf=x(i-1)-xu(i-1)
      else
          phic=u(i,j)
          phid=u(i-1,j)
          phiu=u(i+1,j)
          rd=dxpwu(i)
          ru=dxpwu(i+1)
          rf=xu(i)-x(i-1)
      endif
      if(cw.ge.0.and.(i.eq.istp1.and.j.lt.jst)) then
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cw*(phiho-phic)
c-----north face
      if(j.eq.njm1) goto 31
      if(cn.ge.0) then
          phic=u(i,j)
          phid=u(i,j+1)
          phiu=u(i,j-1)
          rd=dynp(j)
          ru=dynp(j-1)
          rf=yv(j+1)-y(j) 
      else
          phic=u(i,j+1)
          phid=u(i,j)
          phiu=u(i,j+2)
          rd=dynp(j)
          ru=dynp(j+1)
          rf=y(j+1)-yv(j+1)
      endif
      if(j.eq.2.and.cn.ge.0) then
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(j.eq.nj-2.and.cn.lt.0) then
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cn.ge.0.and.(j.eq.jst.and.i.le.ist)) then
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum-cn*(phiho-phic)
31    continue
c-----south face
      if(j.eq.2.or.(j.eq.jst.and.i.le.ist)) goto 41
      if(cs.ge.0) then
          phic=u(i,j-1)
          phid=u(i,j)
          phiu=u(i,j-2)
          rd=dyps(j)
          ru=dyps(j-1)
          rf=yv(j)-y(j-1)
      else
          phic=u(i,j)
          phid=u(i,j-1)
          phiu=u(i,j+1)
          rd=dyps(j)
          ru=dyps(j+1)
          rf=y(j)-yv(j) 
      endif
      if(j.eq.3.and.cs.ge.0) then  
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(j.eq.njm1.and.cs.lt.0) then  
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cs.ge.0.and.(j.eq.jstp1.and.i.le.ist)) then  
           if(indbu.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cs*(phiho-phic)
41    continue
c--------------------------
80          du(i,j)=areaew
            sse=vis(i,j)*(u(i+1,j)-u(i,j))/dxepu(i)*sns(j)
            ssw=vis(i-1,j)*(u(i,j)-u(i-1,j))/dxpwu(i)*sns(j)
            ssn=visn*(v(i,j+1)-v(i-1,j+1))
            sss=viss*(v(i,j)-v(i-1,j))
            ssu=sse-ssw+ssn-sss
            su(i,j)=ssu+cpo*u(i,j)+du(i,j)*(p(i-1,j)-p(i,j))+sum
	    sp(i,j)=-cp
	  endif
	enddo   
      enddo



c
c-----problem modifications  
c
      call modu
c
c-----final coeff. assembly and residual source calculations  
c
      resoru=0.0
      do i=3,nim1
	do j=2,njm1
	  if(i.gt.ist.or.j.ge.jst) then
	    ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j) 
	    du(i,j)=du(i,j)/ap(i,j)
	    resor=an(i,j)*u(i,j+1)+as(i,j)*u(i,j-1)+ae(i,j)*u(i+1,j)
     &            +aw(i,j)*u(i-1,j)-ap(i,j)*u(i,j)+su(i,j)
	    vol=sew(i)*sns(j)
	    sorvol=great*vol
	    if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
	    resoru=resoru+abs(resor)
c-----under relaxation
	    ap(i,j)=ap(i,j)/urfu
	    su(i,j)=su(i,j)+(1.-urfu)*ap(i,j)*u(i,j)
	    du(i,j)=du(i,j)*urfu
	  endif
	enddo
      enddo
c-----solution of difference equations  
      do n=1,nswpu
	call lisolv(3,ist,jst,njm1,u)
	call lisolv(istp1,nim1,2,njm1,u)
      enddo
      return     
      end
c     ---------------------------------------------------------------
c                            subroutine calcv
c     ---------------------------------------------------------------
      subroutine calcv
      include 'wall.inc'
c     .................

c-----assembly of coefficients  
      do i=2,nim1
	do j=3,njm1
	  if(i.ge.ist.or.j.gt.jst) then
c-----compute areas and volume
	    arean=sew(i)
	    areas=sew(i)
	    areaew=snsv(j)
	    vol=sew(i)*snsv(j)
c-----calculate convection coefficients
            gn=densit*v(i,j+1)
            gp=densit*v(i,j)
            gs=densit*v(i,j-1)
            ge=densit*u(i+1,j)
            gse=densit*u(i+1,j-1)
            gw=densit*u(i,j)
            gsw=densit*u(i,j-1)
	    cn=0.5*(gn+gp)*arean
	    cs=0.5*(gs+gp)*areas
	    ce=(ge*dypsv(j)+gse*dynpv(j))/(dypsv(j)+dynpv(j))*areaew
	    cw=(gw*dypsv(j)+gsw*dynpv(j))/(dypsv(j)+dynpv(j))*areaew
c-----calculate diffusion coefficients

            vise=((sns(j)*vis(i,j-1)+sns(j-1)*vis(i,j))
     &             /(sns(j)+sns(j-1))*sew(i+1)
     &             +(sns(j)*vis(i+1,j-1)+sns(j-1)*vis(i+1,j))
     &             /(sns(j)+sns(j-1))*sew(i))
     &             /(sew(i)+sew(i+1))
            visw=((sns(j)*vis(i,j-1)+sns(j-1)*vis(i,j))
     &             /(sns(j)+sns(j-1))*sew(i-1)
     &             +(sns(j)*vis(i-1,j-1)+sns(j-1)*vis(i-1,j))
     &             /(sns(j)+sns(j-1))*sew(i))
     &             /(sew(i)+sew(i-1))

	    dn=vis(i,j)*arean/dynpv(j)
	    ds=vis(i,j-1)*areas/dypsv(j)
	    de=vise*areaew/dxep(i)
	    dw=visw*areaew/dxpw(i)
c-----calculate coefficients of source terms
	    smp=cn-cs+ce-cw
	    cp=dmax1(0.0,smp)
	    cpo=cp
c-----assemble main coefficients
c           an(i,j)=dn*dmax1(0,1-0.5*abs(cn/dn))+dmax1(-cn,0)
c           as(i,j)=ds*dmax1(0,1-0.5*abs(cs/ds))+dmax1(cs,0)
c           ae(i,j)=de*dmax1(0,1-0.5*abs(ce/de))+dmax1(-ce,0)
c           aw(i,j)=dw*dmax1(0,1-0.5*abs(cw/dw))+dmax1(cw,0)

	    an(i,j)=dmax1(0,-cn)+dn
	    as(i,j)=dmax1(0,cs)+ds
	    ae(i,j)=dmax1(0,-ce)+de
	    aw(i,j)=dmax1(0,cw)+dw
c-----second order operation
c---------------------------
	    sum=0
	    if(indv.eq.1.or.iter.lt.200) goto 80
c-----east face
      if(ce.ge.0) then
           phic=v(i,j)
           phid=v(i+1,j)
           phiu=v(i-1,j)
           rd=dxep(i)
           ru=dxep(i-1)
           rf=xu(i+1)-x(i)
      else 
           phic=v(i+1,j)
           phid=v(i,j)
           phiu=v(i+2,j)
           rd=dxep(i)
           ru=dxep(i+1)
           rf=x(i+1)-xu(i+1)
      endif 
      if(ce.ge.0.and.(i.eq.ist.and.j.le.jst)) then 
           if(indbv.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum-ce*(phiho-phic)
11    continue
c-----west face
      if(i.eq.ist.and.j.le.jst) goto 21
      if(cw.ge.0) then
           phic=v(i-1,j)
           phid=v(i,j)
           phiu=v(i-2,j)
           rd=dxpw(i)
           ru=dxpw(i-1)
           rf=xu(i)-x(i-1)
      else
           phic=v(i,j)
           phid=v(i-1,j)
           phiu=v(i+1,j)
           rd=dxpw(i)
           ru=dxpw(i+1)
           rf=x(i)-xu(i)
      endif
      if(i.eq.2) then
           if(indbv.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cw.ge.0.and.(i.eq.istp1.and.j.le.jst)) then
           if(indbv.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cw*(phiho-phic)
21    continue
c-----north face
      if(cn.ge.0) then
           phic=v(i,j)
           phid=v(i,j+1)
           phiu=v(i,j-1)
           rd=dynpv(j)
           ru=dynpv(j-1)
           rf=y(j)-yv(j)
      else
           phic=v(i,j+1)
           phid=v(i,j)
           phiu=v(i,j+2)
           rd=dynpv(j)
           ru=dynpv(j+1)
           rf=yv(j+1)-y(j)
      endif
      if(j.eq.njm1.and.cn.lt.0) then 
           if(indbv.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else      
           call calc
      endif
      sum=sum-cn*(phiho-phic)
31    continue
c-----south face
      if(cs.ge.0) then
           phic=v(i,j-1)
           phid=v(i,j)
           phiu=v(i,j-2)
           rd=dypsv(j)
           ru=dypsv(j-1)
           rf=y(j-1)-yv(j-1)
      else
           phic=v(i,j)
           phid=v(i,j-1)
           phiu=v(i,j+1)
           rd=dypsv(j)
           ru=dypsv(j+1)
           rf=yv(j)-y(j-1)
      endif
      if(j.eq.3.and.cs.ge.0) then
           if(indbv.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cs.ge.0.and.(j.eq.jstp1.and.i.lt.ist)) then
           if(indbv.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cs*(phiho-phic)
41    continue
c--------------------------c
80          dv(i,j)=0.5*(arean+areas)
            sse=vise*(u(i+1,j)-u(i+1,j-1))
            ssw=visw*(u(i,j)-u(i,j-1))
            ssn=vis(i,j)*(v(i,j+1)-v(i,j))/dynpv(j)*sew(i)
            sss=vis(i,j-1)*(v(i,j)-v(i,j-1))/dypsv(j)*sew(i)
            ssu=sse-ssw+ssn-sss
            su(i,j)=ssu+cpo*v(i,j)+dv(i,j)*(p(i,j-1)-p(i,j))+sum
	    sp(i,j)=-cp
	  endif
	enddo
      enddo
c
c-----problem modifications 
c
      call modv
c
c-----final coeff. assembly and residual source calculations  
c
      resorv=0.0
      do i=2,nim1
	do j=3,njm1
	  if(i.ge.ist.or.j.gt.jst) then
	    ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
	    dv(i,j)=dv(i,j)/ap(i,j)
	    resor=an(i,j)*v(i,j+1)+as(i,j)*v(i,j-1)+ae(i,j)*v(i+1,j)
     &            +aw(i,j)*v(i-1,j)-ap(i,j)*v(i,j)+su(i,j)
	    vol=sew(i)*sns(j)
	    sorvol=great*vol
	    if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
	    resorv=resorv+abs(resor)
c-----under relaxation
	    ap(i,j)=ap(i,j)/urfv
	    su(i,j)=su(i,j)+(1.-urfv)*ap(i,j)*v(i,j)
	    dv(i,j)=dv(i,j)*urfv
	  endif
	enddo
      enddo
c-----solution of difference equations  
      do n=1,nswpv
	call lisolv(2,istm1,jstp1,njm1,v)
	call lisolv(ist,nim1,3,njm1,v)
      enddo
      return
      end
c     ---------------------------------------------------------------
c                             subroutine calcp
c     ---------------------------------------------------------------
      subroutine calcp
      include 'wall.inc'
c     .................

c-----assembly of coefficients   
      resorm=0.0
      do i=2,nim1
	do j=2,njm1
	  if(i.ge.ist.or.j.ge.jst) then
c-----compute areas and volume
	    arean=sew(i)
	    areas=sew(i)
	    areaew=sns(j)
	    vol=sns(j)*sew(i)
c-----calculate coefficients
            denn=densit
            dens=densit
            dene=densit
            denw=densit
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
	  endif
	enddo
      enddo
c
c-----problem modifications  
c
      call modp
c
c-----final coeff. assembly and residual source calculations  
c
      do i=2,nim1
	do j=2,njm1
	  if(i.ge.ist.or.j.ge.jst) then
	    ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
	  endif
	enddo
      enddo
c-----solution of difference equations  
      do n=1,nswpp
	call lisolv(2,istm1,jst,njm1,pp)
	call lisolv(ist,nim1,2,njm1,pp)
      enddo
c-----velocity correction
      do i=2,nim1
	do j=2,njm1
	  if(i.ge.ist.or.j.ge.jst) then
	    u(i,j)=u(i,j)+du(i,j)*(pp(i-1,j)-pp(i,j))
	    v(i,j)=v(i,j)+dv(i,j)*(pp(i,j-1)-pp(i,j))
	  endif
	enddo
      enddo
c-----pressure correction with under relaxation                                                                                                                                                                                                                
      ppref=pp(ipref,jpref)
      do i=2,nim1
	do j=2,njm1
	  if(i.ge.ist.or.j.ge.jst) then
	    p(i,j)=p(i,j)+urfp*(pp(i,j)-ppref)
	    pp(i,j)=0.0
	  endif
	enddo
      enddo
      return
      end
c     ---------------------------------------------------------------
c                            subroutine calcte
c     ---------------------------------------------------------------
      subroutine calcte
      include 'wall.inc'
c     .................

c-----assembly of coefficients
      prte=1.0
      call modte
      do i=2,nim1
       do j=3,nj-2
	if(i.gt.ist.or.j.gt.jst.or.(i.eq.ist.and.j.eq.jst)) then
c-----compute areas and volume
	  arean=sew(i)
	  areas=sew(i)
	  areaew=sns(j)
	  vol=sns(j)*sew(i)
c-----calculate convection coefficients
          gn=densit*v(i,j+1)
          gs=densit*v(i,j)
          ge=densit*u(i+1,j)
          gw=densit*u(i,j)
	  cn=gn*arean
	  cs=gs*areas
	  ce=ge*areaew
	  cw=gw*areaew
c-----calculate diffusion coefficients
	  gamn=(sns(j)*vis(i,j+1)+sns(j+1)*vis(i,j))/
     &         (sns(j)+sns(j+1))/prte
	  gams=(sns(j)*vis(i,j-1)+sns(j-1)*vis(i,j))/
     &         (sns(j)+sns(j-1))/prte
	  game=(sew(i)*vis(i+1,j)+sew(i+1)*vis(i,j))/
     &         (sew(i)+sew(i+1))/prte
	  gamw=(sew(i)*vis(i-1,j)+sew(i-1)*vis(i,j))/
     &         (sew(i)+sew(i-1))/prte
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
	  dudy=((sns(j)*(u(i,j+1)+u(i+1,j+1))/2
     &          +sns(j+1)*(u(i,j)+u(i+1,j))/2)/(sns(j)+sns(j+1))
     &          -(sns(j)*(u(i,j-1)+u(i+1,j-1))/2
     &          +sns(j-1)*(u(i,j)+u(i+1,j))/2)/(sns(j)+sns(j-1)))
     &          /sns(j)
	  dvdx=((sew(i)*(v(i+1,j)+v(i+1,j+1))/2
     &           +sew(i+1)*(v(i,j)+v(i,j+1))/2)/(sew(i)+sew(i+1))
     &           -(sew(i)*(v(i-1,j)+v(i-1,j+1))/2
     &           +sew(i-1)*(v(i,j)+v(i,j+1))/2)/(sew(i)+sew(i-1)))
     &           /sew(i)
      gen(i,j)=(2.*(dudx**2+dvdy**2)+(dudy+dvdx)**2)*(vis(i,j)-viscos)
c-----assemble main coefficients
	    an(i,j)=dmax1(0,-cn)+dn
	    as(i,j)=dmax1(0,cs)+ds
	    ae(i,j)=dmax1(0,-ce)+de
	    aw(i,j)=dmax1(0,cw)+dw
c-----second order operation
c---------------------------
	    sum=0
            if(indk.eq.1.or.iter.lt.400) goto 80
c-----east face
      if(ce.ge.0) then
           phic=te(i,j)
           phid=te(i+1,j)
           phiu=te(i-1,j)
           rd=dxep(i)
           ru=dxep(i-1)
           rf=xu(i+1)-x(i)
           call calc
           sum=sum-ce*(phiho-phic)
      else 
           phic=te(i+1,j)
           phid=te(i,j)
           phiu=te(i+2,j)
           rd=dxep(i)
           ru=dxep(i+1)
           rf=x(i+1)-xu(i+1)
           call calc
           sum=sum-ce*(phiho-phic)
      endif
c-----west face
      if(cw.ge.0) then
           phic=te(i-1,j)
           phid=te(i,j)
           phiu=te(i-2,j)
           rd=dxpw(i)
           ru=dxpw(i-1)
           rf=xu(i)-x(i-1)
      else
           phic=te(i,j)
           phid=te(i-1,j)
           phiu=te(i+1,j)
           rd=dxpw(i)
           ru=dxpw(i+1)
           rf=x(i)-xu(i)
      endif
      if(i.eq.2) then   
           if(indbk.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cw.ge.0.and.(i.eq.istp1.and.j.le.jstm1)) then
           if(indbk.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cw*(phiho-phic)
c-----north face
      if(cn.ge.0) then
           phic=te(i,j)
           phid=te(i,j+1)
           phiu=te(i,j-1)
           rd=dynp(j)
           ru=dynp(j-1)
           rf=yv(j+1)-y(j) 
      else
           phic=te(i,j+1)
           phid=te(i,j)
           phiu=te(i,j+2)
           rd=dynp(j)
           ru=dynp(j+1)
           rf=y(j+1)-yv(j+1)
      endif
      if(j.eq.nj-2.and.cn.lt.0) then     
           if(indbk.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else 
           call calc
      endif
      sum=sum-cn*(phiho-phic)
c-----south face
      if(cs.ge.0) then
           phic=te(i,j-1)
           phid=te(i,j)
           phiu=te(i,j-2)
           rd=dyps(j)
           ru=dyps(j-1)
           rf=yv(j)-y(j-1)
      else
           phic=te(i,j)
           phid=te(i,j-1)
           phiu=te(i,j+1)
           rd=dyps(j)
           ru=dyps(j+1)
           rf=y(j)-yv(j) 
      endif
      if(j.eq.3.and.cs.ge.0) then     
           if(indbk.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cs.ge.0.and.(j.eq.jstp1.and.i.le.istm1)) then
           if(indbk.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cs*(phiho-phic)

c--------------------------

80        su(i,j)=cpo*te(i,j)+sum
	  su(i,j)=su(i,j)+gen(i,j)*vol
	  sp(i,j)=-cp
          if(iter.lt.100) then
            sp(i,j)=sp(i,j)-cd*cmu*densit**2*te(i,j)*vol/vis(i,j)
          else
        sp(i,j)=sp(i,j)-cd*cmu*densit**2*te(i,j)*vol/(vis(i,j)-viscos)
          endif
	endif
       enddo
      enddo
c
c-----problem modifications 
c
c      call modte
c
c-----final coeff. assembly and residual source calculations  
c
      resork=0.0
      do i=2,nim1
       do j=3,nj-2
       if(i.gt.ist.or.j.gt.jst.or.(i.eq.ist.and.j.eq.jst)) then
	 ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
	 resor=an(i,j)*te(i,j+1)+as(i,j)*te(i,j-1)+ae(i,j)*te(i+1,j)
     1         +aw(i,j)*te(i-1,j)-ap(i,j)*te(i,j)+su(i,j)
	 vol=sew(i)*sns(j)
	 sorvol=great*vol
	 if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
	 resork=resork+abs(resor)
c-----under relaxation
	 ap(i,j)=ap(i,j)/urfk
	 su(i,j)=su(i,j)+(1.-urfk)*ap(i,j)*te(i,j)
	endif
       enddo
      enddo
c-----solution of difference equations  
      do n=1,nswpk
	call lisolv(2,istm1,jstp1,nj-2,te)
	call lisolv(ist,ist,jst,nj-2,te)
	call lisolv(istp1,nim1,3,nj-2,te)
      enddo
      return
      end
c     ---------------------------------------------------------------
c                           subroutine calced
c     ---------------------------------------------------------------
      subroutine calced
      include 'wall.inc'
c     .................

c-----assembly of coefficients   
      call moded
      do i=2,nim1
       do j=3,nj-2
	if(i.gt.ist.or.j.gt.jst.or.(i.eq.ist.and.j.eq.jst)) then
c-----compute areas and volume
	 arean=sew(i)
	 areas=sew(i)
	 areaew=sns(j)
	 vol=sns(j)*sew(i)
c-----calculate convection coefficients
         gn=densit*v(i,j+1)
         gs=densit*v(i,j)
         ge=densit*u(i+1,j)
         gw=densit*u(i,j)
	 cn=gn*arean
	 cs=gs*areas
	 ce=ge*areaew
	 cw=gw*areaew
c-----calculate diffusion coefficients
	  gamn=(sns(j)*vis(i,j+1)+sns(j+1)*vis(i,j))/
     &         (sns(j)+sns(j+1))/pred-viscos/pred+viscos
	  gams=(sns(j)*vis(i,j-1)+sns(j-1)*vis(i,j))/
     &         (sns(j)+sns(j-1))/pred-viscos/pred+viscos
	  game=(sew(i)*vis(i+1,j)+sew(i+1)*vis(i,j))/
     &         (sew(i)+sew(i+1))/pred-viscos/pred+viscos
	  gamw=(sew(i)*vis(i-1,j)+sew(i-1)*vis(i,j))/
     &         (sew(i)+sew(i-1))/pred-viscos/pred+viscos
	  dn=gamn*arean/dynp(j)
	  ds=gams*areas/dyps(j)
	  de=game*areaew/dxep(i)
	  dw=gamw*areaew/dxpw(i)
c-----source terms
	 smp=cn-cs+ce-cw
	 cp=dmax1(0.0,smp)
	 cpo=cp
c-----assemble main coefficients
	    an(i,j)=dmax1(0,-cn)+dn
	    as(i,j)=dmax1(0,cs)+ds
	    ae(i,j)=dmax1(0,-ce)+de
	    aw(i,j)=dmax1(0,cw)+dw
c-----second order operation
c---------------------------
	    sum=0
            if(indd.eq.1.or.iter.lt.400) goto 80
c-----east face
      if(ce.ge.0) then
           phic=ed(i,j)
           phid=ed(i+1,j)
           phiu=ed(i-1,j)
           rd=dxep(i)
           ru=dxep(i-1)
           rf=xu(i+1)-x(i)
           call calc
           sum=sum-ce*(phiho-phic)
      else 
           phic=ed(i+1,j)
           phid=ed(i,j)
           phiu=ed(i+2,j)
           rd=dxep(i)
           ru=dxep(i+1)
           rf=x(i+1)-xu(i+1)
           call calc
           sum=sum-ce*(phiho-phic)
      endif
c-----west face
      if(cw.ge.0) then
           phic=ed(i-1,j)
           phid=ed(i,j)
           phiu=ed(i-2,j)
           rd=dxpw(i)
           ru=dxpw(i-1)
           rf=xu(i)-x(i-1)
      else
           phic=ed(i,j)
           phid=ed(i-1,j)
           phiu=ed(i+1,j)
           rd=dxpw(i)
           ru=dxpw(i+1)
           rf=x(i)-xu(i)
      endif
      if(i.eq.2) then   
           if(indbd.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cw.ge.0.and.(i.eq.istp1.and.j.le.jstm1)) then
           if(indbd.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cw*(phiho-phic)
c-----north face
      if(cn.ge.0) then
           phic=ed(i,j)
           phid=ed(i,j+1)
           phiu=ed(i,j-1)
           rd=dynp(j)
           ru=dynp(j-1)
           rf=yv(j+1)-y(j) 
      else
           phic=ed(i,j+1)
           phid=ed(i,j)
           phiu=ed(i,j+2)
           rd=dynp(j)
           ru=dynp(j+1)
           rf=y(j+1)-yv(j+1)
      endif
      if(j.eq.nj-2.and.cn.lt.0) then     
           if(indbd.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else 
           call calc
      endif
      sum=sum-cn*(phiho-phic)
c-----south face
      if(cs.ge.0) then
           phic=ed(i,j-1)
           phid=ed(i,j)
           phiu=ed(i,j-2)
           rd=dyps(j)
           ru=dyps(j-1)
           rf=yv(j)-y(j-1)
      else
           phic=ed(i,j)
           phid=ed(i,j-1)
           phiu=ed(i,j+1)
           rd=dyps(j)
           ru=dyps(j+1)
           rf=y(j)-yv(j) 
      endif
      if(j.eq.3.and.cs.ge.0) then     
           if(indbd.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      elseif(cs.ge.0.and.(j.eq.jstp1.and.i.le.istm1)) then
           if(indbd.eq.1) then
             phiho=phic
           else
             phiho=rf/rd*phid+(1-rf/rd)*phic
           endif
      else
           call calc
      endif
      sum=sum+cs*(phiho-phic)

c--------------------------

80       su(i,j)=cpo*ed(i,j)+sum
         if(iter.lt.100) then
           su(i,j)=su(i,j)+c1*cmu*gen(i,j)*vol*densit*te(i,j)
     &             /(vis(i,j))
         else
           su(i,j)=su(i,j)+c1*cmu*gen(i,j)*vol*densit*te(i,j)
     &             /(vis(i,j)-viscos)
         endif
	 sp(i,j)=-cp
         sp(i,j)=sp(i,j)-c2*densit*ed(i,j)*vol/te(i,j)
	endif
       enddo
      enddo
c
c-----problem modifications  
c
c      call moded
c
c-----final coeff. assembly and residual source calculations 
c
      resord=0.0
      do i=2,nim1
       do j=3,nj-2
        if(i.gt.ist.or.j.gt.jst.or.(i.eq.ist.and.j.eq.jst)) then
	  ap(i,j)=an(i,j)+as(i,j)+ae(i,j)+aw(i,j)-sp(i,j)
	  resor=an(i,j)*ed(i,j+1)+as(i,j)*ed(i,j-1)+ae(i,j)*ed(i+1,j)
     1          +aw(i,j)*ed(i-1,j)-ap(i,j)*ed(i,j)+su(i,j)
	  vol=sew(i)*sns(j)
	  sorvol=great*vol
	  if (-sp(i,j).gt.0.5*sorvol) resor=resor/sorvol
	  resord=resord+abs(resor)
c-----under relaxation
	  ap(i,j)=ap(i,j)/urfd
	  su(i,j)=su(i,j)+(1.-urfd)*ap(i,j)*ed(i,j)
	endif
       enddo
      enddo
c-----solution of difference equations 
      do n=1,nswpd
	call lisolv(2,istm1,jstp1,nj-2,ed)
	call lisolv(ist,ist,jst,nj-2,ed)
	call lisolv(istp1,nim1,3,nj-2,ed)
      enddo
      return
      end
c     ---------------------------------------------------------------
c                            subroutine lisolv
c     ---------------------------------------------------------------
      subroutine lisolv(istart,iend,jstart,jend,phi)
      include 'wall.inc'
      dimension phi(ni,nj),a(nj),b(nj),c(nj),d(nj)
c     ....................      

c-----commence w-e sweep
      do i=istart,iend
c-----commence s-n traverse
	do j=jstart,jend
c-----assemble toma coefficients
          if(j.eq.jend) then
            a(j)=0
          else
            a(j)=an(i,j)
          endif
          if(j.eq.jstart) then
            b(j)=0
          else
            b(j)=as(i,j)
          endif
          if(j.eq.jstart) then
              c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
     &             +as(i,j)*phi(i,j-1)
          elseif(j.eq.jend) then
              c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
     &             +an(i,j)*phi(i,j+1)
          else
              c(j)=ae(i,j)*phi(i+1,j)+aw(i,j)*phi(i-1,j)+su(i,j)
          endif
	  d(j)=ap(i,j)
c-----calculate coefficients of recorrence formula
          if(j.eq.jstart) then
            term=1./d(j)
          else
            term=1./(d(j)-b(j)*a(j-1))
          endif
	  a(j)=a(j)*term
          if(j.eq.jstart) then
            c(j)=c(j)*term
          else
            c(j)=(c(j)+b(j)*c(j-1))*term
          endif
	enddo
c-----obtain new phi's
	do j=jend,jstart,-1
          if(j.eq.jend) then
            phi(i,j)=c(j)
          else
            phi(i,j)=a(j)*phi(i,j+1)+c(j)
          endif
	enddo
      enddo
      return
      end
c     ---------------------------------------------------------------
c                           subroutine promod
c     ---------------------------------------------------------------
      subroutine promod
      include 'wall.inc'
c     .................

c     ---------------------------------
c          property modification  
c     ---------------------------------
c-----no modifications require
      entry modpro
      return

c     ---------------------------------
c            u-mom modification  
c     ---------------------------------
      entry modu
c-----top wall
      cdterm=cmu**0.25
      yp=yv(nj)-y(njm1)
      j=njm1
      do i=3,nim1
	 sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
         denu=densit
	 yplusa=0.5*(yplusn(i)+yplusn(i-1))
	 if (yplusa.le.10.8) then 
	   tmult=viscos/yp
	 else
	   tmult=denu*cdterm*sqrtk*cappa/alog(elog*yplusa)
	 endif
	 sp(i,j)=sp(i,j)-tmult*sewu(i)
	 an(i,j)=0.0
      end do
c-----south wall left
      cdterm=cmu**0.25
      yp=y(jst)-yv(jst)
      j=jst
      do i=3,ist        
	 sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
         denu=densit
	 if(i.eq.ist) then
	   yplusa=ypluss(i-1)
	 else
	   yplusa=0.5*(ypluss(i)+ypluss(i-1))
	 endif
	 if (yplusa.le.10.8) then 
	   tmult=viscos/yp
	 else
	   tmult=denu*cdterm*sqrtk*cappa/alog(elog*yplusa)
	 endif
	 sp(i,j)=sp(i,j)-tmult*sewu(i)
	 if(i.eq.ist) then
	   su(i,j)=su(i,j)+densit*(x(i)-xu(i))*v(i,j)/2*u(i,j)
	 endif
	 as(i,j)=0.0
      end do
c-----south wall right
      cdterm=cmu**0.25
      yp=y(2)
      j=2
      do i=istp1,nim1
	 sqrtk=sqrt(0.5*(te(i,j)+te(i-1,j)))
         denu=densit
	 yplusa=0.5*(ypluss(i)+ypluss(i-1))
	 if (yplusa.le.10.8) then 
	   tmult=viscos/yp
	 else
	   tmult=denu*cdterm*sqrtk*cappa/alog(elog*yplusa)
	 endif
	 sp(i,j)=sp(i,j)-tmult*sewu(i)
	 as(i,j)=0.0
      enddo
c-----outlet
      ardent=0.0
      flow=0.0
      do j=2,njm1
        arden=densit*sns(j)
	ardent=ardent+arden
	flow=flow+arden*u(nim1,j)
      enddo
      uinc=(flowin-flow)/ardent
      do j=2,njm1
	u(ni,j)=u(nim1,j)+uinc
      enddo
      return
c     ---------------------------------
c            v-mom modification  
c     ---------------------------------
      entry modv
c-----west wall
      cdterm=cmu**0.25
      xp=x(ist)
      i=ist
      do j=3,jst         
	 sqrtk=sqrt(0.5*(te(i,j)+te(i,j-1)))
         denv=densit
	 if(j.eq.jst) then
	   xplusa=xplusw(j-1)
	 else
	   xplusa=0.5*(xplusw(j)+xplusw(j-1))
	 endif
	 if (xplusa.le.10.8) then
	   tmult=viscos/xp
	 else
	   tmult=denv*cdterm*sqrtk*cappa/alog(elog*xplusa)
	 endif
	 sp(i,j)=sp(i,j)-tmult*snsv(j)
	 if(j.eq.jst) then
           su(i,j)=su(i,j)+densit*(y(j)-yv(j))*u(i,j)/2*v(i,j)
	 endif
	 aw(i,j)=0.0
      enddo
c-----outlet
      do j=3,njm1
        v(ni,j)=v(nim1,j)+dxep(nim1)/dxpw(nim1)*(v(nim1,j)-v(ni-2,j))
      enddo
      return
c     ---------------------------------
c           pressure modification  
c     ---------------------------------
c-----no modification requires      
      entry modp
      return
c     ---------------------------------
c             energy modification  
c     ---------------------------------
c-----no modification requires      
      entry modt
      return
c     ---------------------------------
c            k-trans modification  
c     ---------------------------------
      entry modte
c-----top boundry
      cdterm=sqrt(cmu)
      yp=yv(nj)-y(njm1)
      j=njm1
      do i=2,nim1
	 up=(u(i,j)+u(i+1,j))/2
	 ueff=abs(up)
	 utaw=sqrt(ueff*viscos/yp/densit)
	 yplusn(i)=utaw*yp*densit/viscos
	 if(yplusn(i).le.10.8) then
	   te(i,j)=utaw**2/cdterm*(yplusn(i)/10.8)**2
	 else
	   root=cappa*ueff
	   if(iter.gt.100) then
	     call newton(yp,ueff,root)
	   endif
	   utaw=root
	   yplusn(i)=utaw*yp*densit/viscos
	   te(i,j)=utaw**2/cdterm
	 endif
      enddo
c-----south wall left
      cdterm=sqrt(cmu)
      yp=y(jst)-yv(jst)
      j=jst
      do i=2,istm1
	 up=(u(i,j)+u(i+1,j))/2
	 ueff=abs(up)
	 utaw=sqrt(ueff*viscos/yp/densit)
	 ypluss(i)=utaw*yp*densit/viscos
	 if(ypluss(i).le.10.8) then
	   te(i,j)=utaw**2/cdterm*(ypluss(i)/10.8)**2
	 else
	   root=cappa*ueff
	   if(iter.gt.100) then
	     call newton(yp,ueff,root)
	   endif
	   utaw=root
	   ypluss(i)=utaw*yp*densit/viscos
	   te(i,j)=utaw**2/cdterm
	 endif
      enddo
c-----south wall right
      cdterm=sqrt(cmu)
      yp=y(2)
      j=2
      do i=ist,nim1
	 up=(u(i,j)+u(i+1,j))/2
	 ueff=abs(up)
	 utaw=sqrt(ueff*viscos/yp/densit)
	 ypluss(i)=utaw*yp*densit/viscos
	 if(ypluss(i).le.10.8) then
	   te(i,j)=utaw**2/cdterm*(ypluss(i)/10.8)**2
	 else
	   root=cappa*ueff
	   if(iter.gt.100) then
	     call newton(yp,ueff,root)
	   endif
	   utaw=root
	   ypluss(i)=utaw*yp*densit/viscos
	   te(i,j)=utaw**2/cdterm
	 endif
      enddo
c-----west wall
      cdterm=sqrt(cmu)
      xp=x(ist)
      i=ist
      do j=2,jstm1
	 vp=(v(i,j)+v(i,j+1))/2
	 veff=abs(vp)
	 vtaw=sqrt(veff*viscos/xp/densit)
	 xplusw(j)=vtaw*xp*densit/viscos
	 if(xplusw(j).le.10.8) then
	   te(i,j)=vtaw**2/cdterm*(xplusw(j)/10.8)**2
	 else
	   root=cappa*veff
	   if(iter.gt.100) then
	     call newton(xp,veff,root)
	   endif
	   vtaw=root
	   xplusw(j)=vtaw*xp*densit/viscos
	   te(i,j)=vtaw**2/cdterm
	 endif
      enddo
c-----outlet
      do j=2,njm1
      te(ni,j)=te(nim1,j)+dxep(nim1)/dxpw(nim1)*(te(nim1,j)-te(ni-2,j))
      enddo
      return
c     ---------------------------------
c            d-trans modification  
c     ---------------------------------
      entry moded
c-----top boundry
      yp=yv(nj)-y(njm1)
      j=njm1
      term=(cmu**0.75)/(cappa*yp)
      do i=2,nim1
       if(yplusn(i).le.10.8) then
	 ret=sqrt(te(i,j))*yp*densit/viscos
	 ed(i,j)=term*te(i,j)**1.5*(1+5.3/ret)
       else
	 ed(i,j)=term*te(i,j)**1.5
       endif
      enddo
c-----south wall left
      yp=y(jst)-yv(jst)
      j=jst
      term=(cmu**0.75)/(cappa*yp)
      do i=2,istm1
       if(ypluss(i).le.10.8) then
	 ret=sqrt(te(i,j))*yp*densit/viscos
	 ed(i,j)=term*te(i,j)**1.5*(1+5.3/ret)
       else
	 ed(i,j)=term*te(i,j)**1.5
       endif
      enddo
c-----south wall right
      yp=y(2)
      j=2
      term=(cmu**0.75)/(cappa*yp)
      do i=ist,nim1
       if(ypluss(i).le.10.8) then
	 ret=sqrt(te(i,j))*yp*densit/viscos
	 ed(i,j)=term*te(i,j)**1.5*(1+5.3/ret)
       else
	 ed(i,j)=term*te(i,j)**1.5
       endif
      enddo
c-----west wall
      xp=x(ist)
      i=ist
      term=(cmu**.75)/(cappa*xp)
      do j=2,jstm1
       if(xplusw(j).le.10.8) then
	 ret=sqrt(te(i,j))*xp*densit/viscos
	 ed(i,j)=term*te(i,j)**1.5*(1+5.3/ret)
       else
	 ed(i,j)=term*te(i,j)**1.5
       endif
      enddo
c-----outlet 
      do j=2,njm1
      ed(ni,j)=ed(nim1,j)+dxep(nim1)/dxpw(nim1)*(ed(nim1,j)-ed(ni-2,j))
      enddo
      return
      end
c     ---------------------------------------------------------------
c                           subroutine newton
c     ---------------------------------------------------------------         
      subroutine newton(yp,ueff,root)     
      include 'wall.inc'
      dnu=viscos/densit
      x1=root
10    f=alog(elog*yp*x1/dnu)/cappa
      x2=ueff/f
      dx=abs(x2-x1)
      if(dx.gt.0.01) then
	x1=x2
	goto 10
      else
	root=x2
      endif
      return
      end
c     ---------------------------------------------------------------
c                           subroutine output
c     ---------------------------------------------------------------         
      subroutine output
      include 'wall.inc'
c     .................

c-----tecplot format
      write(1,*)'variables=x,y,u,v,p,k,d,re'
      write(1,10)ni-2,nj-2
10    format('zone',3x,'i=',i3,3x,'j=',i3)   
c-----interpolating      
      do j=2,njm1
	 do i=2,nim1
	    u(i,j)=(u(i,j)+u(i+1,j))/2
	    v(i,j)=(v(i,j)+v(i,j+1))/2
            vis(i,j)=densit*te(i,j)**2*cmu/ed(i,j)
            uv=sqrt(u(i,j)**2+v(i,j)**2)
            celre(i,j)=densit*uv*0.5*(sns(j)+sew(i))/vis(i,j)
	 enddo
      enddo
      do j=2,njm1
	 do i=2,nim1
	    write(1,20) x(i),y(j),u(i,j),v(i,j),p(i,j),
     &                      te(i,j),ed(i,j),celre(i,j)
	 enddo
      enddo
20    format(8(3x,e15.7))
c-----wall function
      write(2,*)'top wall'
      do i=2,nim1
	 write(2,*)yplusn(i)
      enddo
      write(2,*)'south wall on left'
      do i=2,istm1
	 write(2,*)ypluss(i)
      enddo
      write(2,*)'south wall on right'
      do i=ist,nim1
	 write(2,*)ypluss(i)
      enddo
      write(2,*)'west wall'
      do j=2,jstm1
	write(2,*)xplusw(j)
      enddo   
      return
      end
c     ---------------------------------------------------------------
c                           subroutine calc
c     ---------------------------------------------------------------
      subroutine calc
      include 'wall.inc'
c     .................
c      if(abs(phic-phiu).ge.0.5*abs(phid-phiu)) then
c         phiho=phic+rf/ru*(phic-phiu)
        c=phic
        a=(ru*phid+rd*phiu-(ru+rd)*c)/(rd**2*ru+ru**2*rd)
        b=(phid-a*rd**2-c)/rd
        phiho=a*rf**2+b*rf+c
c      elseif(abs(phic-phiu).lt.0.5*abs(phid-phiu)) then
c         phiho=rf/rd*phid+(1-rf/rd)*phic
c      endif
      p1=phic
      p2=phid
      p3=phiu+3*(phic-phiu)
      if(p1.Gt.p2.and.p1.Gt.p3) then
	  p4=dmax1(dmin1(p1,p2),dmin1(p1,p3))
      else
	  p4=dmin1(dmax1(p1,p2),dmax1(p1,p3))
      endif 
      if(p1.Gt.p4.and.p1.Gt.phiho) then
	  phiho=dmax1(dmin1(p1,p4),dmin1(p1,phiho))
      else 
	  phiho=dmin1(dmax1(p1,p4),dmax1(p1,phiho)) 
      endif
       return
       end
c     ---------------------------------------------------------------
c                           END OF PROGRAM
c     ---------------------------------------------------------------   

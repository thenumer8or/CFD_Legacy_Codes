      implicit double precision(a-h,o-z)
      parameter(ni=52,nj=42)
      common
     1/uvel/resoru,nswpu,urfu,dxepu(ni),dxpwu(ni),sewu(ni)
     1/vvel/resorv,nswpv,urfv,dynpv(nj),dypsv(nj),snsv(nj),rcv(nj)
     1/pcor/resorm,nswpp,urfp,du(ni,nj),dv(ni,nj),ipref,jpref
     1/ten/resork,nswpk,urfk
     1/tdis/resore,nswpd,urfe
     1/var/u(ni,nj),v(ni,nj),p(ni,nj),pp(ni,nj),te(ni,nj),ed(ni,nj)
     1/all/nim1,njm1,great,ihp
     1/geom/indcos,x(ni),y(nj),dxep(ni),dxpw(ni),dynp(nj),dyps(nj),
     2      sns(nj),sew(ni),xu(ni),yv(nj),r(nj),rv(nj)
     1/flupr/urfvis,viscos,densit,prandt,den(ni,nj),vis(ni,nj)
     1/kase t1/uin,tein,edin,flowin
     1/divid/jst,jstp1,jstm1,ista,istap1,istam1,istb,istbp1,istbm1
     1/turb/gen(ni,nj),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/wallf/ypluss(ni),xplusw(nj),xpluse(nj),
     2       taus(ni),tauw(nj),taue(nj)
     1/coef/ap(ni,nj),an(ni,nj),as(ni,nj),ae(ni,nj),aw(ni,nj),
     2      su(ni,nj),sp(ni,nj),sukd(ni,nj),spkd(ni,nj)                                             

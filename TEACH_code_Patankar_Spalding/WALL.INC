      parameter(ni=90,nj=28,nim1=ni-1,njm1=nj-1)
      implicit double precision(a-h,o-z)
      common
     1/uvel/xu(ni),dxepu(ni),dxpwu(ni),sewu(ni)
     2/vvel/yv(nj),dynpv(nj),dypsv(nj),snsv(nj)
     3/pcor/du(ni,nj),dv(ni,nj)
     4/urf/urfu,urfv,urfp,urfk,urfd,urfvis
     5/sor/resoru,resorv,resorm,resork,resord
     6/var/u(ni,nj),v(ni,nj),p(ni,nj),pp(ni,nj),te(ni,nj),ed(ni,nj)
     7/swee/nswpu,nswpv,nswpp,nswpk,nswpd
     8/all/indu,indbu,indv,indbv,indk,indbk,indd,indbd
     8/sss/iter,maxit,ipref,jpref
     9/geo/h,w1,w2,dl,dr
     1/cont/great,sormax
     1/geom/x(ni),y(nj),dxep(ni),dxpw(ni),dynp(nj),dyps(nj),
     2      sns(nj),sew(ni)
     1/flupr/viscos,densit,prandt,vis(ni,nj),celre(ni,nj)
     1/kase/uin,tein,edin,flowin
     1/divid/jst,jstp1,jstm1,ist,istp1,istm1
     1/turb/gen(ni,nj),cd,cmu,c1,c2,cappa,elog,pred,prte
     1/wallf/yplusn(ni),ypluss(ni),xplusw(nj)
     1/coef/ap(ni,nj),an(ni,nj),as(ni,nj),ae(ni,nj),aw(ni,nj),
     2      su(ni,nj),sp(ni,nj)                                             
     1/hor/phic,phid,phiu,phiho,ru,rd,rf
     

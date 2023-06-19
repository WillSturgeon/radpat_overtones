      subroutine forward_model_mineos(
     1 phvel_all,grvel_all,lorder_all,attn_all,per_all,
     1 jcomin,epsin,wgravin,lminin,lmaxin,wminin,wmaxin,nminin,nmaxin,
     1 model_file,outputs_dir,premnm)

      implicit real*8(a-h,o-z)

      integer*4 mk
      parameter (mk=650)
      integer*4 pk
      parameter (pk=450)
      integer*4 wk
      parameter (wk=3000)

      character*256  model_file,outputs_dir
      character*256  out_plain_file,out_bin_file
      character*256  dbase_name,eigenasc,kernelasc
      real*4      rad(mk)
      real*4      omega
      real*4      kkappa(mk),kmu(mk)
      real*8      alpha(mk),beta(mk)
      real*4      kalpha(mk),kbeta(mk)
      real*4      rhobar,bigg,tau
      real*4      fl
      real*4      phvel_all(pk),grvel_all(pk)
      real*4      attn_all(pk),per_all(pk)
      integer*4   lorder_all(pk)
      integer*4   jcomin,wgravin,lminin,lmaxin
      integer*4   wminin,wmaxin,nminin,nmaxin,premnm,pm
      real*4      epsin

      real*4    per_eigen,phvel_eigen,grvel_eigen,attn_eigen
      integer*4 norder_eigen,lorder_eigen,eigid_eigen,
     +          nraw_eigen,ncol_eigen,npar_eigen,foff_eigen,
     +          commid_eigen
      character*2 datatype_eigen
      character*64 dir_eigen
      character*32 dfile_eigen
      character*17 lddate_eigen
      character*1 typeo_eigen
      character*35 path
      real*4 buf
      common/c_eigen/norder_eigen,lorder_eigen,
     +      eigid_eigen,per_eigen,phvel_eigen,grvel_eigen,
     +      attn_eigen,nraw_eigen,ncol_eigen,npar_eigen,
     +      foff_eigen,commid_eigen,typeo_eigen,
     +      datatype_eigen,dir_eigen,dfile_eigen,lddate_eigen

      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/vpv(mk),vph(mk),vsv(mk),vsh(mk),eta(mk),wrk(mk*10)
      common/will/cvel,wmhz,tcom,gcom,qmod,cvel_all(wk),wmhz_all(wk),
     +  tcom_all(wk),gcom_all(wk),qmod_all(wk)
      common/c_buf/nn,ll,ww,qq,gc,buf(6*mk)
      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/

      jcom=jcomin
      if (jcom.lt.0.or.jcom.gt.5) then 
        print*,"Invalid jcom"
      endif
      eps=epsin
      wgrav=wgravin
      lmin=lminin
      lmax=lmaxin
      wmin=wminin
      wmax=wmaxin
      nmin=nminin
      nmax=nmaxin


      open(7,file=model_file,status='old',form='formatted',iostat=iret)
      call model(7,8) 
      close(7)
      alpha=vpv/tau
      beta=vsv/tau
      ifreq=1
      call wtable(8,3,ifreq,lmin,lmax,wmin,wmax,nmin,nmax)

      pm=premnm
      do n=nmin,nmax
       do l=lmin,lmax
            do i=1,pm
            if (jcomin==3) then
            buf((pm*2)+1+i)=buf((pm*2)+1+i)/sqrt(l*(l+1.)) 
            buf((pm*3)+1+i)=buf((pm*3)+1+i)/sqrt(l*(l+1.))
            elseif (jcomin==2) then
            buf(i)=buf(i)/sqrt(l*(l+1.)) 
            buf(pm+1+i)=buf(pm+1+i)/sqrt(l*(l+1.)) 
            endif
            enddo
       enddo
      enddo

      end subroutine
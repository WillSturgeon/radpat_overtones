      program main
    
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
      integer*4   wminin,wmaxin,nminin,nmaxin
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
      real*4 abuf,buf
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
      common/c_buf/nn,ll,ww,qq,gc,buf(6*mk)
      common/will/cvel,wmhz,tcom,gcom,qmod,cvel_all(wk),wmhz_all(wk),
     +  tcom_all(wk),gcom_all(wk),qmod_all(wk)
      common/eifx/a(14,mk),dum(mk)
      dimension abuf(6*mk+5)
      equivalence (nn,abuf)

      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/
c---- declarations for rad_pat
      real*4 lambda, sdpth, azimuth
      real*4 part1,part3,part3_az
      dimension fmom(6),azep_i(360)
      integer alpha_rp
      complex part2,MEs(360),part2_az,MEs_az,MEs_azepin
      complex MEs_azepin_plus1,MEs_azepin_minus1
      real*8 X(40),Y(40),U,V,B(40),C(40),D(40)
      real*8 Ueigen_sdpth,Udoteigen_sdpth,Veigen_sdpth,Vdoteigen_sdpth
      real*8 Weigen_sdpth,Wdoteigen_sdpth
      real*8, dimension(:), allocatable :: Ueigen,Udoteigen   
      real*8, dimension(:), allocatable :: Veigen,Vdoteigen
      real*8, dimension(:), allocatable :: Weigen,Wdoteigen

c-----------------------------------------------------------------------
c---- declarations for model inputs etc
      character*40    premstr
      integer         premifanis,premifdeck,premnm,premnoc
      real*4          premnic
      parameter       (maxmod=700)
      real*8          prrad(maxmod), prrho(maxmod)
      real*8          prvpv(maxmod), prvsv(maxmod)
      real*8          prqka(maxmod), prqmu(maxmod)
      real*8          prvph(maxmod), prvsh(maxmod)
      real*8          preta(maxmod)
      integer         i,mm,k
      real*4          moho, mineos_moho
c      real*4          vs(maxlayer)
      real*8          radius(300)
      character*256   path1,modelin,path2,path3,path4

c---- declarations for GCMT catalogue 
      parameter (max_events=49526)
      character*16 evnam,evnam2
      character*4 hypo(max_events)
      character*10 refdat(max_events)
      character*10 reftime(max_events)
      real*4 lat(max_events)
      real*4 long(max_events)
      real*4 deptha(max_events)
      real*4 mb(max_events)
      real*4 ms(max_events)
      character*24 geoloc(max_events)
      character*16 CMTname(max_events),dd(max_events)
      character*1 CMTpreletter(max_events)
      character*58 centroidparam(max_events)
      integer exponent(max_events),exponentx
      real*4 Mrr(max_events),Mtt(max_events),Mpp(max_events)
      real*4 Mrt(max_events),Mrp(max_events),Mtp(max_events)
      real*4 Mrr_err(max_events),Mtt_err(max_events),Mpp_err(max_events)
      real*4 Mrt_err(max_events),Mrp_err(max_events),Mtp_err(max_events)
      character*3 version(max_events)
      real*4 Mmm(max_events),Mmm_form(max_events)
      real*4 xm(6),rp_max_azimuth(1)
c      real*4 eventt(401501),network(401501),station(401501),az_earth(401501)
      integer CMT_mm(max_events),CMT_dd(max_events),CMT_yy(max_events)
      character*1 CMTpostletter(max_events),slash(max_events)
      character*2 datQQ(max_events),datYY(max_events),datMM(max_events)
      character*2 datDD(max_events),timeHH(max_events)
      character*2 timeMM(max_events),timeSS(max_events)
      character*2 evMM(1),evDD(1),evYY(1)
      character*1 evLETTER(1),timeMS(max_events)
      character*6 ndkdata(max_events),evformat(1)
      character*6 ndkdata_MMDDYY(max_events),evformat_MMDDYY(1)
      character*7 evformatLETTER(1)
      character*16 short_CMTname(max_events)
      logical foundevent
      character*50 input_model,outdir
      character*1 jcominarg,nmaxinarg
      character*3 lmaxinarg
      character*10 azep_inarg
      real*8 azep_in
      character*6 station
      character*4 net
      character*10 evlat,evlon,stlat,stlon,Ampfl,error,mima,nsim
      integer azep_inx,linecount

c-------------------------------------------------------
c-------------------- inputs ---------------------------
c-- example: ./radpat_overtones prem_noocean 3 0 098 202201281114A 70 /data/will/rad_pat/
c First, make sure the right number of inputs have been provided
      IF(COMMAND_ARGUMENT_COUNT().NE.18)THEN
      WRITE(*,*)'ERROR, INCORRECT N. INPUT ARGUMENTS, STOPPING'
      STOP
      ENDIF

      CALL GET_COMMAND_ARGUMENT(1,input_model)
      CALL GET_COMMAND_ARGUMENT(2,jcominarg)
      CALL GET_COMMAND_ARGUMENT(3,nmaxinarg)
      CALL GET_COMMAND_ARGUMENT(4,lmaxinarg)
      CALL GET_COMMAND_ARGUMENT(5,evnam)
      CALL GET_COMMAND_ARGUMENT(6,azep_inarg)
      CALL GET_COMMAND_ARGUMENT(7,outdir)
      CALL GET_COMMAND_ARGUMENT(8,evnam2)
      CALL GET_COMMAND_ARGUMENT(9,station)
      CALL GET_COMMAND_ARGUMENT(10,net)
      CALL GET_COMMAND_ARGUMENT(11,stlat)
      CALL GET_COMMAND_ARGUMENT(12,stlon)
      CALL GET_COMMAND_ARGUMENT(13,evlat)
      CALL GET_COMMAND_ARGUMENT(14,evlon)
      CALL GET_COMMAND_ARGUMENT(15,Ampfl)
      CALL GET_COMMAND_ARGUMENT(16,error)
      CALL GET_COMMAND_ARGUMENT(17,mima)
      CALL GET_COMMAND_ARGUMENT(18,nsim)

      read(jcominarg,*)jcomin
      read(nmaxinarg,*)nmaxin
      read(lmaxinarg,*)lmaxin
      read(azep_inarg,*)azep_in

c      write(*,*)'input model = ',input_model
c      write(*,*)'jcom = ',jcomin
c      write(*,*)'nmax= ',nmaxin
c      write(*,*)'lmax = ',lmaxin
c      write(*,*)'evnam = ',evnam
c      write(*,*)'input azimuth = ',azep_in
c      write(*,*)'output directory =',outdir
c      write(*,*)'station = ',station
c      write(*,*)'net = ',net
c      write(*,*)'stlat = ',stlat
c      write(*,*)'stlon = ',stlon
c      write(*,*)'evlat = ',evlat
c      write(*,*)'evlon = ',evlon
c      write(*,*)'A/A0 = ',Ampfl
c      write(*,*)'error = ',error
c      write(*,*)'mima =',mima
c      write(*,*)'nsim = ',nsim
      
c-----  load input model

      model_file=trim('../models/'//trim(input_model))
      modelin=trim(input_model)

      linecount=0

      open(112,file=model_file,status="old",access="sequential")
      read(112,'(a)') premstr
      read(112,*) premifanis, premtref, premifdeck
      read(112,*) premnm, premnic, premnoc   
          do i=1,premnm 
          linecount=linecount+1
          read(112,201)prrad(i),prrho(i),prvpv(i),prvsv(i),
     1 prqka(i),prqmu(i),prvph(i),prvsh(i),preta(i)
          end do
      close(112)

        premnm=linecount

      open(112,file=model_file,status="old",access="sequential")
      read(112,'(a)') premstr
      read(112,*) premifanis, premtref, premifdeck
      read(112,*) premnm, premnic, premnoc   
          do i=1,premnm 
          read(112,201)prrad(i),prrho(i),prvpv(i),prvsv(i),
     1 prqka(i),prqmu(i),prvph(i),prvsh(i),preta(i)
          end do
      close(112)

  201 format(f8.0,3f9.2,2f9.1,2f9.2,f9.5)
c-----------------------------------------------------------------------
c------------------------- input parameters ----------------------------
c-----------------------------------------------------------------------

c ------- jcom - 1=radial, 2=toroidal, 3=spheroidal, 4=inner core toroidal
c ------- eps - 10−7 for periods > 10 s. 10−12 − 10−10 for periods between 5-10 s
c ------- wgrav - frequency in millihertz (mHz) above which gravitational terms are neglected; this gives about a factor of 3 increase in speed.

      jcomin=jcomin
      epsin=1e-7
      wgravin=10
      lminin=lmaxin
      lmaxin=lmaxin
      wminin=0
      wmaxin=166.0
      nminin=nmaxin
      nmaxin=nmaxin
c-----------------------------------------------------------------------

      call forward_model_mineos(
     1 phvel_all,grvel_all,lorder_all,attn_all,per_all,jcomin,epsin,
     1 wgravin,lminin,lmaxin,wminin,wmaxin,nminin,nmaxin,
     1 model_file,outputs_dir,premnm)

c-----------------------------------------------------------------------
cc open new GCMT file and read (WS)
      open(201,file='../data/jan76_dec17.ndk')
       do i=1,max_events
c        do i=1,20
         read(201,504)hypo(i),datQQ(i),datYY(i),slash(i),datMM(i),
     1   slash(i),datDD(i),timeHH(i),slash(i),timeMM(i),slash(i),
     1   timeSS(i),slash(i),timeMS(i),lat(i),long(i),deptha(i),mb(i),
     1   ms(i),geoloc(1)
c         write(*,504)hypo(i),datQQ(i),datYY(i),slash(i),datMM(i),
c     1   slash(i),datDD(i),timeHH(i),slash(i),timeMM(i),slash(i),
c     1   timeSS(i),slash(i),timeMS(i),lat(i),long(i),deptha(i),mb(i),
c     1   ms(i),geoloc(1)       
         read(201,503)CMTpreletter(i),CMTname(i)
c         write(*,*)CMTpreletter(i)," ",CMTname(i)
         read(201,*)centroidparam(i)
c         write(*,*)centroidparam(i)
         read(201,502) exponent(i),Mrr(i),Mrr_err(i),Mtt(i),Mtt_err(i),
     1   Mpp(i),Mpp_err(i),Mrt(i),Mrt_err(i),Mrp(i),Mrp_err(i),Mtp(i),
     1   Mtp_err(i)
c         write(*,502) exponent(i),Mrr(i),Mrr_err(i),Mtt(i),Mtt_err(i),
c     1   Mpp(i),Mpp_err(i),Mrt(i),Mrt_err(i),Mrp(i),Mrp_err(i),Mtp(i),
c     1   Mtp_err(i)
         read(201,*)version(i)
c         write(*,*)version(i)
       enddo
      close(201)

  501 format(A4,1x,A10,1x,A10,1x,f6.2,1x,f7.2,1x,f5.1,1x,
     1 f3.1,1x,f3.1,1x,A26)
  502 format(i2.0,1x,f6.3,1x,f5.3,1x,f6.3,1x,f5.3,
     1 1x,f6.3,1x,f5.3,1x,f6.3,1x,f5.3,1x,f6.3,1x,
     1 f5.3,1x,f6.3,1x,f5.3)
  503 format(A1,A13) 
c above A7 should be A13 for the whole name, but A7 just picks out the date, the rest includes the time
c  504 format(A4,1x,i2,i2,A1,i2,A1,i2,1x,i2,A1,i2,A1,i2,A1,i1,1x,f6.2,1x,f7.2,1x,f5.1,1x,f3.1,1x,f3.1,1x,A26)    
  504 format(A4,1x,A2,A2,A1,A2,A1,A2,1x,A2,A1,A2,A1,A2,A1,
     1 A1,1x,f6.2,1x,f7.2,1x,f5.1,1x,f3.1,1x,f3.1,1x,A26)  
      
C reads event name from the standard input:

cc splitting up the input event name, translating it to the format is YY/MM/DD   
      evMM=evnam(1:2)
      evDD=evnam(3:4)
      evYY=evnam(5:6)
      evLETTER=evnam(7:7)
      evformat=evYY//evMM//evDD
      evformat_MMDDYY=evMM//evDD//evYY
      evformatLETTER=evYY//evMM//evDD//evLETTER   

c cc forming date from ndk file, in same format as input event
      do i=1,max_events
      ndkdata(i)=datYY(i)//datMM(i)//datDD(i)
      ndkdata_MMDDYY(i)=datMM(i)//datDD(i)//datYY(i)
c      write(*,*) 'ndk YY ',datYY(i),' ndk MM ',datMM(i),' ndk DD ',datDD(i),' ndkdata ',ndkdata(i),' nkdata_v2 ',ndkdata_MMDDYY(i)
      end do

cc ensure the evnam (from command line) is equal to the CMTname from the input GCMT file (WS)
      do i=1,max_events
c      write(*,*) '---- ',CMTname(i), evnam
c      write(*,*) "start of loop CMTname(i) ",CMTname(i)," ndkdata(i) ",ndkdata(i)," ndkdata_MMDDYY ",ndkdata_MMDDYY(i)," evnam ",evnam," evformat ",evformat," evformat_MMDDYY ", evformat_MMDDYY

        if (trim(evnam).eq.trim(CMTname(i))) then 
        j=i
c        write(*,*) 'condition found - exit'
        evnam=CMTname(j)
        foundevent=.true.
        exit
        else
        foundevent=.false.
        end if
      end do
      
c        write(*,*) 'after loop 1 evnam ',evnam
      
      if (foundevent.eqv..false.) then

      do i=1,max_events
        if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'A') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','A'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'a') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','a'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
         j=i
         exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'B') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','B'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'b') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','b'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'C') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','C'
c        evnam=CMTname(j)
c        write(*,*) '++++++ ',evnam, CMTname(j)
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'c') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','c'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq."D") then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','D'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit
  
        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'d') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','d'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'E') then 
c        write(*,*) 'E====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','E'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam, CMTname(i)
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'e') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','e'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'F') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','F'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'f') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','f'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'G') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','G'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'g') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','g'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'H') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','H'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'h') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','h'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'I') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','I'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'i') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','i'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'J') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','J'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'j') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','j'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'K') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','K'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i     
        exit                               

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'k') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','k'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'L') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','L'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i   
        exit                                   

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'l') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','l'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'M') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','M'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'m') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','m'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'N') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','N'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'n') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','n'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'O') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','O'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'o') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','o'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'P') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','P'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'p') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','p'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'Q') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','Q'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'q') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','q'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'R') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','R'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'r') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','r'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'S') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','S'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'s') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','s'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'T') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','T'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'t') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','t'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'U') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','U'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'u') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','u'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'V') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','V'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'v') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','v'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'W') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','W'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'w') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','w'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'X') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','X'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'x') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','x'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'Y') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','Y'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'y') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','y'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'Z') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','Z'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata(i).eq.evformat(1).AND.evLETTER(1).eq.'z') then 
c        write(*,*) '====== ',ndkdata(i),' ',evformat(1),' ',evLETTER(1),' ','z'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

c-----------------------------------------------------------------------------------------------------------------------------------

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'A') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','A'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'a') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','a'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'B') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','B'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'b') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','b'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'C') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','C'
c        evnam=CMTname(j)
c        write(*,*) '++++++ ',evnam, CMTname(j)
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'c') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','c'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq."D") then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','D'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit
  
        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1 AND.evLETTER(1).eq.'d') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','d'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'E') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','E'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam, CMTname(i)
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'e') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','e'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'F') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','F'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat_MMDDYY(1).
     1  AND.evLETTER(1).eq.'f') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat_MMDDYY(1),' ',evLETTER(1),' ','f'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'G') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','G'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'g') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','g'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'H') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','H'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'h') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','h'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'I') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','I'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'i') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','i'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'J') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','J'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'j') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','j'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'K') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','K'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i 
        exit                                    

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'k') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','k'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'L') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','L'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i      
        exit                                

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'l') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','l'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'M') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','M'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
        exit

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'m') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','m'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'N') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','N'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'n') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','n'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'O') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','O'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'o') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','o'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'P') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','P'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'p') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','p'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'Q') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','Q'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'q') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','q'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'R') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','R'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'r') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','r'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'S') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','S'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'s') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','s'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'T') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','T'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'t') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','t'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'U') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','U'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'u') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','u'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'V') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','V'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'v') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','v'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'W') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','W'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'w') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','w'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'X') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','X'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'x') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','x'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'Y') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','Y'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'y') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','y'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'Z') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','Z'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i

        else if (ndkdata_MMDDYY(i).eq.evformat(1).
     1  AND.evLETTER(1).eq.'z') then 
c        write(*,*) '====== ',ndkdata_MMDDYY(i),' ',evformat(1),' ',evLETTER(1),' ','z'
c        evnam=CMTname(i)
c        write(*,*) '++++++ ',evnam
        j=i
c-----------------------------------------------------------------------------------------------------------------------------------        
         end if
       end do
       end if

      evnam=CMTname(j)
      write(*,*)"the CMT event name is ",evnam,j
      sdpth=deptha(j)
      exponentx=exponent(j)

      xm(1)=Mrr(j)
      xm(2)=Mtt(j)
      xm(3)=Mpp(j)
      xm(4)=Mrt(j)
      xm(5)=Mrp(j)
      xm(6)=Mtp(j)
c the order is Mrr, Mtt, Mpp, Mrt, Mrp, Mtp

c      write(*,*) 'evnam new ',evnam
c      write(*,*) 'source depth ',sdpth
c      write(*,*)'exponent ',exponentx

c eigenfunctions at the surface - will need to change this to the correct source depth
c I had to multiply by -1 so that is matches the modes in the mode files....

      if (jcomin==3) then
            Ueigen=buf(1:premnm)*(-1)
            Udoteigen=buf(premnm+1:premnm*2)*(-1)
            Veigen=buf(premnm*2+1:premnm*3)*(-1)
            Vdoteigen=buf(premnm*3+1:premnm*4)*(-1)
c            Peigen=buf(premnm*5)*(-1)
c            Pdoteigen=buf(premnm*6)*(-1)
c            write(*,*)'U --- ',Ueigen
c            write(*,*)'Udot --- ',Udoteigen
c            write(*,*)'V --- ',Veigen
c            write(*,*)'Vdot --- ',Vdoteigen
c            write(*,*)'P --- ',Peigen
c            write(*,*)'Pdot --- ',Pdoteigen
      elseif (jcomin==2) then
            Weigen=buf(1:premnm)*(-1)
            Wdoteigen=buf(premnm+1:premnm*2)*(-1)
c            write(*,*)'W --- ',Weigen
c            write(*,*)"Wdot --- ",Wdoteigen
      end if
      
c tweaking radius to ensure there are no repeated values for the interpolation
      do i=premnm-40,premnm
      if (prrad(i)==prrad(i+1)) then
c      write(*,*) i,'radius duplicate'
      prrad(i)=prrad(i)-1
      endif
      enddo 

c---- check is sdpth already exists in input model
      do j=1,premnm
      if (6371000-(sdpth*1000) == prrad(j)) then
c      write(*,*) 'source depth exists in input model'
      if (jcomin==3) then 
      Ueigen_sdpth=Ueigen(j)
      Udoteigen_sdpth=Udoteigen(j)
      Veigen_sdpth=Veigen(j)
      Vdoteigen_sdpth=Vdoteigen(j)
c      write(*,*) '1.Ueigen_sdpth ',Ueigen_sdpth
c      write(*,*) '1.Udoteigen_sdpth ',Udoteigen_sdpth
c      write(*,*) '1.Veigen_sdpth ',Veigen_sdpth
c      write(*,*) '1.Vdoteigen_sdpth ',Vdoteigen_sdpth
      elseif (jcomin==2) then
      Weigen_sdpth=Weigen(j)
      Wdoteigen_sdpth=Wdoteigen(j)
c      write(*,*) '1.Weigen_sdpth ',Weigen_sdpth
c      write(*,*) '1.Wdoteigen_sdpth ',Wdoteigen_sdpth
      endif
      goto 11
      else
c      write(*,*)'need to interpolate'
      endif 
      enddo

c else perform linear interepolation between the layers above and below
      do i=1,premnm
        if (prrad(i)<=(6371000-(sdpth*1000))) then
c        write(*,*) (prrad(i)), 6371000-(sdpth*1000)
        nlines=nlines+1
        endif
      enddo

      if (jcomin==3) then

       Ueigen_sdpth=(Ueigen(nlines)*(prrad(nlines+1)
     1 -(6371000-(sdpth*1000)))+(Ueigen(nlines+1)
     1 *((6371000-(sdpth*1000))-prrad(nlines))))/
     1 (prrad(nlines+1)-prrad(nlines))

       Udoteigen_sdpth=(Udoteigen(nlines)*(prrad(nlines+1)
     1 -(6371000-(sdpth*1000)))+(Udoteigen(nlines+1)
     1 *((6371000-(sdpth*1000))-prrad(nlines))))/
     1 (prrad(nlines+1)-prrad(nlines))

       Veigen_sdpth=(Veigen(nlines)*(prrad(nlines+1)
     1 -(6371000-(sdpth*1000)))+(Veigen(nlines+1)
     1 *((6371000-(sdpth*1000))-prrad(nlines))))/
     1 (prrad(nlines+1)-prrad(nlines))

       Vdoteigen_sdpth=(Vdoteigen(nlines)*(prrad(nlines+1)
     1 -(6371000-(sdpth*1000)))+(Vdoteigen(nlines+1)
     1 *((6371000-(sdpth*1000))-prrad(nlines))))/
     1 (prrad(nlines+1)-prrad(nlines))

c      write(*,*)'2.Ueigen_spdth ',Ueigen_sdpth
c      write(*,*)'2.Udoteigen_spdth ',Udoteigen_sdpth
c      write(*,*)'2.Veigen_spdth ',Veigen_sdpth
c      write(*,*)'2.Vdoteigen_spdth ',Vdoteigen_sdpth

      elseif (jcomin==2) then
        Weigen_sdpth=(Weigen(nlines)*(prrad(nlines+1)
     1 -(6371000-(sdpth*1000)))+(Weigen(nlines+1)
     1 *((6371000-(sdpth*1000))-prrad(nlines))))/
     1 (prrad(nlines+1)-prrad(nlines))

       Wdoteigen_sdpth=(Wdoteigen(nlines)*(prrad(nlines+1)
     1 -(6371000-(sdpth*1000)))+(Wdoteigen(nlines+1)
     1 *((6371000-(sdpth*1000))-prrad(nlines))))/
     1 (prrad(nlines+1)-prrad(nlines))

c      write(*,*)'2.Weigen_spdth ',Weigen_sdpth
c      write(*,*)'2.Wdoteigen_spdth ',Wdoteigen_sdpth
      endif

   11 continue

      alpha_rp=0 ! surface wave orbit 
      a=6371.0 ! earth radius 
      lambda=lmaxin+0.5 ! a=1
      aziumth=0

C Moment tensor components !why -30???? ask Ana.
      do i=1,6
      fmom(i)=xm(i)*(10.0**(exponentx-30)) 
c      write(*,*) 'Ms ',fmom(i)
      enddo

c------------ Rayleigh waves
c M:Es* (MEs) - Table A1 - Ferreira & Woodhouse 2006.
      if (jcomin==3) then
      azep_inx=NINT(azep_in)
      do i=1,360

      azep_i(i)=i*pi/180.
      azep=azep_i(i)

      part1=-((Udoteigen_sdpth*fmom(1))+
     1 (0.5*((2*Ueigen_sdpth)-((lambda**2)*Veigen_sdpth)))
     1 *(fmom(2)+fmom(3)))

      part2=-cmplx(0,((lambda)*((fmom(4)*cos(pi-azep))
     1 +(fmom(5)*sin(pi-azep)))))*
     1 (Vdoteigen_sdpth+(Ueigen_sdpth-Veigen_sdpth))

      part3=((lambda**2)*Veigen_sdpth)*
     1 ((0.5*(fmom(2)-fmom(3))*
     1 cos(2*(pi-azep))) +(fmom(6)*sin(2*(pi-azep))))

      MEs(i)=abs(part1+part2+part3)

      enddo

       MEs_azepin=MEs(azep_inx)
       MEs_azepin_plus1=MEs(azep_inx+1)
       MEs_azepin_minus1=MEs(azep_inx-1)

c       write(*,*)'MEs in ',abs(MEs_azepin)
c       write(*,*)'MEs in+1 ',abs(MEs_azepin_plus1)
c       write(*,*)'MEs in-1 ',abs(MEs_azepin_minus1)

c------ write final file output Rayleigh

c      write(path1,'(A,A,A,A,A,I1.1,A,I3.3,A)') trim(outdir),
c     1 trim(evnam),"_",trim(input_model),"_n",nmaxin,"_l",
c     1 lmaxin,"_Rayleigh.txt"

c      open(1,file=path1,status='unknown',access='sequential',
c     1  position='append')
c      do i=1,360
c       azep_i(i)=i*pi/180.
c       azep=azep_i(i)

c --- this was (pi-azep)*180/pi - not sure why exactly
c       write(1,*) azep*(180/pi), abs(MEs(i))*1000 

c      enddo
c      close(1)

      write(path3,'(A,I1.1,A,I3.3,A)') trim(outdir),
     1 nmaxin,"S",lmaxin,"_nodality.txt"

      open(3,file=path3,status='unknown',access='sequential',
     1  position='append')
       write(3,*)azep_inx,abs(MEs_azepin)*1000,
     1 abs(MEs_azepin_plus1)*1000,
     1 abs(MEs_azepin_minus1)*1000,
     1 station,net,stlat,stlon,evlat,evlon,
     1 Ampfl,error,mima,nsim
      close(3)

      elseif (jcomin==2) then

c --------- for Love waves
c M:Es* (MEs) - Table A1 - Ferreira & Woodhouse 2006.
      do i=1,360
      azep_inx=NINT(azep_in)

      azep_i(i)=i*pi/180.
      azep=azep_i(i)

      part1=-cmplx(0,((lambda)*
     1 ((fmom(4)*sin(pi-azep))-(fmom(5)*cos(pi-azep))*
     1 (Wdoteigen_sdpth - Weigen_sdpth))))

      part2=((lambda**2)*Weigen_sdpth)*
     1 ((0.5*(fmom(2)-fmom(3))*
     1 sin(2*(pi-azep))) -(fmom(6)*cos(2*(pi-azep))))

      MEs(i)=abs(part1+part2)
      enddo

       MEs_azepin=MEs(azep_inx)
       MEs_azepin_plus1=MEs(azep_inx+1)
       MEs_azepin_minus1=MEs(azep_inx-1)

c------ write final file output Love

c      write(path2,'(A,A,A,A,A,I1.1,A,I3.3,A)') trim(outdir),
c     1 trim(evnam),"_",trim(input_model),"_n",nmaxin,"_l",
c     1 lmaxin,"_Love.txt"

c      open(2,file=path2,status='unknown',access='sequential',
c     1  position='append')
c      do i=1,360
c       azep_i(i)=i*pi/180.
c       azep=azep_i(i)
cc ---- was (pi-azep)*180/pi, not sure why 
c       write(2,*) azep*(180/pi), abs(MEs(i))*1000
c      enddo
c      close(2)

      write(path4,'(A,I1.1,A,I3.3,A)') trim(outdir),
     1 nmaxin,"T",lmaxin,"_nodality.txt"

      open(4,file=path4,status='unknown',access='sequential',
     1  position='append')
       write(4,*)azep_inx,abs(MEs_azepin)*1000,
     1 abs(MEs_azepin_plus1)*1000,
     1 abs(MEs_azepin_minus1)*1000,
     1 station,net,stlat,stlon,evlat,evlon,
     1 Ampfl,error,mima,nsim
      close(4)

      endif

      end program
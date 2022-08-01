!-------------------
!*** mls part
!*** meshless part
!-------------------
!------------------------------
!*** general triangular solve
!*** multiple right-hand-sides
!*** chk0
!------------------------------
  subroutine mls_gentriangsolve(upper,nrhs,n,r,b,x)
    implicit real(8) (a-h,o-z)
    logical::upper
    integer::n
    real(8),dimension(n,n)::r
    real(8),dimension(n,nrhs)::b,x
    if(upper)then
!-------------------------------------------
!*** solve r.x=b where r is upper triangular
!-------------------------------------------
       do ir=1,nrhs
          do id=n,1,-1
             x(id,ir)=b(id,ir)
             do jd=id+1,n
                x(id,ir)=x(id,ir)-r(id,jd)*x(jd,ir)
             end do
             x(id,ir)=x(id,ir)/r(id,id)
          end do
       end do
    else
!---------------------------------------------
!*** solve r^t.x=b where r is upper triangular
!---------------------------------------------
       do ir=1,nrhs
          do id=1,n
             x(id,ir)=b(id,ir)
             do jd=1,id-1
                x(id,ir)=x(id,ir)-r(jd,id)*x(jd,ir)
             end do
             x(id,ir)=x(id,ir)/r(id,id)
          end do
       end do
    end if
  end subroutine mls_gentriangsolve

!--------------------------
!*** performs a naive qr
!*** decomposition
!*** chk0
!--------------------------
  subroutine mls_projqr(n,u,a,ua)
    implicit real(8) (a-h,o-z)
    integer::n
    real(8),dimension(n)::u,a,ua
    ua=u*dotprod(n,u,a)/dotprod(n,u,u)
  end subroutine mls_projqr

!--------------------------
!*** performs a naive qr
!*** decomposition
!*** chk0
!--------------------------
  subroutine mls_naiveqr(n,m,at,r)
    implicit real(8) (a-h,o-z)
    integer::n,m
    real(8),dimension(n)::temp
    real(8),dimension(n,m)::e
    real(8),dimension(n,m)::a
    real(8),dimension(m,m)::r
    real(8),dimension(n,m)::at
!-------------------
!*** sets r to zero
!-------------------
    r=0.0d00    
    do im=1,m
       e(1:n,im)=at(1:n,im)
       do jm=1,im-1
          call mls_projqr(n,e(1:n,jm),at(1:n,im),temp(1:n))
          e(1:n,im)=e(1:n,im)-temp
       end do
       call nrmali(n,e(1:n,im))
    end do
    do im=1,m
       a(1:n,im)=0.0d00
       do jm=1,im
          a(1:n,im)=a(1:n,im)+dotprod(n,e(1:n,jm),at(1:n,im))*e(1:n,jm)
       end do
    end do
    do im=1,m
       do jm=im,m
          r(im,jm)=dotprod(n,e(1:n,im),a(1:n,jm))
       end do
    end do
  end subroutine mls_naiveqr

!--------------------
!*** from w, p and b
!*** determines u2
!*** chk0 (both)
!--------------------  
  subroutine mls_determu2(m,n,w,p,u2)
    implicit real(8) (a-h,o-z)
    integer::m,n
    real(8),dimension(n)::w
    real(8),dimension(m,n)::b,u1,u2
    real(8),dimension(m,n)::p
    real(8),dimension(n,m)::at
    real(8),dimension(m,m)::r,atemp
    do im=1,m
       do in=1,n
          at(in,im)=sqrt(w(in))*p(im,in)
       end do
    end do
    do in=1,n
       do im=1,m
          b(im,in)=p(im,in)*w(in)
       end do
    end do
    call mls_naiveqr(n,m,at,r)
    call mls_gentriangsolve(.false.,n,m,r,b,u1)
    call mls_gentriangsolve(.true.,n,m,r,u1,u2)
  end subroutine mls_determu2

!---------------------------------------
!*** determines generic shape functions
!*** and derivatives
!*** chk0
!       call mls_sfder(ndi,n,m,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m),u2(1:m,1:n),ffloc,dffloc,dff2loc)
!---------------------------------------
  subroutine mls_sfder(ndi,n,m,polyn,dpolyn,dpolyn2,u2,ff,dff,dff2)
    implicit real(8) (a-h,o-z)
    integer::m,n
    real(8),dimension(m)::polyn
    real(8),dimension(ndi,m)::dpolyn
    real(8),dimension(ndi,ndi,m)::dpolyn2
    real(8),dimension(m,n)::u2
    real(8),dimension(n)::ff
    real(8),dimension(ndi,n)::dff
    real(8),dimension(ndi,ndi,n)::dff2
    do in=1,n
       ff(in)=dotprod(m,polyn(1:m),u2(1:m,in))
       do id=1,ndi
          dff(id,in)=dotprod(m,dpolyn(id,1:m),u2(1:m,in))
          do jd=1,ndi
             dff2(id,jd,in)=dotprod(m,dpolyn2(id,jd,1:m),u2(1:m,in))
          end do
       end do
    end do
  end subroutine mls_sfder

  subroutine mls_polynquad1d(d,x,xbar,polyn,dpolyn,dpolyn2)
    implicit none
    DOUBLE PRECISION v(40),d,x(1),xbar(1),polyn(*),dpolyn(1,*),dpolyn2(1,1,*)
    v(34)=1d0/d
    v(33)=x(1)-xbar(1)
    v(32)=1d0/d**2
    v(35)=2d0*v(32)
    v(31)=1d0/d**3
    v(23)=(v(33)*v(33))
    polyn(1)=1d0
    polyn(2)=v(33)*v(34)
    polyn(3)=v(23)*v(32)
    polyn(4)=v(31)*v(33)**3
    dpolyn(1,1)=0d0
    dpolyn(1,2)=v(34)
    dpolyn(1,3)=v(33)*v(35)
    dpolyn(1,4)=3d0*v(23)*v(31)
    dpolyn2(1,1,1)=0d0
    dpolyn2(1,1,2)=0d0
    dpolyn2(1,1,3)=v(35)
    dpolyn2(1,1,4)=6d0*v(31)*v(33)
  end subroutine mls_polynquad1d

  subroutine mls_polynquad2d(d,x,xbar,polyn,dpolyn,dpolyn2)
    implicit none
    DOUBLE PRECISION v(124),d,x(2),xbar(2),polyn(*),dpolyn(2,*),dpolyn2(2,2,*)
    v(115)=1d0/d
    v(113)=x(2)-xbar(2)
    v(112)=x(1)-xbar(1)
    v(111)=1d0/d**2
    v(116)=v(111)*v(113)
    v(110)=1d0/d**3
    v(119)=3d0*v(110)
    v(114)=v(110)*v(113)
    v(91)=2d0*v(112)
    v(87)=(v(112)*v(112))
    v(117)=v(110)*v(87)
    v(97)=2d0*v(113)
    v(89)=(v(113)*v(113))
    v(118)=v(110)*v(89)
    v(95)=v(114)*v(91)
    v(102)=2d0*v(111)
    v(104)=2d0*v(114)
    v(105)=v(110)*v(91)
    polyn(1)=1d0
    polyn(2)=v(112)*v(115)
    polyn(3)=v(113)*v(115)
    polyn(4)=v(111)*v(87)
    polyn(5)=v(111)*v(89)
    polyn(6)=v(112)*v(116)
    polyn(7)=v(110)*v(112)**3
    polyn(8)=v(110)*v(113)**3
    polyn(9)=v(114)*v(87)
    polyn(10)=v(112)*v(118)
    dpolyn(1,1)=0d0
    dpolyn(1,2)=v(115)
    dpolyn(1,3)=0d0
    dpolyn(1,4)=v(111)*v(91)
    dpolyn(1,5)=0d0
    dpolyn(1,6)=v(116)
    dpolyn(1,7)=3d0*v(117)
    dpolyn(1,8)=0d0
    dpolyn(1,9)=v(95)
    dpolyn(1,10)=v(118)
    dpolyn(2,1)=0d0
    dpolyn(2,2)=0d0
    dpolyn(2,3)=v(115)
    dpolyn(2,4)=0d0
    dpolyn(2,5)=v(111)*v(97)
    dpolyn(2,6)=v(111)*v(112)
    dpolyn(2,7)=0d0
    dpolyn(2,8)=3d0*v(118)
    dpolyn(2,9)=v(117)
    dpolyn(2,10)=v(95)
    dpolyn2(1,1,1)=0d0
    dpolyn2(1,1,2)=0d0
    dpolyn2(1,1,3)=0d0
    dpolyn2(1,1,4)=v(102)
    dpolyn2(1,1,5)=0d0
    dpolyn2(1,1,6)=0d0
    dpolyn2(1,1,7)=v(119)*v(91)
    dpolyn2(1,1,8)=0d0
    dpolyn2(1,1,9)=v(104)
    dpolyn2(1,1,10)=0d0
    dpolyn2(1,2,1)=0d0
    dpolyn2(1,2,2)=0d0
    dpolyn2(1,2,3)=0d0
    dpolyn2(1,2,4)=0d0
    dpolyn2(1,2,5)=0d0
    dpolyn2(1,2,6)=v(111)
    dpolyn2(1,2,7)=0d0
    dpolyn2(1,2,8)=0d0
    dpolyn2(1,2,9)=v(105)
    dpolyn2(1,2,10)=v(104)
    dpolyn2(2,1,1)=0d0
    dpolyn2(2,1,2)=0d0
    dpolyn2(2,1,3)=0d0
    dpolyn2(2,1,4)=0d0
    dpolyn2(2,1,5)=0d0
    dpolyn2(2,1,6)=v(111)
    dpolyn2(2,1,7)=0d0
    dpolyn2(2,1,8)=0d0
    dpolyn2(2,1,9)=v(105)
    dpolyn2(2,1,10)=v(104)
    dpolyn2(2,2,1)=0d0
    dpolyn2(2,2,2)=0d0
    dpolyn2(2,2,3)=0d0
    dpolyn2(2,2,4)=0d0
    dpolyn2(2,2,5)=v(102)
    dpolyn2(2,2,6)=0d0
    dpolyn2(2,2,7)=0d0
    dpolyn2(2,2,8)=v(119)*v(97)
    dpolyn2(2,2,9)=0d0
    dpolyn2(2,2,10)=v(105)
  end subroutine mls_polynquad2d

  subroutine mls_polynquad3d(d,x,xbar,polyn,dpolyn,dpolyn2)
    implicit none
    DOUBLE PRECISION v(340),d,x(3),xbar(3),polyn(*),dpolyn(3,*),dpolyn2(3,3,*)
    v(334)=1d0/d
    v(333)=x(3)-xbar(3)
    v(332)=x(2)-xbar(2)
    v(331)=x(1)-xbar(1)
    v(330)=1d0/d**2
    v(329)=1d0/d**3
    v(335)=3d0*v(329)
    v(295)=2d0*v(331)
    v(286)=(v(331)*v(331))
    v(316)=v(329)*v(332)
    v(305)=2d0*v(332)
    v(289)=(v(332)*v(332))
    v(313)=2d0*v(333)
    v(309)=v(329)*v(333)
    v(300)=v(309)*v(332)
    v(292)=(v(333)*v(333))
    v(297)=v(330)*v(332)
    v(298)=v(330)*v(333)
    v(301)=v(295)*v(316)
    v(302)=v(295)*v(309)
    v(303)=v(289)*v(329)
    v(304)=v(292)*v(329)
    v(307)=v(330)*v(331)
    v(311)=v(286)*v(329)
    v(312)=v(305)*v(309)
    v(318)=2d0*v(330)
    v(320)=v(305)*v(329)
    v(321)=2d0*v(309)
    v(322)=v(295)*v(329)
    v(324)=v(329)*v(331)
    polyn(1)=1d0
    polyn(2)=v(331)*v(334)
    polyn(3)=v(332)*v(334)
    polyn(4)=v(333)*v(334)
    polyn(5)=v(286)*v(330)
    polyn(6)=v(289)*v(330)
    polyn(7)=v(292)*v(330)
    polyn(8)=v(307)*v(332)
    polyn(9)=v(307)*v(333)
    polyn(10)=v(298)*v(332)
    polyn(11)=v(329)*v(331)**3
    polyn(12)=v(329)*v(332)**3
    polyn(13)=v(329)*v(333)**3
    polyn(14)=v(300)*v(331)
    polyn(15)=v(311)*v(332)
    polyn(16)=v(311)*v(333)
    polyn(17)=v(289)*v(324)
    polyn(18)=v(303)*v(333)
    polyn(19)=v(292)*v(324)
    polyn(20)=v(304)*v(332)
    dpolyn(1,1)=0d0
    dpolyn(1,2)=v(334)
    dpolyn(1,3)=0d0
    dpolyn(1,4)=0d0
    dpolyn(1,5)=v(295)*v(330)
    dpolyn(1,6)=0d0
    dpolyn(1,7)=0d0
    dpolyn(1,8)=v(297)
    dpolyn(1,9)=v(298)
    dpolyn(1,10)=0d0
    dpolyn(1,11)=3d0*v(311)
    dpolyn(1,12)=0d0
    dpolyn(1,13)=0d0
    dpolyn(1,14)=v(300)
    dpolyn(1,15)=v(301)
    dpolyn(1,16)=v(302)
    dpolyn(1,17)=v(303)
    dpolyn(1,18)=0d0
    dpolyn(1,19)=v(304)
    dpolyn(1,20)=0d0
    dpolyn(2,1)=0d0
    dpolyn(2,2)=0d0
    dpolyn(2,3)=v(334)
    dpolyn(2,4)=0d0
    dpolyn(2,5)=0d0
    dpolyn(2,6)=v(305)*v(330)
    dpolyn(2,7)=0d0
    dpolyn(2,8)=v(307)
    dpolyn(2,9)=0d0
    dpolyn(2,10)=v(298)
    dpolyn(2,11)=0d0
    dpolyn(2,12)=3d0*v(303)
    dpolyn(2,13)=0d0
    dpolyn(2,14)=v(309)*v(331)
    dpolyn(2,15)=v(311)
    dpolyn(2,16)=0d0
    dpolyn(2,17)=v(301)
    dpolyn(2,18)=v(312)
    dpolyn(2,19)=0d0
    dpolyn(2,20)=v(304)
    dpolyn(3,1)=0d0
    dpolyn(3,2)=0d0
    dpolyn(3,3)=0d0
    dpolyn(3,4)=v(334)
    dpolyn(3,5)=0d0
    dpolyn(3,6)=0d0
    dpolyn(3,7)=v(313)*v(330)
    dpolyn(3,8)=0d0
    dpolyn(3,9)=v(307)
    dpolyn(3,10)=v(297)
    dpolyn(3,11)=0d0
    dpolyn(3,12)=0d0
    dpolyn(3,13)=3d0*v(304)
    dpolyn(3,14)=v(316)*v(331)
    dpolyn(3,15)=0d0
    dpolyn(3,16)=v(311)
    dpolyn(3,17)=0d0
    dpolyn(3,18)=v(303)
    dpolyn(3,19)=v(302)
    dpolyn(3,20)=v(312)
    dpolyn2(1,1,1)=0d0
    dpolyn2(1,1,2)=0d0
    dpolyn2(1,1,3)=0d0
    dpolyn2(1,1,4)=0d0
    dpolyn2(1,1,5)=v(318)
    dpolyn2(1,1,6)=0d0
    dpolyn2(1,1,7)=0d0
    dpolyn2(1,1,8)=0d0
    dpolyn2(1,1,9)=0d0
    dpolyn2(1,1,10)=0d0
    dpolyn2(1,1,11)=v(295)*v(335)
    dpolyn2(1,1,12)=0d0
    dpolyn2(1,1,13)=0d0
    dpolyn2(1,1,14)=0d0
    dpolyn2(1,1,15)=v(320)
    dpolyn2(1,1,16)=v(321)
    dpolyn2(1,1,17)=0d0
    dpolyn2(1,1,18)=0d0
    dpolyn2(1,1,19)=0d0
    dpolyn2(1,1,20)=0d0
    dpolyn2(1,2,1)=0d0
    dpolyn2(1,2,2)=0d0
    dpolyn2(1,2,3)=0d0
    dpolyn2(1,2,4)=0d0
    dpolyn2(1,2,5)=0d0
    dpolyn2(1,2,6)=0d0
    dpolyn2(1,2,7)=0d0
    dpolyn2(1,2,8)=v(330)
    dpolyn2(1,2,9)=0d0
    dpolyn2(1,2,10)=0d0
    dpolyn2(1,2,11)=0d0
    dpolyn2(1,2,12)=0d0
    dpolyn2(1,2,13)=0d0
    dpolyn2(1,2,14)=v(309)
    dpolyn2(1,2,15)=v(322)
    dpolyn2(1,2,16)=0d0
    dpolyn2(1,2,17)=v(320)
    dpolyn2(1,2,18)=0d0
    dpolyn2(1,2,19)=0d0
    dpolyn2(1,2,20)=0d0
    dpolyn2(1,3,1)=0d0
    dpolyn2(1,3,2)=0d0
    dpolyn2(1,3,3)=0d0
    dpolyn2(1,3,4)=0d0
    dpolyn2(1,3,5)=0d0
    dpolyn2(1,3,6)=0d0
    dpolyn2(1,3,7)=0d0
    dpolyn2(1,3,8)=0d0
    dpolyn2(1,3,9)=v(330)
    dpolyn2(1,3,10)=0d0
    dpolyn2(1,3,11)=0d0
    dpolyn2(1,3,12)=0d0
    dpolyn2(1,3,13)=0d0
    dpolyn2(1,3,14)=v(316)
    dpolyn2(1,3,15)=0d0
    dpolyn2(1,3,16)=v(322)
    dpolyn2(1,3,17)=0d0
    dpolyn2(1,3,18)=0d0
    dpolyn2(1,3,19)=v(321)
    dpolyn2(1,3,20)=0d0
    dpolyn2(2,1,1)=0d0
    dpolyn2(2,1,2)=0d0
    dpolyn2(2,1,3)=0d0
    dpolyn2(2,1,4)=0d0
    dpolyn2(2,1,5)=0d0
    dpolyn2(2,1,6)=0d0
    dpolyn2(2,1,7)=0d0
    dpolyn2(2,1,8)=v(330)
    dpolyn2(2,1,9)=0d0
    dpolyn2(2,1,10)=0d0
    dpolyn2(2,1,11)=0d0
    dpolyn2(2,1,12)=0d0
    dpolyn2(2,1,13)=0d0
    dpolyn2(2,1,14)=v(309)
    dpolyn2(2,1,15)=v(322)
    dpolyn2(2,1,16)=0d0
    dpolyn2(2,1,17)=v(320)
    dpolyn2(2,1,18)=0d0
    dpolyn2(2,1,19)=0d0
    dpolyn2(2,1,20)=0d0
    dpolyn2(2,2,1)=0d0
    dpolyn2(2,2,2)=0d0
    dpolyn2(2,2,3)=0d0
    dpolyn2(2,2,4)=0d0
    dpolyn2(2,2,5)=0d0
    dpolyn2(2,2,6)=v(318)
    dpolyn2(2,2,7)=0d0
    dpolyn2(2,2,8)=0d0
    dpolyn2(2,2,9)=0d0
    dpolyn2(2,2,10)=0d0
    dpolyn2(2,2,11)=0d0
    dpolyn2(2,2,12)=3d0*v(320)
    dpolyn2(2,2,13)=0d0
    dpolyn2(2,2,14)=0d0
    dpolyn2(2,2,15)=0d0
    dpolyn2(2,2,16)=0d0
    dpolyn2(2,2,17)=v(322)
    dpolyn2(2,2,18)=v(321)
    dpolyn2(2,2,19)=0d0
    dpolyn2(2,2,20)=0d0
    dpolyn2(2,3,1)=0d0
    dpolyn2(2,3,2)=0d0
    dpolyn2(2,3,3)=0d0
    dpolyn2(2,3,4)=0d0
    dpolyn2(2,3,5)=0d0
    dpolyn2(2,3,6)=0d0
    dpolyn2(2,3,7)=0d0
    dpolyn2(2,3,8)=0d0
    dpolyn2(2,3,9)=0d0
    dpolyn2(2,3,10)=v(330)
    dpolyn2(2,3,11)=0d0
    dpolyn2(2,3,12)=0d0
    dpolyn2(2,3,13)=0d0
    dpolyn2(2,3,14)=v(324)
    dpolyn2(2,3,15)=0d0
    dpolyn2(2,3,16)=0d0
    dpolyn2(2,3,17)=0d0
    dpolyn2(2,3,18)=v(320)
    dpolyn2(2,3,19)=0d0
    dpolyn2(2,3,20)=v(321)
    dpolyn2(3,1,1)=0d0
    dpolyn2(3,1,2)=0d0
    dpolyn2(3,1,3)=0d0
    dpolyn2(3,1,4)=0d0
    dpolyn2(3,1,5)=0d0
    dpolyn2(3,1,6)=0d0
    dpolyn2(3,1,7)=0d0
    dpolyn2(3,1,8)=0d0
    dpolyn2(3,1,9)=v(330)
    dpolyn2(3,1,10)=0d0
    dpolyn2(3,1,11)=0d0
    dpolyn2(3,1,12)=0d0
    dpolyn2(3,1,13)=0d0
    dpolyn2(3,1,14)=v(316)
    dpolyn2(3,1,15)=0d0
    dpolyn2(3,1,16)=v(322)
    dpolyn2(3,1,17)=0d0
    dpolyn2(3,1,18)=0d0
    dpolyn2(3,1,19)=v(321)
    dpolyn2(3,1,20)=0d0
    dpolyn2(3,2,1)=0d0
    dpolyn2(3,2,2)=0d0
    dpolyn2(3,2,3)=0d0
    dpolyn2(3,2,4)=0d0
    dpolyn2(3,2,5)=0d0
    dpolyn2(3,2,6)=0d0
    dpolyn2(3,2,7)=0d0
    dpolyn2(3,2,8)=0d0
    dpolyn2(3,2,9)=0d0
    dpolyn2(3,2,10)=v(330)
    dpolyn2(3,2,11)=0d0
    dpolyn2(3,2,12)=0d0
    dpolyn2(3,2,13)=0d0
    dpolyn2(3,2,14)=v(324)
    dpolyn2(3,2,15)=0d0
    dpolyn2(3,2,16)=0d0
    dpolyn2(3,2,17)=0d0
    dpolyn2(3,2,18)=v(320)
    dpolyn2(3,2,19)=0d0
    dpolyn2(3,2,20)=v(321)
    dpolyn2(3,3,1)=0d0
    dpolyn2(3,3,2)=0d0
    dpolyn2(3,3,3)=0d0
    dpolyn2(3,3,4)=0d0
    dpolyn2(3,3,5)=0d0
    dpolyn2(3,3,6)=0d0
    dpolyn2(3,3,7)=v(318)
    dpolyn2(3,3,8)=0d0
    dpolyn2(3,3,9)=0d0
    dpolyn2(3,3,10)=0d0
    dpolyn2(3,3,11)=0d0
    dpolyn2(3,3,12)=0d0
    dpolyn2(3,3,13)=v(313)*v(335)
    dpolyn2(3,3,14)=0d0
    dpolyn2(3,3,15)=0d0
    dpolyn2(3,3,16)=0d0
    dpolyn2(3,3,17)=0d0
    dpolyn2(3,3,18)=0d0
    dpolyn2(3,3,19)=v(322)
    dpolyn2(3,3,20)=v(320)
  end subroutine mls_polynquad3d

!---------------------------------
!*** mls_polynquadbase determines
!*** poly. base + derivatives
!*** 1d to 3d
!*** chk0
!---------------------------------  
  subroutine mls_polyncubbase(ndi,dtemp,x,xbar,polyn,dpolyn,dpolyn2,m)
    implicit real(8) (a-h,o-z)
!*** check this parameter for updates (mpol)
    real(8)::d,dtemp
    integer,parameter::mpol=20
    real(8),dimension(ndi)::x,xbar
    REAL(8),DIMENSION(*)::polyn
    REAL(8),DIMENSION(ndi,*)::dpolyn
    REAL(8),DIMENSION(ndi,ndi,*)::dpolyn2
    d=dtemp
    select case(ndi)
    case(1)
       m=4
       call mls_polynquad1d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case(2)
       m=10
       call mls_polynquad2d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case(3)      
       m=20      
       call mls_polynquad3d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case default
       stop "case not found in mls_polyncubbase"
    end select
  end subroutine mls_polyncubbase

!---------------------------------
!*** mls_polynquadbase determines
!*** poly. base + derivatives
!*** 1d to 3d
!*** chk0
!---------------------------------  
  subroutine mls_polynqbase(ndi,dtemp,x,xbar,polyn,dpolyn,dpolyn2,m)
    implicit real(8) (a-h,o-z)
!*** check this parameter for updates (mpol)
    real(8)::d,dtemp
    integer,parameter::mpol=20
    real(8),dimension(ndi)::x,xbar
    REAL(8),DIMENSION(*)::polyn
    REAL(8),DIMENSION(ndi,*)::dpolyn
    REAL(8),DIMENSION(ndi,ndi,*)::dpolyn2
    d=dtemp
    select case(ndi)
    case(1)
       m=3
       call mls_polynquad1d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case(2)
       m=6
       call mls_polynquad2d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case(3)      
       m=10       
       call mls_polynquad3d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case default
       stop "case not found in mls_polynqbase"
    end select
  end subroutine mls_polynqbase

!---------------------------------
!*** mls_polynquadbase determines
!*** poly. base + derivatives
!*** 1d to 3d
!*** chk0
!---------------------------------  
  subroutine mls_polynlinbase(ndi,dtemp,x,xbar,polyn,dpolyn,dpolyn2,m)
    implicit real(8) (a-h,o-z)
!*** check this parameter for updates (mpol)
    real(8)::d,dtemp
    integer,parameter::mpol=20
    real(8),dimension(ndi)::x,xbar
    REAL(8),DIMENSION(*)::polyn
    REAL(8),DIMENSION(ndi,*)::dpolyn
    REAL(8),DIMENSION(ndi,ndi,*)::dpolyn2
    d=dtemp
    select case(ndi)
    case(1)
       m=2
       call mls_polynquad1d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case(2)
       m=3
       call mls_polynquad2d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case(3)      
       m=4       
       call mls_polynquad3d(d,x,xbar,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m))
    case default
       stop "case not found in mls_polynlinbase"
    end select
  end subroutine mls_polynlinbase

!----------------------------
!*** weight function for mls
!*** ndi
!*** chk0
!----------------------------
  subroutine mls_weight(ndi,d,tol,xi,x,w)
    implicit real(8) (a-h,o-z)
    real(8),dimension(ndi)::xi,x
    real(8),parameter::pi=4.0d00*atan(1.0d00)
    real(8),parameter::bottom=1.0d-16
    s=rnorm2(ndi,xi-x)
!!    wmax=1.0d00/(tol**2.0d00+tol**4.0d00)
!!    wbot=bottom*wmax
!!    w=max((1.0d00/((s*s)/(d*d)+tol*tol))-(1.0d00/(1.0d00+tol*tol)),wbot)
    w=(1.0d00/((s*s)/(d*d)+tol*tol))-(1.0d00/(1.0d00+tol*tol))
!!    w=max(exp(-s/(1.0d-6*d)),1.0d-15)
!!    w=sqrt(1.0d00/(tol*pi))*exp(-((s/d)**2)/tol)
!!    w=1.0d00/((s/d)+tol)-1.0d00/(1.0d00+tol)
!!    w=max(w,0.0d00)
    w=max(w,0.0d00)
!    w=1.0d00
!    w=1.0d00
!!    w=(1.0d00/((s*s)/(d*d)+tol*tol))-(1.0d00/(1.0d00+tol*tol))
!!    w=max(w,1.0d-20)
!    w=1.0d00
  end subroutine mls_weight

!--------------------------------------------------
!*** determines u2 evaluated in a given coordinate
!    call u2atacoordinate(d,tol,ndi,n,xn,x,u2)
!--------------------------------------------------
  subroutine mls_u2atacoordinate(imeshless,d,tol,ndi,n,xbar,xn,x,u2)
    implicit real(8) (a-h,o-z)
    real(8)::d,tol
    integer,parameter::mpol=20
    real(8),dimension(n)::w
    real(8),dimension(mpol,n)::p
    real(8),dimension(mpol,n)::u2
    real(8),dimension(ndi,n)::xn
    real(8),dimension(ndi)::x,xbar
    real(8),dimension(3,mpol)::trash
    real(8),dimension(3,3,mpol)::trash2
!-------------------------------
!*** now we must define matrix p
!*** derivatives are discarded
!-------------------------------
    select case(imeshless)
    case(-1)
       do in=1,n
          call mls_polynlinbase(ndi,d,xn(1:ndi,in),xbar(1:ndi),p(1:mpol,in),trash(1:ndi,1:mpol),trash2(1:ndi,1:ndi,1:mpol),m)
       end do
    case(-2)
       do in=1,n
          call mls_polynqbase(ndi,d,xn(1:ndi,in),xbar(1:ndi),p(1:mpol,in),trash(1:ndi,1:mpol),trash2(1:ndi,1:ndi,1:mpol),m)
       end do
    case(-3)
       do in=1,n
          call mls_polyncubbase(ndi,d,xn(1:ndi,in),xbar(1:ndi),p(1:mpol,in),trash(1:ndi,1:mpol),trash2(1:ndi,1:ndi,1:mpol),m)
       end do
    case default
       stop "wrong request to u2atacoordinate"
    end select
!------------------------
!*** now we must define w
!------------------------
    do in=1,n
       call mls_weight(ndi,d,tol,xn(1:ndi,in),x(1:ndi),w(in))
    end do
!------------------------------
!*** finalize it with specifics
!------------------------------
    call mls_determu2(m,n,w(1:n),p(1:m,1:n),u2(1:m,1:n))
  end subroutine mls_u2atacoordinate

!-------------------------
!*** gets influence nodes
!*** for a given point
!*** chk0
!-------------------------
  subroutine mls_getsnodes(rmax,x,ndi,ntot,xtot,n,listn)
    implicit real(8)(a-h,o-z)
    real(8)::rmax
    real(8),dimension(ndi,ntot)::xtot
    real(8),dimension(ndi)::x
    integer::ndi
    integer,dimension(:),allocatable::listn
    n=0
    do itot=1,ntot
       if(rnorm2(ndi,xtot(1:ndi,itot)-x(1:ndi)).le.rmax)then
          n=n+1
       end if
    end do
    allocate(listn(n))
!--------------------------------------
!*** now inserts the nodes in the list
!--------------------------------------
    n=0
    do itot=1,ntot
       if(rnorm2(ndi,xtot(1:ndi,itot)-x(1:ndi)).le.rmax)then
          n=n+1
          listn(n)=itot
       end if
    end do
  end subroutine mls_getsnodes

!---------------------------
!*** gets precisely n nodes
!---------------------------
  subroutine mls_getsnnodes(rmax,x,ndi,ntot,xtot,n,listn)
    implicit real(8)(a-h,o-z)
    real(8)::rmax
    real(8),dimension(ndi,ntot)::xtot
    real(8),dimension(ndi)::x
    real(8),dimension(ntot)::ranking
    integer,dimension(ntot)::per
    integer::ndi
    integer,dimension(:),allocatable::listn
    if(n.le.0)then
       call mls_getsnodes(rmax,x,ndi,ntot,xtot,n,listn)
    else
       do ino=1,ntot
          ranking(ino)=rnorm2(ndi,x(1:ndi)-xtot(1:ndi,ino))
       end do
!--------------------
!*** sort by ranking
!--------------------
       call sort(ntot,ranking,per)
       allocate(listn(n))
       do i=1,n
          listn(i)=per(i)
       end do
    end if
  end subroutine mls_getsnnodes

  subroutine mls_cauchygreen(ndi,n,xtot,xdef,dff,f,cauchygreen)
    implicit real(8) (a-h,o-z)
    real(8),dimension(ndi,n)::xtot,xdef
    real(8),dimension(ndi,ndi)::c,f,fo,co
    real(8),dimension(ndi,n)::dff
    real(8),dimension(ndi*(ndi+1)/2)::cauchygreen
!-------------------------
!*** deformation gradient
!*** before correction
!-------------------------
    f=0.0d00
    do in=1,n
       do jd=1,ndi
          do id=1,ndi
             f(id,jd)=f(id,jd)+dff(jd,in)*xdef(id,in)
          end do
       end do
    end do
!------------------------------
!*** right cauchy-green tensor
!------------------------------
    call matmat(ndi,ndi,ndi,f,f,c,2)
    do id=1,ndi*(ndi+1)/2
       call apomat(i1,i2,id,ndi)
       cauchygreen(id)=c(i1,i2)
    end do
  end subroutine mls_cauchygreen

  subroutine mls_strain(ndi,n,xtot,xdef,dff,f,strain)
    implicit real(8) (a-h,o-z)
    real(8),dimension(ndi,n)::xtot,xdef
    real(8),dimension(ndi,ndi)::c,f,fo,co
    real(8),dimension(ndi,n)::dff
    real(8),dimension(ndi*(ndi+1)/2)::strain,straino
!-------------------------
!*** deformation gradient
!*** before correction
!-------------------------
    f=0.0d00
    do in=1,n
       do jd=1,ndi
          do id=1,ndi
             f(id,jd)=f(id,jd)+dff(jd,in)*xdef(id,in)
          end do
       end do
    end do
!------------------------------
!*** right cauchy-green tensor
!------------------------------
    call matmat(ndi,ndi,ndi,f,f,c,2)
    do id=1,ndi*(ndi+1)/2
       call apomat(i1,i2,id,ndi)
       strain(id)=0.5d00*(c(i1,i2)-deltak(i1,i2))
    end do
    do id=ndi+1,ndi*(ndi+1)/2
       strain(id)=2.0d00*strain(id)
    end do
  end subroutine mls_strain

!----------------------------------------------------
!*** determines the deformation gradient
!*** and green-lagrange strain in engineering form
!*** mls stuff
!*** chk0,1
!----------------------------------------------------
  subroutine mls_strain2(ndi,n,grad,grad2,e,e2)
    implicit real(8) (a-h,o-z)
!----------------------------------------------------
    real(8),dimension(ndi,ndi)::grad
    real(8),dimension(ndi,ndi,ndi)::grad2
!----------------------------------------------------
    real(8),dimension(ndi*(ndi+1)/2)::e
    real(8),dimension(ndi*(ndi+1)/2,ndi)::e2
!----------------------------------------------------
    nvoigt=ndi*(ndi+1)/2
    do ij=1,nvoigt
       call apomat(i,j,ij,ndi)
!e
       e(ij)=0.0d00
       do k=1,ndi
          e(ij)=e(ij)+0.5d00*grad(k,i)*grad(k,j)
       end do
       e(ij)=e(ij)-0.5d00*deltak(i,j)
!e2
       do l=1,ndi
          e2(ij,l)=0.0d00
          do k=1,ndi
             e2(ij,l)=e2(ij,l)+0.5d00*grad2(k,i,l)*grad(k,j)+0.5d00*grad(k,i)*grad2(k,j,l)
          end do
       end do
    end do
  end subroutine mls_strain2

  subroutine mls_2dforce(nk,f,s,nucleus)
    implicit none
    double precision v(27),nk(2),f(2,2),s(3),nucleus(2)
    v(22)=nk(2)*s(2)
    v(21)=nk(1)*s(1)
    nucleus(1)=(f(1,2)*nk(1)+f(1,1)*nk(2))*s(3)+f(1,1)*v(21)+f(1,2)*v(22)
    nucleus(2)=(f(2,2)*nk(1)+f(2,1)*nk(2))*s(3)+f(2,1)*v(21)+f(2,2)*v(22)
  end subroutine mls_2dforce

  subroutine mls_2d(nk,nl,f,s,ds,kernel)
    implicit none
    double precision v(67),nk(2),nl(2),f(2,2),s(3),ds(3,3),kernel(2,2)
    v(53)=nk(2)*(nl(2)*s(2)+nl(1)*s(3))+nk(1)*(nl(1)*s(1)+nl(2)*s(3))
    v(52)=f(2,2)*nk(1)+f(2,1)*nk(2)
    v(51)=f(1,2)*nk(1)+f(1,1)*nk(2)
    v(50)=f(2,2)*nk(2)
    v(49)=f(1,2)*nk(2)
    v(48)=f(2,1)*nk(1)
    v(47)=f(1,1)*nk(1)
    v(46)=f(2,2)*nl(1)+f(2,1)*nl(2)
    v(45)=f(1,2)*nl(1)+f(1,1)*nl(2)
    v(44)=f(2,2)*nl(2)
    v(60)=ds(2,2)*v(44)
    v(59)=ds(1,2)*v(44)+ds(1,3)*v(46)
    v(43)=f(1,2)*nl(2)
    v(54)=ds(1,2)*v(43)+ds(1,3)*v(45)
    v(42)=f(2,1)*nl(1)
    v(62)=ds(3,1)*v(42)+ds(3,2)*v(44)+ds(3,3)*v(46)
    v(61)=ds(2,1)*v(42)+ds(2,3)*v(46)
    v(58)=ds(1,1)*v(42)
    v(41)=f(1,1)*nl(1)
    v(56)=ds(3,1)*v(41)+ds(3,2)*v(43)+ds(3,3)*v(45)
    v(55)=ds(2,1)*v(41)+ds(2,3)*v(45)
    v(38)=v(47)*v(58)
    v(39)=v(49)*v(60)
    v(57)=v(38)+v(39)
    kernel(1,1)=v(53)+v(47)*(ds(1,1)*v(41)+v(54))+v(49)*(ds(2,2)*v(43)+v(55))+v(51)*v(56)
    kernel(1,2)=v(57)+v(47)*v(59)+v(49)*v(61)+v(51)*v(62)
    kernel(2,1)=v(48)*v(54)+v(50)*v(55)+v(52)*v(56)+v(57)
    kernel(2,2)=v(53)+v(48)*(v(58)+v(59))+v(50)*(v(60)+v(61))+v(52)*v(62)
  end subroutine mls_2d

  subroutine mls_3dforce(nk,f,s,nucleus)
    implicit none
    double precision v(51),nk(3),f(3,3),s(6),nucleus(3)
    v(46)=nk(3)*s(3)
    v(45)=nk(2)*s(2)
    v(44)=nk(1)*s(1)
    nucleus(1)=(f(1,2)*nk(1)+f(1,1)*nk(2))*s(4)+(f(1,3)*nk(1)+f(1,1)*nk(3))*s(5)+(f(1,3)*nk(2)+f(1,2)*nk(3))*s(6)+f(1,1)*v&
         &(44)+f(1,2)*v(45)+f(1,3)*v(46)
    nucleus(2)=(f(2,2)*nk(1)+f(2,1)*nk(2))*s(4)+(f(2,3)*nk(1)+f(2,1)*nk(3))*s(5)+(f(2,3)*nk(2)+f(2,2)*nk(3))*s(6)+f(2,1)*v&
         &(44)+f(2,2)*v(45)+f(2,3)*v(46)
    nucleus(3)=(f(3,2)*nk(1)+f(3,1)*nk(2))*s(4)+(f(3,3)*nk(1)+f(3,1)*nk(3))*s(5)+(f(3,3)*nk(2)+f(3,2)*nk(3))*s(6)+f(3,1)*v&
         &(44)+f(3,2)*v(45)+f(3,3)*v(46)
  end subroutine mls_3dforce

  subroutine mls_3d(nk,nl,f,s,ds,kernel)
    implicit none
    double precision v(173),nk(3),nl(3),f(3,3),s(6),ds(6,6),kernel(3,3)
    v(144)=nk(1)*(nl(1)*s(1)+nl(2)*s(4)+nl(3)*s(5))+nk(3)*(nl(3)*s(3)+nl(1)*s(5)+nl(2)*s(6))+nk(2)*(nl(2)*s(2)+nl(1)*s(4)&
         &+nl(3)*s(6))
    v(143)=f(3,3)*nk(2)+f(3,2)*nk(3)
    v(142)=f(2,3)*nk(2)+f(2,2)*nk(3)
    v(141)=f(1,3)*nk(2)+f(1,2)*nk(3)
    v(140)=f(3,3)*nk(1)+f(3,1)*nk(3)
    v(139)=f(2,3)*nk(1)+f(2,1)*nk(3)
    v(138)=f(1,3)*nk(1)+f(1,1)*nk(3)
    v(137)=f(3,2)*nk(1)+f(3,1)*nk(2)
    v(136)=f(2,2)*nk(1)+f(2,1)*nk(2)
    v(135)=f(1,2)*nk(1)+f(1,1)*nk(2)
    v(134)=f(3,3)*nk(3)
    v(133)=f(2,3)*nk(3)
    v(132)=f(1,3)*nk(3)
    v(131)=f(3,2)*nk(2)
    v(130)=f(2,2)*nk(2)
    v(129)=f(1,2)*nk(2)
    v(128)=f(3,1)*nk(1)
    v(127)=f(2,1)*nk(1)
    v(126)=f(1,1)*nk(1)
    v(125)=f(3,3)*nl(2)+f(3,2)*nl(3)
    v(124)=f(2,3)*nl(2)+f(2,2)*nl(3)
    v(123)=f(1,3)*nl(2)+f(1,2)*nl(3)
    v(122)=f(3,3)*nl(1)+f(3,1)*nl(3)
    v(121)=f(2,3)*nl(1)+f(2,1)*nl(3)
    v(120)=f(1,3)*nl(1)+f(1,1)*nl(3)
    v(119)=f(3,2)*nl(1)+f(3,1)*nl(2)
    v(118)=f(2,2)*nl(1)+f(2,1)*nl(2)
    v(117)=f(1,2)*nl(1)+f(1,1)*nl(2)
    v(116)=f(3,3)*nl(3)
    v(115)=f(2,3)*nl(3)
    v(168)=ds(3,3)*v(115)
    v(114)=f(1,3)*nl(3)
    v(147)=ds(3,3)*v(114)
    v(113)=f(3,2)*nl(2)
    v(160)=ds(1,2)*v(113)+ds(1,3)*v(116)+ds(1,4)*v(119)+ds(1,5)*v(122)+ds(1,6)*v(125)
    v(112)=f(2,2)*nl(2)
    v(167)=ds(2,2)*v(112)
    v(154)=ds(1,2)*v(112)+ds(1,3)*v(115)+ds(1,4)*v(118)+ds(1,5)*v(121)+ds(1,6)*v(124)
    v(111)=f(1,2)*nl(2)
    v(148)=ds(1,2)*v(111)+ds(1,3)*v(114)+ds(1,4)*v(117)+ds(1,5)*v(120)+ds(1,6)*v(123)
    v(146)=ds(2,2)*v(111)
    v(110)=f(3,1)*nl(1)
    v(165)=ds(6,1)*v(110)+ds(6,2)*v(113)+ds(6,3)*v(116)+ds(6,4)*v(119)+ds(6,5)*v(122)+ds(6,6)*v(125)
    v(164)=ds(5,1)*v(110)+ds(5,2)*v(113)+ds(5,3)*v(116)+ds(5,4)*v(119)+ds(5,5)*v(122)+ds(5,6)*v(125)
    v(163)=ds(4,1)*v(110)+ds(4,2)*v(113)+ds(4,3)*v(116)+ds(4,4)*v(119)+ds(4,5)*v(122)+ds(4,6)*v(125)
    v(162)=ds(3,1)*v(110)+ds(3,2)*v(113)+ds(3,4)*v(119)+ds(3,5)*v(122)+ds(3,6)*v(125)
    v(161)=ds(2,1)*v(110)+ds(2,3)*v(116)+ds(2,4)*v(119)+ds(2,5)*v(122)+ds(2,6)*v(125)
    v(109)=f(2,1)*nl(1)
    v(166)=ds(1,1)*v(109)
    v(159)=ds(6,1)*v(109)+ds(6,2)*v(112)+ds(6,3)*v(115)+ds(6,4)*v(118)+ds(6,5)*v(121)+ds(6,6)*v(124)
    v(158)=ds(5,1)*v(109)+ds(5,2)*v(112)+ds(5,3)*v(115)+ds(5,4)*v(118)+ds(5,5)*v(121)+ds(5,6)*v(124)
    v(157)=ds(4,1)*v(109)+ds(4,2)*v(112)+ds(4,3)*v(115)+ds(4,4)*v(118)+ds(4,5)*v(121)+ds(4,6)*v(124)
    v(156)=ds(3,1)*v(109)+ds(3,2)*v(112)+ds(3,4)*v(118)+ds(3,5)*v(121)+ds(3,6)*v(124)
    v(155)=ds(2,1)*v(109)+ds(2,3)*v(115)+ds(2,4)*v(118)+ds(2,5)*v(121)+ds(2,6)*v(124)
    v(108)=f(1,1)*nl(1)
    v(153)=ds(6,1)*v(108)+ds(6,2)*v(111)+ds(6,3)*v(114)+ds(6,4)*v(117)+ds(6,5)*v(120)+ds(6,6)*v(123)
    v(152)=ds(5,1)*v(108)+ds(5,2)*v(111)+ds(5,3)*v(114)+ds(5,4)*v(117)+ds(5,5)*v(120)+ds(5,6)*v(123)
    v(151)=ds(4,1)*v(108)+ds(4,2)*v(111)+ds(4,3)*v(114)+ds(4,4)*v(117)+ds(4,5)*v(120)+ds(4,6)*v(123)
    v(150)=ds(3,1)*v(108)+ds(3,2)*v(111)+ds(3,4)*v(117)+ds(3,5)*v(120)+ds(3,6)*v(123)
    v(149)=ds(2,1)*v(108)+ds(2,3)*v(114)+ds(2,4)*v(117)+ds(2,5)*v(120)+ds(2,6)*v(123)
    v(145)=ds(1,1)*v(108)
    v(104)=v(127)*v(145)+v(130)*v(146)+v(133)*v(147)
    v(105)=v(128)*v(145)+v(131)*v(146)+v(134)*v(147)
    v(106)=v(128)*v(166)+v(131)*v(167)+v(134)*v(168)
    kernel(1,1)=v(144)+v(126)*(v(145)+v(148))+v(129)*(v(146)+v(149))+v(132)*(v(147)+v(150))+v(135)*v(151)+v(138)*v(152)+v&
         &(141)*v(153)
    kernel(1,2)=v(104)+v(126)*v(154)+v(129)*v(155)+v(132)*v(156)+v(135)*v(157)+v(138)*v(158)+v(141)*v(159)
    kernel(1,3)=v(105)+v(126)*v(160)+v(129)*v(161)+v(132)*v(162)+v(135)*v(163)+v(138)*v(164)+v(141)*v(165)
    kernel(2,1)=v(104)+v(127)*v(148)+v(130)*v(149)+v(133)*v(150)+v(136)*v(151)+v(139)*v(152)+v(142)*v(153)
    kernel(2,2)=v(144)+v(136)*v(157)+v(139)*v(158)+v(142)*v(159)+v(127)*(v(154)+v(166))+v(130)*(v(155)+v(167))+v(133)*(v&
         &(156)+v(168))
    kernel(2,3)=v(106)+v(127)*v(160)+v(130)*v(161)+v(133)*v(162)+v(136)*v(163)+v(139)*v(164)+v(142)*v(165)
    kernel(3,1)=v(105)+v(128)*v(148)+v(131)*v(149)+v(134)*v(150)+v(137)*v(151)+v(140)*v(152)+v(143)*v(153)
    kernel(3,2)=v(106)+v(128)*v(154)+v(131)*v(155)+v(134)*v(156)+v(137)*v(157)+v(140)*v(158)+v(143)*v(159)
    kernel(3,3)=v(144)+v(128)*(ds(1,1)*v(110)+v(160))+v(131)*(ds(2,2)*v(113)+v(161))+v(134)*(ds(3,3)*v(116)+v(162))+v(137&
         &)*v(163)+v(140)*v(164)+v(143)*v(165)
  end subroutine mls_3d

  subroutine mls_fbarforce(dck,s,nucleus)
    implicit none
    double precision v(54),dck(6,3),s(6),nucleus(3)
    v(49)=s(3)/2d0
    v(48)=s(2)/2d0
    v(47)=s(1)/2d0
    nucleus(1)=dck(4,1)*s(4)+dck(5,1)*s(5)+dck(6,1)*s(6)+dck(1,1)*v(47)+dck(2,1)*v(48)+dck(3,1)*v(49)
    nucleus(2)=dck(4,2)*s(4)+dck(5,2)*s(5)+dck(6,2)*s(6)+dck(1,2)*v(47)+dck(2,2)*v(48)+dck(3,2)*v(49)
    nucleus(3)=dck(4,3)*s(4)+dck(5,3)*s(5)+dck(6,3)*s(6)+dck(1,3)*v(47)+dck(2,3)*v(48)+dck(3,3)*v(49)
  end subroutine mls_fbarforce

  subroutine mls_fbarforce_2d(dck,s,nucleus)
    implicit none
    double precision v(28),dck(3,2),s(3),nucleus(2)
    v(23)=s(2)/2d0
    v(22)=s(1)/2d0
    nucleus(1)=dck(3,1)*s(3)+dck(1,1)*v(22)+dck(2,1)*v(23)
    nucleus(2)=dck(3,2)*s(3)+dck(1,2)*v(22)+dck(2,2)*v(23)
  end subroutine mls_fbarforce_2d

  SUBROUTINE mls_fbar(dck,dcl,d2c,s,ds,kernel)
    IMPLICIT NONE
    DOUBLE PRECISION v(303),dck(6,3),dcl(6,3),d2c(6,3,3),s(6),ds(6,6),kernel(3,3)
    v(280)=s(3)/2d0
    v(279)=s(2)/2d0
    v(278)=s(1)/2d0
    v(277)=dcl(3,3)/2d0
    v(276)=dcl(2,3)/2d0
    v(275)=dcl(1,3)/2d0
    v(274)=ds(4,3)*v(277)
    v(273)=ds(4,2)*v(276)
    v(272)=ds(4,1)*v(275)
    v(298)=dcl(4,3)*ds(4,4)+dcl(5,3)*ds(4,5)+dcl(6,3)*ds(4,6)+v(272)+v(273)+v(274)
    v(271)=dcl(6,3)/2d0
    v(270)=dcl(5,3)/2d0
    v(269)=dcl(4,3)/2d0
    v(268)=dcl(3,3)/4d0
    v(267)=dcl(2,3)/4d0
    v(266)=dcl(1,3)/4d0
    v(265)=ds(1,6)*v(271)
    v(264)=ds(1,5)*v(270)
    v(263)=ds(1,4)*v(269)
    v(262)=ds(1,3)*v(268)
    v(261)=ds(1,2)*v(267)
    v(260)=ds(1,1)*v(266)
    v(297)=v(260)+v(261)+v(262)+v(263)+v(264)+v(265)
    v(259)=dcl(3,2)/2d0
    v(258)=dcl(2,2)/2d0
    v(257)=dcl(1,2)/2d0
    v(256)=ds(4,3)*v(259)
    v(255)=ds(4,2)*v(258)
    v(254)=ds(4,1)*v(257)
    v(292)=dcl(4,2)*ds(4,4)+dcl(5,2)*ds(4,5)+dcl(6,2)*ds(4,6)+v(254)+v(255)+v(256)
    v(253)=dcl(6,2)/2d0
    v(252)=dcl(5,2)/2d0
    v(251)=dcl(4,2)/2d0
    v(250)=dcl(3,2)/4d0
    v(249)=dcl(2,2)/4d0
    v(248)=dcl(1,2)/4d0
    v(247)=ds(1,6)*v(253)
    v(246)=ds(1,5)*v(252)
    v(245)=ds(1,4)*v(251)
    v(244)=ds(1,3)*v(250)
    v(243)=ds(1,2)*v(249)
    v(242)=ds(1,1)*v(248)
    v(291)=v(242)+v(243)+v(244)+v(245)+v(246)+v(247)
    v(241)=dcl(3,1)/2d0
    v(240)=dcl(2,1)/2d0
    v(239)=dcl(1,1)/2d0
    v(238)=ds(4,3)*v(241)
    v(237)=ds(4,2)*v(240)
    v(236)=ds(4,1)*v(239)
    v(286)=dcl(4,1)*ds(4,4)+dcl(5,1)*ds(4,5)+dcl(6,1)*ds(4,6)+v(236)+v(237)+v(238)
    v(235)=dcl(6,1)/2d0
    v(234)=dcl(5,1)/2d0
    v(233)=dcl(4,1)/2d0
    v(232)=dcl(3,1)/4d0
    v(231)=dcl(2,1)/4d0
    v(230)=dcl(1,1)/4d0
    v(229)=ds(1,6)*v(235)
    v(228)=ds(1,5)*v(234)
    v(227)=ds(1,4)*v(233)
    v(226)=ds(1,3)*v(232)
    v(225)=ds(1,2)*v(231)
    v(224)=ds(1,1)*v(230)
    v(285)=v(224)+v(225)+v(226)+v(227)+v(228)+v(229)
    v(148)=ds(2,1)*v(230)
    v(149)=ds(2,2)*v(231)
    v(150)=ds(2,3)*v(232)
    v(151)=ds(2,4)*v(233)
    v(152)=ds(2,5)*v(234)
    v(153)=ds(2,6)*v(235)
    v(281)=v(148)+v(149)+v(150)+v(151)+v(152)+v(153)
    v(154)=ds(3,1)*v(230)
    v(155)=ds(3,2)*v(231)
    v(156)=ds(3,3)*v(232)
    v(157)=ds(3,4)*v(233)
    v(158)=ds(3,5)*v(234)
    v(159)=ds(3,6)*v(235)
    v(282)=v(154)+v(155)+v(156)+v(157)+v(158)+v(159)
    v(163)=ds(5,1)*v(239)
    v(164)=ds(5,2)*v(240)
    v(165)=ds(5,3)*v(241)
    v(283)=dcl(4,1)*ds(5,4)+dcl(5,1)*ds(5,5)+dcl(6,1)*ds(5,6)+v(163)+v(164)+v(165)
    v(166)=ds(6,1)*v(239)
    v(167)=ds(6,2)*v(240)
    v(168)=ds(6,3)*v(241)
    v(284)=dcl(4,1)*ds(6,4)+dcl(5,1)*ds(6,5)+dcl(6,1)*ds(6,6)+v(166)+v(167)+v(168)
    v(175)=ds(2,1)*v(248)
    v(176)=ds(2,2)*v(249)
    v(177)=ds(2,3)*v(250)
    v(178)=ds(2,4)*v(251)
    v(179)=ds(2,5)*v(252)
    v(180)=ds(2,6)*v(253)
    v(287)=v(175)+v(176)+v(177)+v(178)+v(179)+v(180)
    v(181)=ds(3,1)*v(248)
    v(182)=ds(3,2)*v(249)
    v(183)=ds(3,3)*v(250)
    v(184)=ds(3,4)*v(251)
    v(185)=ds(3,5)*v(252)
    v(186)=ds(3,6)*v(253)
    v(288)=v(181)+v(182)+v(183)+v(184)+v(185)+v(186)
    v(190)=ds(5,1)*v(257)
    v(191)=ds(5,2)*v(258)
    v(192)=ds(5,3)*v(259)
    v(289)=dcl(4,2)*ds(5,4)+dcl(5,2)*ds(5,5)+dcl(6,2)*ds(5,6)+v(190)+v(191)+v(192)
    v(193)=ds(6,1)*v(257)
    v(194)=ds(6,2)*v(258)
    v(195)=ds(6,3)*v(259)
    v(290)=dcl(4,2)*ds(6,4)+dcl(5,2)*ds(6,5)+dcl(6,2)*ds(6,6)+v(193)+v(194)+v(195)
    v(202)=ds(2,1)*v(266)
    v(203)=ds(2,2)*v(267)
    v(204)=ds(2,3)*v(268)
    v(205)=ds(2,4)*v(269)
    v(206)=ds(2,5)*v(270)
    v(207)=ds(2,6)*v(271)
    v(293)=v(202)+v(203)+v(204)+v(205)+v(206)+v(207)
    v(208)=ds(3,1)*v(266)
    v(209)=ds(3,2)*v(267)
    v(210)=ds(3,3)*v(268)
    v(211)=ds(3,4)*v(269)
    v(212)=ds(3,5)*v(270)
    v(213)=ds(3,6)*v(271)
    v(294)=v(208)+v(209)+v(210)+v(211)+v(212)+v(213)
    v(217)=ds(5,1)*v(275)
    v(218)=ds(5,2)*v(276)
    v(219)=ds(5,3)*v(277)
    v(295)=dcl(4,3)*ds(5,4)+dcl(5,3)*ds(5,5)+dcl(6,3)*ds(5,6)+v(217)+v(218)+v(219)
    v(220)=ds(6,1)*v(275)
    v(221)=ds(6,2)*v(276)
    v(222)=ds(6,3)*v(277)
    v(296)=dcl(4,3)*ds(6,4)+dcl(5,3)*ds(6,5)+dcl(6,3)*ds(6,6)+v(220)+v(221)+v(222)
    kernel(1,1)=d2c(4,1,1)*s(4)+d2c(5,1,1)*s(5)+d2c(6,1,1)*s(6)+d2c(1,1,1)*v(278)+d2c(2,1,1)*v(279)+d2c(3,1,1)*v(280)+dck(2&
         &,1)*v(281)+dck(3,1)*v(282)+dck(5,1)*v(283)+dck(6,1)*v(284)+dck(1,1)*v(285)+dck(4,1)*v(286)
    kernel(1,2)=d2c(4,1,2)*s(4)+d2c(5,1,2)*s(5)+d2c(6,1,2)*s(6)+d2c(1,1,2)*v(278)+d2c(2,1,2)*v(279)+d2c(3,1,2)*v(280)+dck(2&
         &,1)*v(287)+dck(3,1)*v(288)+dck(5,1)*v(289)+dck(6,1)*v(290)+dck(1,1)*v(291)+dck(4,1)*v(292)
    kernel(1,3)=d2c(4,1,3)*s(4)+d2c(5,1,3)*s(5)+d2c(6,1,3)*s(6)+d2c(1,1,3)*v(278)+d2c(2,1,3)*v(279)+d2c(3,1,3)*v(280)+dck(2&
         &,1)*v(293)+dck(3,1)*v(294)+dck(5,1)*v(295)+dck(6,1)*v(296)+dck(1,1)*v(297)+dck(4,1)*v(298)
    kernel(2,1)=d2c(4,2,1)*s(4)+d2c(5,2,1)*s(5)+d2c(6,2,1)*s(6)+d2c(1,2,1)*v(278)+d2c(2,2,1)*v(279)+d2c(3,2,1)*v(280)+dck(2&
         &,2)*v(281)+dck(3,2)*v(282)+dck(5,2)*v(283)+dck(6,2)*v(284)+dck(1,2)*v(285)+dck(4,2)*v(286)
    kernel(2,2)=d2c(4,2,2)*s(4)+d2c(5,2,2)*s(5)+d2c(6,2,2)*s(6)+d2c(1,2,2)*v(278)+d2c(2,2,2)*v(279)+d2c(3,2,2)*v(280)+dck(2&
         &,2)*v(287)+dck(3,2)*v(288)+dck(5,2)*v(289)+dck(6,2)*v(290)+dck(1,2)*v(291)+dck(4,2)*v(292)
    kernel(2,3)=d2c(4,2,3)*s(4)+d2c(5,2,3)*s(5)+d2c(6,2,3)*s(6)+d2c(1,2,3)*v(278)+d2c(2,2,3)*v(279)+d2c(3,2,3)*v(280)+dck(2&
         &,2)*v(293)+dck(3,2)*v(294)+dck(5,2)*v(295)+dck(6,2)*v(296)+dck(1,2)*v(297)+dck(4,2)*v(298)
    kernel(3,1)=d2c(4,3,1)*s(4)+d2c(5,3,1)*s(5)+d2c(6,3,1)*s(6)+d2c(1,3,1)*v(278)+d2c(2,3,1)*v(279)+d2c(3,3,1)*v(280)+dck(2&
         &,3)*v(281)+dck(3,3)*v(282)+dck(5,3)*v(283)+dck(6,3)*v(284)+dck(1,3)*v(285)+dck(4,3)*v(286)
    kernel(3,2)=d2c(4,3,2)*s(4)+d2c(5,3,2)*s(5)+d2c(6,3,2)*s(6)+d2c(1,3,2)*v(278)+d2c(2,3,2)*v(279)+d2c(3,3,2)*v(280)+dck(2&
         &,3)*v(287)+dck(3,3)*v(288)+dck(5,3)*v(289)+dck(6,3)*v(290)+dck(1,3)*v(291)+dck(4,3)*v(292)
    kernel(3,3)=d2c(4,3,3)*s(4)+d2c(5,3,3)*s(5)+d2c(6,3,3)*s(6)+d2c(1,3,3)*v(278)+d2c(2,3,3)*v(279)+d2c(3,3,3)*v(280)+dck(2&
         &,3)*v(293)+dck(3,3)*v(294)+dck(5,3)*v(295)+dck(6,3)*v(296)+dck(1,3)*v(297)+dck(4,3)*v(298)
  END SUBROUTINE mls_fbar

  SUBROUTINE mls_fbar_2d(dck,dcl,d2c,s,ds,kernel)
    IMPLICIT NONE
    DOUBLE PRECISION v(86),dck(3,2),dcl(3,2),d2c(3,2,2),s(3),ds(3,3),kernel(2,2)
    v(75)=s(2)/2d0
    v(74)=s(1)/2d0
    v(73)=ds(3,2)/2d0
    v(72)=ds(3,1)/2d0
    v(71)=dcl(3,2)/2d0
    v(70)=dcl(2,2)/4d0
    v(69)=dcl(1,2)/4d0
    v(68)=ds(1,3)*v(71)
    v(67)=ds(1,2)*v(70)
    v(66)=ds(1,1)*v(69)
    v(81)=v(66)+v(67)+v(68)
    v(65)=dcl(2,1)*v(73)
    v(64)=dcl(1,1)*v(72)
    v(78)=dcl(3,1)*ds(3,3)+v(64)+v(65)
    v(63)=dcl(3,1)/2d0
    v(62)=dcl(2,1)/4d0
    v(61)=dcl(1,1)/4d0
    v(60)=ds(1,3)*v(63)
    v(59)=ds(1,2)*v(62)
    v(58)=ds(1,1)*v(61)
    v(77)=v(58)+v(59)+v(60)
    v(44)=ds(2,1)*v(61)
    v(45)=ds(2,2)*v(62)
    v(46)=ds(2,3)*v(63)
    v(76)=v(44)+v(45)+v(46)
    v(52)=ds(2,1)*v(69)
    v(53)=ds(2,2)*v(70)
    v(54)=ds(2,3)*v(71)
    v(79)=v(52)+v(53)+v(54)
    v(55)=dcl(1,2)*v(72)
    v(56)=dcl(2,2)*v(73)
    v(80)=dcl(3,2)*ds(3,3)+v(55)+v(56)
    kernel(1,1)=d2c(3,1,1)*s(3)+d2c(1,1,1)*v(74)+d2c(2,1,1)*v(75)+dck(2,1)*v(76)+dck(1,1)*v(77)+dck(3,1)*v(78)
    kernel(1,2)=d2c(3,1,2)*s(3)+d2c(1,1,2)*v(74)+d2c(2,1,2)*v(75)+dck(2,1)*v(79)+dck(3,1)*v(80)+dck(1,1)*v(81)
    kernel(2,1)=d2c(3,2,1)*s(3)+d2c(1,2,1)*v(74)+d2c(2,2,1)*v(75)+dck(2,2)*v(76)+dck(1,2)*v(77)+dck(3,2)*v(78)
    kernel(2,2)=d2c(3,2,2)*s(3)+d2c(1,2,2)*v(74)+d2c(2,2,2)*v(75)+dck(2,2)*v(79)+dck(3,2)*v(80)+dck(1,2)*v(81)
  END SUBROUTINE mls_fbar_2d


  SUBROUTINE mls_classicaldc(nk,f,dc)
    IMPLICIT NONE
    DOUBLE PRECISION v(57),nk(3),f(3,3),dc(6,3)
    v(52)=2d0*nk(3)
    v(51)=2d0*nk(2)
    v(50)=2d0*nk(1)
    dc(1,1)=f(1,1)*v(50)
    dc(1,2)=f(2,1)*v(50)
    dc(1,3)=f(3,1)*v(50)
    dc(2,1)=f(1,2)*v(51)
    dc(2,2)=f(2,2)*v(51)
    dc(2,3)=f(3,2)*v(51)
    dc(3,1)=f(1,3)*v(52)
    dc(3,2)=f(2,3)*v(52)
    dc(3,3)=f(3,3)*v(52)
    dc(4,1)=f(1,2)*nk(1)+f(1,1)*nk(2)
    dc(4,2)=f(2,2)*nk(1)+f(2,1)*nk(2)
    dc(4,3)=f(3,2)*nk(1)+f(3,1)*nk(2)
    dc(5,1)=f(1,3)*nk(1)+f(1,1)*nk(3)
    dc(5,2)=f(2,3)*nk(1)+f(2,1)*nk(3)
    dc(5,3)=f(3,3)*nk(1)+f(3,1)*nk(3)
    dc(6,1)=f(1,3)*nk(2)+f(1,2)*nk(3)
    dc(6,2)=f(2,3)*nk(2)+f(2,2)*nk(3)
    dc(6,3)=f(3,3)*nk(2)+f(3,2)*nk(3)
  END SUBROUTINE mls_classicaldc


  SUBROUTINE mls_classicaldc_2d(nk,f,dc)
    IMPLICIT NONE
    DOUBLE PRECISION v(26),nk(2),f(2,2),dc(3,2)
    v(21)=2d0*nk(2)
    v(20)=2d0*nk(1)
    dc(1,1)=f(1,1)*v(20)
    dc(1,2)=f(2,1)*v(20)
    dc(2,1)=f(1,2)*v(21)
    dc(2,2)=f(2,2)*v(21)
    dc(3,1)=f(1,2)*nk(1)+f(1,1)*nk(2)
    dc(3,2)=f(2,2)*nk(1)+f(2,1)*nk(2)
  END SUBROUTINE mls_classicaldc_2d

  SUBROUTINE mls_classicaldc2(nk,nl,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(78),nk(3),nl(3),d2c(6,3,3)
    v(73)=nk(3)*nl(2)+nk(2)*nl(3)
    v(72)=nk(3)*nl(1)+nk(1)*nl(3)
    v(71)=nk(2)*nl(1)+nk(1)*nl(2)
    v(70)=2d0*nk(3)*nl(3)
    v(69)=2d0*nk(2)*nl(2)
    v(68)=2d0*nk(1)*nl(1)
    d2c(1,1,1)=v(68)
    d2c(1,1,2)=0d0
    d2c(1,1,3)=0d0
    d2c(1,2,1)=0d0
    d2c(1,2,2)=v(68)
    d2c(1,2,3)=0d0
    d2c(1,3,1)=0d0
    d2c(1,3,2)=0d0
    d2c(1,3,3)=v(68)
    d2c(2,1,1)=v(69)
    d2c(2,1,2)=0d0
    d2c(2,1,3)=0d0
    d2c(2,2,1)=0d0
    d2c(2,2,2)=v(69)
    d2c(2,2,3)=0d0
    d2c(2,3,1)=0d0
    d2c(2,3,2)=0d0
    d2c(2,3,3)=v(69)
    d2c(3,1,1)=v(70)
    d2c(3,1,2)=0d0
    d2c(3,1,3)=0d0
    d2c(3,2,1)=0d0
    d2c(3,2,2)=v(70)
    d2c(3,2,3)=0d0
    d2c(3,3,1)=0d0
    d2c(3,3,2)=0d0
    d2c(3,3,3)=v(70)
    d2c(4,1,1)=v(71)
    d2c(4,1,2)=0d0
    d2c(4,1,3)=0d0
    d2c(4,2,1)=0d0
    d2c(4,2,2)=v(71)
    d2c(4,2,3)=0d0
    d2c(4,3,1)=0d0
    d2c(4,3,2)=0d0
    d2c(4,3,3)=v(71)
    d2c(5,1,1)=v(72)
    d2c(5,1,2)=0d0
    d2c(5,1,3)=0d0
    d2c(5,2,1)=0d0
    d2c(5,2,2)=v(72)
    d2c(5,2,3)=0d0
    d2c(5,3,1)=0d0
    d2c(5,3,2)=0d0
    d2c(5,3,3)=v(72)
    d2c(6,1,1)=v(73)
    d2c(6,1,2)=0d0
    d2c(6,1,3)=0d0
    d2c(6,2,1)=0d0
    d2c(6,2,2)=v(73)
    d2c(6,2,3)=0d0
    d2c(6,3,1)=0d0
    d2c(6,3,2)=0d0
    d2c(6,3,3)=v(73)
  END SUBROUTINE mls_classicaldc2


  SUBROUTINE mls_classicaldc2_2d(nk,nl,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(28),nk(2),nl(2),d2c(3,2,2)
    v(23)=nk(2)*nl(1)+nk(1)*nl(2)
    v(22)=2d0*nk(2)*nl(2)
    v(21)=2d0*nk(1)*nl(1)
    d2c(1,1,1)=v(21)
    d2c(1,1,2)=0d0
    d2c(1,2,1)=0d0
    d2c(1,2,2)=v(21)
    d2c(2,1,1)=v(22)
    d2c(2,1,2)=0d0
    d2c(2,2,1)=0d0
    d2c(2,2,2)=v(22)
    d2c(3,1,1)=v(23)
    d2c(3,1,2)=0d0
    d2c(3,2,1)=0d0
    d2c(3,2,2)=v(23)
  END SUBROUTINE mls_classicaldc2_2d


!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           5 Feb 22 20:02:45  *
!**************************************************************
! User     : Full professional version
! Notebook : mls_combinedfbar0
! Evaluation time                 : 3 s     Mode  : Optimal
! Number of formulae              : 2       Method: Automatic
! Subroutine                      : mls_combinedfbar0 size: 203
! Total size of Mathematica  code : 203 subexpressions
! Total size of Fortran code      : 472 bytes

!******************* S U B R O U T I N E **********************
  SUBROUTINE mls_combinedfbar0(c,cb,cs)
    IMPLICIT NONE
    DOUBLE PRECISION v(32),c(6),cb(6),cs(6)
    v(27)=((cb(3)*cb(4)**2+cb(5)*(cb(2)*cb(5)-2d0*cb(4)*cb(6))+cb(1)*(-(cb(2)*cb(3))+cb(6)**2))/(c(3)*c(4)**2+c(5)*(c(2)*c&
         &(5)-2d0*c(4)*c(6))+c(1)*(-(c(2)*c(3))+c(6)**2)))**0.3333333333333333d0
    cs(1)=c(1)*v(27)
    cs(2)=c(2)*v(27)
    cs(3)=c(3)*v(27)
    cs(4)=c(4)*v(27)
    cs(5)=c(5)*v(27)
    cs(6)=c(6)*v(27)
  END SUBROUTINE mls_combinedfbar0

  SUBROUTINE mls_combinedfbar0_2d(c,cb,cs)
    IMPLICIT NONE
    DOUBLE PRECISION v(20),c(3),cb(3),cs(3)
    v(15)=SQRT((cb(1)*cb(2)-cb(3)**2)/(c(1)*c(2)-c(3)**2))
    cs(1)=c(1)*v(15)
    cs(2)=c(2)*v(15)
    cs(3)=c(3)*v(15)
  END SUBROUTINE mls_combinedfbar0_2d

  SUBROUTINE mls_combinedfbar1(c,cb,dc,dcb,dcs)
    IMPLICIT NONE
    DOUBLE PRECISION v(166),c(6),cb(6),dc(6,3),dcb(6,3),dcs(6,3)
    v(161)=cb(2)*dcb(3,3)
    v(160)=cb(3)*dcb(2,3)
    v(159)=-(cb(2)*dcb(5,3))
    v(158)=cb(6)*dcb(4,3)
    v(157)=c(2)*dc(3,3)
    v(156)=c(3)*dc(2,3)
    v(155)=-(c(2)*dc(5,3))
    v(154)=c(6)*dc(4,3)
    v(153)=2d0*dcb(6,3)
    v(152)=2d0*dc(6,3)
    v(151)=cb(2)*dcb(3,2)
    v(150)=cb(3)*dcb(2,2)
    v(149)=-(cb(2)*dcb(5,2))
    v(148)=cb(6)*dcb(4,2)
    v(147)=c(2)*dc(3,2)
    v(146)=c(3)*dc(2,2)
    v(145)=-(c(2)*dc(5,2))
    v(144)=c(6)*dc(4,2)
    v(143)=2d0*dcb(6,2)
    v(142)=cb(4)**2
    v(141)=cb(5)**2
    v(140)=2d0*dc(6,2)
    v(139)=c(4)**2
    v(138)=c(5)**2
    v(135)=cb(2)*dcb(3,1)
    v(134)=cb(3)*dcb(2,1)
    v(133)=-(cb(2)*dcb(5,1))
    v(132)=cb(6)*dcb(4,1)
    v(131)=-(dcb(3,1)*v(142))
    v(130)=-(dcb(2,1)*v(141))
    v(129)=c(2)*dc(3,1)
    v(128)=c(3)*dc(2,1)
    v(127)=-(c(2)*dc(5,1))
    v(126)=c(6)*dc(4,1)
    v(125)=-(dc(3,1)*v(139))
    v(124)=-(dc(2,1)*v(138))
    v(123)=2d0*dcb(6,1)
    v(122)=cb(4)*cb(6)
    v(121)=cb(3)*cb(4)
    v(120)=2d0*dc(6,1)
    v(119)=c(4)*c(6)
    v(118)=c(3)*c(4)
    v(117)=cb(2)*cb(3)-cb(6)**2
    v(116)=-(cb(5)*cb(6))+v(121)
    v(115)=-(cb(2)*cb(5))+v(122)
    v(114)=c(2)*c(3)-c(6)**2
    v(113)=-(c(5)*c(6))+v(118)
    v(112)=-(c(2)*c(5))+v(119)
    v(78)=c(5)*v(112)-c(4)*v(113)+c(1)*v(114)
    v(88)=1d0/v(78)**2
    v(79)=cb(5)*v(115)-cb(4)*v(116)+cb(1)*v(117)
    v(137)=-(v(79)*v(88))
    v(77)=v(79)/v(78)
    v(87)=1d0/v(77)**0.6666666666666666d0
    v(136)=v(87)/3d0
    v(71)=v(77)**0.3333333333333333d0
    v(86)=v(136)*((dc(5,1)*v(112)-dc(4,1)*v(113)+dc(1,1)*v(114)-dc(4,1)*v(118)+dc(5,1)*v(119)+v(124)+v(125)+c(5)*(c(4)*v&
         &(120)+v(126)+v(127))+c(1)*(-(c(6)*v(120))+v(128)+v(129)))*v(137)+(dcb(5,1)*v(115)-dcb(4,1)*v(116)+dcb(1,1)*v(117)-dcb(4&
         &,1)*v(121)+dcb(5,1)*v(122)+v(130)+v(131)+cb(5)*(cb(4)*v(123)+v(132)+v(133))+cb(1)*(-(cb(6)*v(123))+v(134)+v(135)))/v(78&
         &))
    v(90)=v(136)*(v(137)*(dc(5,2)*v(112)-dc(4,2)*v(113)+dc(1,2)*v(114)-dc(4,2)*v(118)+dc(5,2)*v(119)-dc(2,2)*v(138)-dc(3,2&
         &)*v(139)+c(5)*(c(4)*v(140)+v(144)+v(145))+c(1)*(-(c(6)*v(140))+v(146)+v(147)))+(dcb(5,2)*v(115)-dcb(4,2)*v(116)+dcb(1,2&
         &)*v(117)-dcb(4,2)*v(121)+dcb(5,2)*v(122)-dcb(2,2)*v(141)-dcb(3,2)*v(142)+cb(5)*(cb(4)*v(143)+v(148)+v(149))+cb(1)*(-(cb&
         &(6)*v(143))+v(150)+v(151)))/v(78))
    v(92)=v(136)*(v(137)*(dc(5,3)*v(112)-dc(4,3)*v(113)+dc(1,3)*v(114)-dc(4,3)*v(118)+dc(5,3)*v(119)-dc(2,3)*v(138)-dc(3,3&
         &)*v(139)+c(5)*(c(4)*v(152)+v(154)+v(155))+c(1)*(-(c(6)*v(152))+v(156)+v(157)))+(dcb(5,3)*v(115)-dcb(4,3)*v(116)+dcb(1,3&
         &)*v(117)-dcb(4,3)*v(121)+dcb(5,3)*v(122)-dcb(2,3)*v(141)-dcb(3,3)*v(142)+cb(5)*(cb(4)*v(153)+v(158)+v(159))+cb(1)*(-(cb&
         &(6)*v(153))+v(160)+v(161)))/v(78))
    dcs(1,1)=dc(1,1)*v(71)+c(1)*v(86)
    dcs(1,2)=dc(1,2)*v(71)+c(1)*v(90)
    dcs(1,3)=dc(1,3)*v(71)+c(1)*v(92)
    dcs(2,1)=dc(2,1)*v(71)+c(2)*v(86)
    dcs(2,2)=dc(2,2)*v(71)+c(2)*v(90)
    dcs(2,3)=dc(2,3)*v(71)+c(2)*v(92)
    dcs(3,1)=dc(3,1)*v(71)+c(3)*v(86)
    dcs(3,2)=dc(3,2)*v(71)+c(3)*v(90)
    dcs(3,3)=dc(3,3)*v(71)+c(3)*v(92)
    dcs(4,1)=dc(4,1)*v(71)+c(4)*v(86)
    dcs(4,2)=dc(4,2)*v(71)+c(4)*v(90)
    dcs(4,3)=dc(4,3)*v(71)+c(4)*v(92)
    dcs(5,1)=dc(5,1)*v(71)+c(5)*v(86)
    dcs(5,2)=dc(5,2)*v(71)+c(5)*v(90)
    dcs(5,3)=dc(5,3)*v(71)+c(5)*v(92)
    dcs(6,1)=dc(6,1)*v(71)+c(6)*v(86)
    dcs(6,2)=dc(6,2)*v(71)+c(6)*v(90)
    dcs(6,3)=dc(6,3)*v(71)+c(6)*v(92)
  END SUBROUTINE mls_combinedfbar1

  SUBROUTINE mls_combinedfbar1_2d(c,cb,dc,dcb,dcs)
    IMPLICIT NONE
    DOUBLE PRECISION v(62),c(3),cb(3),dc(3,2),dcb(3,2),dcs(3,2)
    v(57)=c(1)*dc(2,2)
    v(56)=c(2)*dc(1,2)
    v(55)=cb(1)*dcb(2,2)
    v(54)=cb(2)*dcb(1,2)
    v(53)=(-2d0)*c(3)
    v(51)=(-2d0)*cb(3)
    v(49)=c(2)*dc(1,1)+c(1)*dc(2,1)+dc(3,1)*v(53)
    v(48)=cb(2)*dcb(1,1)+cb(1)*dcb(2,1)+dcb(3,1)*v(51)
    v(47)=cb(1)*cb(2)-cb(3)**2
    v(46)=c(1)*c(2)-c(3)**2
    v(36)=1d0/v(46)**2
    v(52)=-(v(36)*v(47))
    v(31)=v(47)/v(46)
    v(35)=1d0/sqrt(v(31))
    v(50)=v(35)/2d0
    v(28)=sqrt(v(31))
    v(34)=v(50)*(v(48)/v(46)+v(49)*v(52))
    v(38)=v(50)*((dcb(3,2)*v(51)+v(54)+v(55))/v(46)+v(52)*(dc(3,2)*v(53)+v(56)+v(57)))
    dcs(1,1)=dc(1,1)*v(28)+c(1)*v(34)
    dcs(1,2)=dc(1,2)*v(28)+c(1)*v(38)
    dcs(2,1)=dc(2,1)*v(28)+c(2)*v(34)
    dcs(2,2)=dc(2,2)*v(28)+c(2)*v(38)
    dcs(3,1)=dc(3,1)*v(28)+c(3)*v(34)
    dcs(3,2)=dc(3,2)*v(28)+c(3)*v(38)
  END SUBROUTINE mls_combinedfbar1_2d

  SUBROUTINE mls_combinedfbar2_2d(c1,c0,dck1,dcl1,dck0,dcl0,d2c1,d2c0,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(202),c1(3),c0(3),dck1(3,2),dcl1(3,2),dck0(3,2),dcl0(3,2),d2c1(3,2,2),d2c0(3,2,2),d2c(3,2,2&
         &)
    v(197)=dck1(1,1)*dcl1(2,2)
    v(196)=dck1(2,1)*dcl1(1,2)
    v(195)=dck0(1,1)*dcl0(2,2)
    v(194)=dck0(2,1)*dcl0(1,2)
    v(193)=c0(1)*d2c0(2,1,2)
    v(192)=c0(2)*d2c0(1,1,2)
    v(191)=(-2d0)*dcl1(3,2)
    v(188)=(-2d0)*dcl0(3,2)
    v(187)=dck1(3,2)*v(191)
    v(186)=dck1(1,2)*dcl1(2,2)
    v(185)=dck1(2,2)*dcl1(1,2)
    v(184)=dck0(3,2)*v(188)
    v(183)=dck0(1,2)*dcl0(2,2)
    v(182)=dck0(2,2)*dcl0(1,2)
    v(181)=c0(1)*d2c0(2,2,2)
    v(180)=c0(2)*d2c0(1,2,2)
    v(178)=dck1(1,1)*dcl1(2,1)
    v(177)=dck1(2,1)*dcl1(1,1)
    v(176)=dck0(1,1)*dcl0(2,1)
    v(175)=dck0(2,1)*dcl0(1,1)
    v(174)=c0(1)*d2c0(2,1,1)
    v(173)=c0(2)*d2c0(1,1,1)
    v(172)=(-2d0)*dcl1(3,1)
    v(169)=(-2d0)*dcl0(3,1)
    v(168)=dck1(3,2)*v(172)
    v(167)=dck1(1,2)*dcl1(2,1)
    v(166)=dck1(2,2)*dcl1(1,1)
    v(165)=dck0(3,2)*v(169)
    v(164)=dck0(1,2)*dcl0(2,1)
    v(163)=dck0(2,2)*dcl0(1,1)
    v(162)=c0(1)*d2c0(2,2,1)
    v(161)=c0(2)*d2c0(1,2,1)
    v(160)=c0(1)*dcl0(2,2)
    v(159)=c0(2)*dcl0(1,2)
    v(158)=c0(1)*dcl0(2,1)
    v(157)=c0(2)*dcl0(1,1)
    v(155)=c1(1)*dcl1(2,2)
    v(154)=c1(2)*dcl1(1,2)
    v(153)=c1(1)*dcl1(2,1)
    v(152)=c1(2)*dcl1(1,1)
    v(150)=d2c1(3,2,2)
    v(149)=d2c1(3,2,1)
    v(148)=d2c1(3,1,2)
    v(147)=d2c1(3,1,1)
    v(146)=d2c1(2,2,2)
    v(145)=d2c1(2,2,1)
    v(144)=d2c1(2,1,2)
    v(143)=d2c1(2,1,1)
    v(142)=d2c1(1,2,2)
    v(141)=d2c1(1,2,1)
    v(140)=d2c1(1,1,2)
    v(139)=d2c1(1,1,1)
    v(137)=c0(1)*dck0(2,2)
    v(136)=c0(2)*dck0(1,2)
    v(135)=(-2d0)*c0(3)
    v(134)=c0(2)*dck0(1,1)+c0(1)*dck0(2,1)+dck0(3,1)*v(135)
    v(133)=c1(1)*dck1(2,2)
    v(132)=c1(2)*dck1(1,2)
    v(131)=(-2d0)*c1(3)
    v(130)=c1(2)*dck1(1,1)+c1(1)*dck1(2,1)+dck1(3,1)*v(131)
    v(129)=c0(1)*c0(2)-c0(3)**2
    v(128)=c1(1)*c1(2)-c1(3)**2
    v(94)=1d0/v(128)**3
    v(156)=(-2d0)*v(94)
    v(80)=1d0/v(128)**2
    v(138)=-(v(129)*v(80))
    v(75)=v(129)/v(128)
    v(101)=1d0/v(75)**0.15d1
    v(179)=(-0.5d0)*v(101)
    v(79)=1d0/sqrt(v(75))
    v(151)=v(79)/2d0
    v(72)=sqrt(v(75))
    v(111)=dck1(3,2)*v(131)+v(132)+v(133)
    v(107)=v(134)/v(128)+v(130)*v(138)
    v(110)=dck0(3,2)*v(135)+v(136)+v(137)
    v(112)=v(110)/v(128)+v(111)*v(138)
    v(78)=v(107)*v(151)
    v(82)=v(112)*v(151)
    v(89)=dcl1(3,1)*v(131)+v(152)+v(153)
    v(90)=dcl1(3,2)*v(131)+v(154)+v(155)
    v(91)=-(v(80)*v(89))
    v(92)=-(v(80)*v(90))
    v(93)=v(156)*v(89)
    v(171)=-(v(129)*v(93))
    v(95)=v(156)*v(90)
    v(190)=-(v(129)*v(95))
    v(96)=dcl0(3,1)*v(135)+v(157)+v(158)
    v(170)=-(v(80)*v(96))
    v(97)=dcl0(3,2)*v(135)+v(159)+v(160)
    v(189)=-(v(80)*v(97))
    v(98)=v(129)*v(91)+v(96)/v(128)
    v(99)=v(129)*v(92)+v(97)/v(128)
    v(100)=v(179)*v(98)
    v(113)=(v(100)*v(112)+v(79)*((d2c0(3,2,1)*v(135)+v(161)+v(162)+v(163)+v(164)+v(165))/v(128)+v(138)*(c1(2)*v(141)+c1(1&
         &)*v(145)+v(131)*v(149)+v(166)+v(167)+v(168))+v(111)*v(170)+v(111)*v(171)+v(110)*v(91)))/2d0
    v(108)=(v(100)*v(107)+v(79)*(v(130)*v(170)+v(130)*v(171)+(d2c0(3,1,1)*v(135)+dck0(3,1)*v(169)+v(173)+v(174)+v(175)+v&
         &(176))/v(128)+v(138)*(c1(2)*v(139)+c1(1)*v(143)+v(131)*v(147)+dck1(3,1)*v(172)+v(177)+v(178))+v(134)*v(91)))/2d0
    v(102)=v(179)*v(99)
    v(114)=(v(102)*v(112)+v(79)*((d2c0(3,2,2)*v(135)+v(180)+v(181)+v(182)+v(183)+v(184))/v(128)+v(138)*(c1(2)*v(142)+c1(1&
         &)*v(146)+v(131)*v(150)+v(185)+v(186)+v(187))+v(111)*v(189)+v(111)*v(190)+v(110)*v(92)))/2d0
    v(109)=(v(102)*v(107)+v(79)*(v(130)*v(189)+v(130)*v(190)+(d2c0(3,1,2)*v(135)+dck0(3,1)*v(188)+v(192)+v(193)+v(194)+v&
         &(195))/v(128)+v(138)*(c1(2)*v(140)+c1(1)*v(144)+v(131)*v(148)+dck1(3,1)*v(191)+v(196)+v(197))+v(134)*v(92)))/2d0
    v(103)=v(151)*v(98)
    v(104)=v(151)*v(99)
    d2c(1,1,1)=dck1(1,1)*v(103)+c1(1)*v(108)+v(139)*v(72)+dcl1(1,1)*v(78)
    d2c(1,1,2)=dck1(1,1)*v(104)+c1(1)*v(109)+v(140)*v(72)+dcl1(1,2)*v(78)
    d2c(1,2,1)=dck1(1,2)*v(103)+c1(1)*v(113)+v(141)*v(72)+dcl1(1,1)*v(82)
    d2c(1,2,2)=dck1(1,2)*v(104)+c1(1)*v(114)+v(142)*v(72)+dcl1(1,2)*v(82)
    d2c(2,1,1)=dck1(2,1)*v(103)+c1(2)*v(108)+v(143)*v(72)+dcl1(2,1)*v(78)
    d2c(2,1,2)=dck1(2,1)*v(104)+c1(2)*v(109)+v(144)*v(72)+dcl1(2,2)*v(78)
    d2c(2,2,1)=dck1(2,2)*v(103)+c1(2)*v(113)+v(145)*v(72)+dcl1(2,1)*v(82)
    d2c(2,2,2)=dck1(2,2)*v(104)+c1(2)*v(114)+v(146)*v(72)+dcl1(2,2)*v(82)
    d2c(3,1,1)=dck1(3,1)*v(103)+c1(3)*v(108)+v(147)*v(72)+dcl1(3,1)*v(78)
    d2c(3,1,2)=dck1(3,1)*v(104)+c1(3)*v(109)+v(148)*v(72)+dcl1(3,2)*v(78)
    d2c(3,2,1)=dck1(3,2)*v(103)+c1(3)*v(113)+v(149)*v(72)+dcl1(3,1)*v(82)
    d2c(3,2,2)=dck1(3,2)*v(104)+c1(3)*v(114)+v(150)*v(72)+dcl1(3,2)*v(82)
  END SUBROUTINE mls_combinedfbar2_2d

!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           9 Feb 22 10:52:10  *
!**************************************************************
! User     : Full professional version
! Notebook : mls_combinedfbar2
! Evaluation time                 : 25 s    Mode  : Optimal
! Number of formulae              : 416     Method: Automatic
! Subroutine                      : mls_combinedfbar2 size: 13062
! Total size of Mathematica  code : 13062 subexpressions
! Total size of Fortran code      : 24811 bytes

!******************* S U B R O U T I N E **********************
  SUBROUTINE mls_combinedfbar2(c1,c0,dck1,dcl1,dck0,dcl0,d2c1,d2c0,d2c)
    IMPLICIT NONE
    DOUBLE PRECISION v(785),c1(6),c0(6),dck1(6,3),dcl1(6,3),dck0(6,3),dcl0(6,3),d2c1(6,3,3),d2c0(6,3,3),d2c(6,3,3&
         &)
    v(780)=dck0(2,1)*dcl0(3,3)
    v(779)=dck0(3,1)*dcl0(2,3)
    v(778)=-(dck0(5,1)*dcl0(6,3))
    v(777)=-(dck0(6,1)*dcl0(5,3))
    v(776)=dck0(3,1)*dcl0(4,3)
    v(775)=dck0(4,1)*dcl0(3,3)
    v(774)=dck0(4,1)*dcl0(6,3)
    v(773)=-(dck0(2,1)*dcl0(5,3))
    v(772)=dck0(6,1)*dcl0(4,3)
    v(771)=-(dck0(5,1)*dcl0(2,3))
    v(770)=dck1(2,1)*dcl1(3,3)
    v(769)=dck1(3,1)*dcl1(2,3)
    v(768)=-(dck1(5,1)*dcl1(6,3))
    v(767)=-(dck1(6,1)*dcl1(5,3))
    v(766)=dck1(3,1)*dcl1(4,3)
    v(765)=dck1(4,1)*dcl1(3,3)
    v(764)=dck1(4,1)*dcl1(6,3)
    v(763)=-(dck1(2,1)*dcl1(5,3))
    v(762)=dck1(6,1)*dcl1(4,3)
    v(761)=-(dck1(5,1)*dcl1(2,3))
    v(760)=dck0(2,2)*dcl0(3,3)
    v(759)=dck0(3,2)*dcl0(2,3)
    v(758)=-(dck0(5,2)*dcl0(6,3))
    v(757)=-(dck0(6,2)*dcl0(5,3))
    v(756)=dck0(3,2)*dcl0(4,3)
    v(755)=dck0(4,2)*dcl0(3,3)
    v(754)=dck0(4,2)*dcl0(6,3)
    v(753)=-(dck0(2,2)*dcl0(5,3))
    v(752)=dck0(6,2)*dcl0(4,3)
    v(751)=-(dck0(5,2)*dcl0(2,3))
    v(750)=dck1(2,2)*dcl1(3,3)
    v(749)=dck1(3,2)*dcl1(2,3)
    v(748)=-(dck1(5,2)*dcl1(6,3))
    v(747)=-(dck1(6,2)*dcl1(5,3))
    v(746)=dck1(3,2)*dcl1(4,3)
    v(745)=dck1(4,2)*dcl1(3,3)
    v(744)=dck1(4,2)*dcl1(6,3)
    v(743)=-(dck1(2,2)*dcl1(5,3))
    v(742)=dck1(6,2)*dcl1(4,3)
    v(741)=-(dck1(5,2)*dcl1(2,3))
    v(740)=(-2d0)*dcl0(6,3)
    v(739)=(-2d0)*dcl1(6,3)
    v(736)=dck0(6,3)*v(740)
    v(735)=dck0(2,3)*dcl0(3,3)
    v(734)=dck0(3,3)*dcl0(2,3)
    v(733)=-(dck0(5,3)*dcl0(6,3))
    v(732)=-(dck0(6,3)*dcl0(5,3))
    v(731)=dck0(3,3)*dcl0(4,3)
    v(730)=dck0(4,3)*dcl0(3,3)
    v(729)=dck0(4,3)*dcl0(6,3)
    v(728)=-(dck0(2,3)*dcl0(5,3))
    v(727)=dck0(6,3)*dcl0(4,3)
    v(726)=-(dck0(5,3)*dcl0(2,3))
    v(725)=dck1(6,3)*v(739)
    v(724)=dck1(2,3)*dcl1(3,3)
    v(723)=dck1(3,3)*dcl1(2,3)
    v(722)=-(dck1(5,3)*dcl1(6,3))
    v(721)=-(dck1(6,3)*dcl1(5,3))
    v(720)=dck1(3,3)*dcl1(4,3)
    v(719)=dck1(4,3)*dcl1(3,3)
    v(718)=dck1(4,3)*dcl1(6,3)
    v(717)=-(dck1(2,3)*dcl1(5,3))
    v(716)=dck1(6,3)*dcl1(4,3)
    v(715)=-(dck1(5,3)*dcl1(2,3))
    v(714)=dck0(2,1)*dcl0(3,2)
    v(713)=dck0(3,1)*dcl0(2,2)
    v(712)=-(dck0(5,1)*dcl0(6,2))
    v(711)=-(dck0(6,1)*dcl0(5,2))
    v(710)=dck0(3,1)*dcl0(4,2)
    v(709)=dck0(4,1)*dcl0(3,2)
    v(708)=dck0(4,1)*dcl0(6,2)
    v(707)=-(dck0(2,1)*dcl0(5,2))
    v(706)=dck0(6,1)*dcl0(4,2)
    v(705)=-(dck0(5,1)*dcl0(2,2))
    v(704)=dck1(2,1)*dcl1(3,2)
    v(703)=dck1(3,1)*dcl1(2,2)
    v(702)=-(dck1(5,1)*dcl1(6,2))
    v(701)=-(dck1(6,1)*dcl1(5,2))
    v(700)=dck1(3,1)*dcl1(4,2)
    v(699)=dck1(4,1)*dcl1(3,2)
    v(698)=dck1(4,1)*dcl1(6,2)
    v(697)=-(dck1(2,1)*dcl1(5,2))
    v(696)=dck1(6,1)*dcl1(4,2)
    v(695)=-(dck1(5,1)*dcl1(2,2))
    v(694)=dck0(2,2)*dcl0(3,2)
    v(693)=dck0(3,2)*dcl0(2,2)
    v(692)=-(dck0(5,2)*dcl0(6,2))
    v(691)=-(dck0(6,2)*dcl0(5,2))
    v(690)=dck0(3,2)*dcl0(4,2)
    v(689)=dck0(4,2)*dcl0(3,2)
    v(688)=dck0(4,2)*dcl0(6,2)
    v(687)=-(dck0(2,2)*dcl0(5,2))
    v(686)=dck0(6,2)*dcl0(4,2)
    v(685)=-(dck0(5,2)*dcl0(2,2))
    v(684)=dck1(2,2)*dcl1(3,2)
    v(683)=dck1(3,2)*dcl1(2,2)
    v(682)=-(dck1(5,2)*dcl1(6,2))
    v(681)=-(dck1(6,2)*dcl1(5,2))
    v(680)=dck1(3,2)*dcl1(4,2)
    v(679)=dck1(4,2)*dcl1(3,2)
    v(678)=dck1(4,2)*dcl1(6,2)
    v(677)=-(dck1(2,2)*dcl1(5,2))
    v(676)=dck1(6,2)*dcl1(4,2)
    v(675)=-(dck1(5,2)*dcl1(2,2))
    v(674)=(-2d0)*dcl0(6,2)
    v(673)=(-2d0)*dcl1(6,2)
    v(670)=dck0(6,3)*v(674)
    v(669)=dck0(2,3)*dcl0(3,2)
    v(668)=dck0(3,3)*dcl0(2,2)
    v(667)=-(dck0(5,3)*dcl0(6,2))
    v(666)=-(dck0(6,3)*dcl0(5,2))
    v(665)=dck0(3,3)*dcl0(4,2)
    v(664)=dck0(4,3)*dcl0(3,2)
    v(663)=dck0(4,3)*dcl0(6,2)
    v(662)=-(dck0(2,3)*dcl0(5,2))
    v(661)=dck0(6,3)*dcl0(4,2)
    v(660)=-(dck0(5,3)*dcl0(2,2))
    v(659)=dck1(6,3)*v(673)
    v(658)=dck1(2,3)*dcl1(3,2)
    v(657)=dck1(3,3)*dcl1(2,2)
    v(656)=-(dck1(5,3)*dcl1(6,2))
    v(655)=-(dck1(6,3)*dcl1(5,2))
    v(654)=dck1(3,3)*dcl1(4,2)
    v(653)=dck1(4,3)*dcl1(3,2)
    v(652)=dck1(4,3)*dcl1(6,2)
    v(651)=-(dck1(2,3)*dcl1(5,2))
    v(650)=dck1(6,3)*dcl1(4,2)
    v(649)=-(dck1(5,3)*dcl1(2,2))
    v(647)=dck0(2,1)*dcl0(3,1)
    v(646)=dck0(3,1)*dcl0(2,1)
    v(645)=-(dck0(5,1)*dcl0(6,1))
    v(644)=-(dck0(6,1)*dcl0(5,1))
    v(643)=dck0(3,1)*dcl0(4,1)
    v(642)=dck0(4,1)*dcl0(3,1)
    v(641)=dck0(4,1)*dcl0(6,1)
    v(640)=-(dck0(2,1)*dcl0(5,1))
    v(639)=dck0(6,1)*dcl0(4,1)
    v(638)=-(dck0(5,1)*dcl0(2,1))
    v(637)=dck1(2,1)*dcl1(3,1)
    v(636)=dck1(3,1)*dcl1(2,1)
    v(635)=-(dck1(5,1)*dcl1(6,1))
    v(634)=-(dck1(6,1)*dcl1(5,1))
    v(633)=dck1(3,1)*dcl1(4,1)
    v(632)=dck1(4,1)*dcl1(3,1)
    v(631)=dck1(4,1)*dcl1(6,1)
    v(630)=-(dck1(2,1)*dcl1(5,1))
    v(629)=dck1(6,1)*dcl1(4,1)
    v(628)=-(dck1(5,1)*dcl1(2,1))
    v(627)=dck0(2,2)*dcl0(3,1)
    v(626)=dck0(3,2)*dcl0(2,1)
    v(625)=-(dck0(5,2)*dcl0(6,1))
    v(624)=-(dck0(6,2)*dcl0(5,1))
    v(623)=dck0(3,2)*dcl0(4,1)
    v(622)=dck0(4,2)*dcl0(3,1)
    v(621)=dck0(4,2)*dcl0(6,1)
    v(620)=-(dck0(2,2)*dcl0(5,1))
    v(619)=dck0(6,2)*dcl0(4,1)
    v(618)=-(dck0(5,2)*dcl0(2,1))
    v(617)=dck1(2,2)*dcl1(3,1)
    v(616)=dck1(3,2)*dcl1(2,1)
    v(615)=-(dck1(5,2)*dcl1(6,1))
    v(614)=-(dck1(6,2)*dcl1(5,1))
    v(613)=dck1(3,2)*dcl1(4,1)
    v(612)=dck1(4,2)*dcl1(3,1)
    v(611)=dck1(4,2)*dcl1(6,1)
    v(610)=-(dck1(2,2)*dcl1(5,1))
    v(609)=dck1(6,2)*dcl1(4,1)
    v(608)=-(dck1(5,2)*dcl1(2,1))
    v(607)=(-2d0)*dcl0(6,1)
    v(606)=(-2d0)*dcl1(6,1)
    v(603)=dck0(6,3)*v(607)
    v(602)=dck0(2,3)*dcl0(3,1)
    v(601)=dck0(3,3)*dcl0(2,1)
    v(600)=-(dck0(5,3)*dcl0(6,1))
    v(599)=-(dck0(6,3)*dcl0(5,1))
    v(598)=dck0(3,3)*dcl0(4,1)
    v(597)=dck0(4,3)*dcl0(3,1)
    v(596)=dck0(4,3)*dcl0(6,1)
    v(595)=-(dck0(2,3)*dcl0(5,1))
    v(594)=dck0(6,3)*dcl0(4,1)
    v(593)=-(dck0(5,3)*dcl0(2,1))
    v(592)=dck1(6,3)*v(606)
    v(591)=dck1(2,3)*dcl1(3,1)
    v(590)=dck1(3,3)*dcl1(2,1)
    v(589)=-(dck1(5,3)*dcl1(6,1))
    v(588)=-(dck1(6,3)*dcl1(5,1))
    v(587)=dck1(3,3)*dcl1(4,1)
    v(586)=dck1(4,3)*dcl1(3,1)
    v(585)=dck1(4,3)*dcl1(6,1)
    v(584)=-(dck1(2,3)*dcl1(5,1))
    v(583)=dck1(6,3)*dcl1(4,1)
    v(582)=-(dck1(5,3)*dcl1(2,1))
    v(581)=c0(2)*dcl0(3,3)
    v(580)=c0(3)*dcl0(2,3)
    v(579)=c0(2)*dcl0(3,2)
    v(578)=c0(3)*dcl0(2,2)
    v(577)=c0(2)*dcl0(3,1)
    v(576)=c0(3)*dcl0(2,1)
    v(575)=c0(4)*dcl0(3,3)+c0(3)*dcl0(4,3)-c0(6)*dcl0(5,3)-c0(5)*dcl0(6,3)
    v(574)=c0(4)*dcl0(3,2)+c0(3)*dcl0(4,2)-c0(6)*dcl0(5,2)-c0(5)*dcl0(6,2)
    v(573)=c0(4)*dcl0(3,1)+c0(3)*dcl0(4,1)-c0(6)*dcl0(5,1)-c0(5)*dcl0(6,1)
    v(572)=-(c0(5)*dcl0(2,3))+c0(6)*dcl0(4,3)-c0(2)*dcl0(5,3)+c0(4)*dcl0(6,3)
    v(571)=-(c0(5)*dcl0(2,2))+c0(6)*dcl0(4,2)-c0(2)*dcl0(5,2)+c0(4)*dcl0(6,2)
    v(570)=-(c0(5)*dcl0(2,1))+c0(6)*dcl0(4,1)-c0(2)*dcl0(5,1)+c0(4)*dcl0(6,1)
    v(568)=c1(2)*dcl1(3,3)
    v(567)=c1(3)*dcl1(2,3)
    v(566)=c1(2)*dcl1(3,2)
    v(565)=c1(3)*dcl1(2,2)
    v(564)=c1(2)*dcl1(3,1)
    v(563)=c1(3)*dcl1(2,1)
    v(562)=c1(4)*dcl1(3,3)+c1(3)*dcl1(4,3)-c1(6)*dcl1(5,3)-c1(5)*dcl1(6,3)
    v(561)=c1(4)*dcl1(3,2)+c1(3)*dcl1(4,2)-c1(6)*dcl1(5,2)-c1(5)*dcl1(6,2)
    v(560)=c1(4)*dcl1(3,1)+c1(3)*dcl1(4,1)-c1(6)*dcl1(5,1)-c1(5)*dcl1(6,1)
    v(559)=-(c1(5)*dcl1(2,3))+c1(6)*dcl1(4,3)-c1(2)*dcl1(5,3)+c1(4)*dcl1(6,3)
    v(558)=-(c1(5)*dcl1(2,2))+c1(6)*dcl1(4,2)-c1(2)*dcl1(5,2)+c1(4)*dcl1(6,2)
    v(557)=-(c1(5)*dcl1(2,1))+c1(6)*dcl1(4,1)-c1(2)*dcl1(5,1)+c1(4)*dcl1(6,1)
    v(555)=d2c0(6,3,3)
    v(554)=d2c0(6,3,2)
    v(553)=d2c0(6,3,1)
    v(552)=d2c0(6,2,3)
    v(551)=d2c0(6,2,2)
    v(550)=d2c0(6,2,1)
    v(549)=d2c0(6,1,3)
    v(548)=d2c0(6,1,2)
    v(547)=d2c0(6,1,1)
    v(546)=d2c0(5,3,3)
    v(545)=d2c0(5,3,2)
    v(544)=d2c0(5,3,1)
    v(543)=d2c0(5,2,3)
    v(542)=d2c0(5,2,2)
    v(541)=d2c0(5,2,1)
    v(540)=d2c0(5,1,3)
    v(539)=d2c0(5,1,2)
    v(538)=d2c0(5,1,1)
    v(537)=d2c0(4,3,3)
    v(536)=d2c0(4,3,2)
    v(535)=d2c0(4,3,1)
    v(534)=d2c0(4,2,3)
    v(533)=d2c0(4,2,2)
    v(532)=d2c0(4,2,1)
    v(531)=d2c0(4,1,3)
    v(530)=d2c0(4,1,2)
    v(529)=d2c0(4,1,1)
    v(528)=d2c0(3,3,3)
    v(527)=d2c0(3,3,2)
    v(526)=d2c0(3,3,1)
    v(525)=d2c0(3,2,3)
    v(524)=d2c0(3,2,2)
    v(523)=d2c0(3,2,1)
    v(522)=d2c0(3,1,3)
    v(521)=d2c0(3,1,2)
    v(520)=d2c0(3,1,1)
    v(519)=d2c0(2,3,3)
    v(518)=d2c0(2,3,2)
    v(517)=d2c0(2,3,1)
    v(516)=d2c0(2,2,3)
    v(515)=d2c0(2,2,2)
    v(514)=d2c0(2,2,1)
    v(513)=d2c0(2,1,3)
    v(512)=d2c0(2,1,2)
    v(511)=d2c0(2,1,1)
    v(510)=d2c1(6,3,3)
    v(509)=d2c1(6,3,2)
    v(508)=d2c1(6,3,1)
    v(507)=d2c1(6,2,3)
    v(506)=d2c1(6,2,2)
    v(505)=d2c1(6,2,1)
    v(504)=d2c1(6,1,3)
    v(503)=d2c1(6,1,2)
    v(502)=d2c1(6,1,1)
    v(501)=d2c1(5,3,3)
    v(500)=d2c1(5,3,2)
    v(499)=d2c1(5,3,1)
    v(498)=d2c1(5,2,3)
    v(497)=d2c1(5,2,2)
    v(496)=d2c1(5,2,1)
    v(495)=d2c1(5,1,3)
    v(494)=d2c1(5,1,2)
    v(493)=d2c1(5,1,1)
    v(492)=d2c1(4,3,3)
    v(491)=d2c1(4,3,2)
    v(490)=d2c1(4,3,1)
    v(489)=d2c1(4,2,3)
    v(488)=d2c1(4,2,2)
    v(487)=d2c1(4,2,1)
    v(486)=d2c1(4,1,3)
    v(485)=d2c1(4,1,2)
    v(484)=d2c1(4,1,1)
    v(483)=d2c1(3,3,3)
    v(482)=d2c1(3,3,2)
    v(481)=d2c1(3,3,1)
    v(480)=d2c1(3,2,3)
    v(479)=d2c1(3,2,2)
    v(478)=d2c1(3,2,1)
    v(477)=d2c1(3,1,3)
    v(476)=d2c1(3,1,2)
    v(475)=d2c1(3,1,1)
    v(474)=d2c1(2,3,3)
    v(473)=d2c1(2,3,2)
    v(472)=d2c1(2,3,1)
    v(471)=d2c1(2,2,3)
    v(470)=d2c1(2,2,2)
    v(469)=d2c1(2,2,1)
    v(468)=d2c1(2,1,3)
    v(467)=d2c1(2,1,2)
    v(466)=d2c1(2,1,1)
    v(465)=d2c1(1,3,3)
    v(464)=d2c1(1,3,2)
    v(463)=d2c1(1,3,1)
    v(462)=d2c1(1,2,3)
    v(461)=d2c1(1,2,2)
    v(460)=d2c1(1,2,1)
    v(459)=d2c1(1,1,3)
    v(458)=d2c1(1,1,2)
    v(457)=d2c1(1,1,1)
    v(456)=c0(2)*dck0(3,3)
    v(455)=c0(3)*dck0(2,3)
    v(454)=c0(4)*dck0(3,3)+c0(3)*dck0(4,3)-c0(6)*dck0(5,3)-c0(5)*dck0(6,3)
    v(453)=-(c0(5)*dck0(2,3))+c0(6)*dck0(4,3)-c0(2)*dck0(5,3)+c0(4)*dck0(6,3)
    v(451)=c0(2)*dck0(3,2)
    v(450)=c0(3)*dck0(2,2)
    v(449)=(-2d0)*c0(6)
    v(448)=c0(4)*dck0(3,2)+c0(3)*dck0(4,2)-c0(6)*dck0(5,2)-c0(5)*dck0(6,2)
    v(447)=-(c0(5)*dck0(2,2))+c0(6)*dck0(4,2)-c0(2)*dck0(5,2)+c0(4)*dck0(6,2)
    v(446)=c0(3)*dck0(2,1)+c0(2)*dck0(3,1)+dck0(6,1)*v(449)
    v(445)=c0(4)*dck0(3,1)+c0(3)*dck0(4,1)-c0(6)*dck0(5,1)-c0(5)*dck0(6,1)
    v(444)=-(c0(5)*dck0(2,1))+c0(6)*dck0(4,1)-c0(2)*dck0(5,1)+c0(4)*dck0(6,1)
    v(443)=c1(2)*dck1(3,3)
    v(442)=c1(3)*dck1(2,3)
    v(441)=c1(4)*dck1(3,3)+c1(3)*dck1(4,3)-c1(6)*dck1(5,3)-c1(5)*dck1(6,3)
    v(440)=-(c1(5)*dck1(2,3))+c1(6)*dck1(4,3)-c1(2)*dck1(5,3)+c1(4)*dck1(6,3)
    v(439)=c1(2)*dck1(3,2)
    v(438)=c1(3)*dck1(2,2)
    v(437)=(-2d0)*c1(6)
    v(436)=c1(4)*dck1(3,2)+c1(3)*dck1(4,2)-c1(6)*dck1(5,2)-c1(5)*dck1(6,2)
    v(435)=-(c1(5)*dck1(2,2))+c1(6)*dck1(4,2)-c1(2)*dck1(5,2)+c1(4)*dck1(6,2)
    v(434)=c1(3)*dck1(2,1)+c1(2)*dck1(3,1)+dck1(6,1)*v(437)
    v(433)=c1(4)*dck1(3,1)+c1(3)*dck1(4,1)-c1(6)*dck1(5,1)-c1(5)*dck1(6,1)
    v(432)=-(c1(5)*dck1(2,1))+c1(6)*dck1(4,1)-c1(2)*dck1(5,1)+c1(4)*dck1(6,1)
    v(431)=c0(2)*c0(3)-c0(6)**2
    v(430)=c0(3)*c0(4)-c0(5)*c0(6)
    v(429)=-(c0(2)*c0(5))+c0(4)*c0(6)
    v(428)=c1(2)*c1(3)-c1(6)**2
    v(427)=c1(3)*c1(4)-c1(5)*c1(6)
    v(426)=-(c1(2)*c1(5))+c1(4)*c1(6)
    v(261)=c1(5)*v(426)-c1(4)*v(427)+c1(1)*v(428)
    v(310)=1d0/v(261)**3
    v(569)=(-2d0)*v(310)
    v(271)=1d0/v(261)**2
    v(262)=c0(5)*v(429)-c0(4)*v(430)+c0(1)*v(431)
    v(452)=-(v(262)*v(271))
    v(260)=v(262)/v(261)
    v(329)=1d0/v(260)**0.16666666666666669d1
    v(648)=(-2d0/3d0)*v(329)
    v(270)=1d0/v(260)**0.6666666666666666d0
    v(556)=v(270)/3d0
    v(254)=v(260)**0.3333333333333333d0
    v(339)=dck1(5,1)*v(426)-dck1(4,1)*v(427)+dck1(1,1)*v(428)+c1(5)*v(432)-c1(4)*v(433)+c1(1)*v(434)
    v(347)=dck1(6,2)*v(437)+v(438)+v(439)
    v(351)=c1(1)*v(347)+dck1(5,2)*v(426)-dck1(4,2)*v(427)+dck1(1,2)*v(428)+c1(5)*v(435)-c1(4)*v(436)
    v(359)=dck1(6,3)*v(437)+v(442)+v(443)
    v(363)=c1(1)*v(359)+dck1(5,3)*v(426)-dck1(4,3)*v(427)+dck1(1,3)*v(428)+c1(5)*v(440)-c1(4)*v(441)
    v(338)=dck0(5,1)*v(429)-dck0(4,1)*v(430)+dck0(1,1)*v(431)+c0(5)*v(444)-c0(4)*v(445)+c0(1)*v(446)
    v(343)=v(338)/v(261)+v(339)*v(452)
    v(352)=dck0(6,2)*v(449)+v(450)+v(451)
    v(350)=c0(1)*v(352)+dck0(5,2)*v(429)-dck0(4,2)*v(430)+dck0(1,2)*v(431)+c0(5)*v(447)-c0(4)*v(448)
    v(355)=v(350)/v(261)+v(351)*v(452)
    v(364)=dck0(6,3)*v(449)+v(455)+v(456)
    v(362)=c0(1)*v(364)+dck0(5,3)*v(429)-dck0(4,3)*v(430)+dck0(1,3)*v(431)+c0(5)*v(453)-c0(4)*v(454)
    v(367)=v(362)/v(261)+v(363)*v(452)
    v(269)=v(343)*v(556)
    v(273)=v(355)*v(556)
    v(275)=v(367)*v(556)
    v(300)=dcl1(6,1)*v(437)+v(563)+v(564)
    v(301)=dcl1(6,2)*v(437)+v(565)+v(566)
    v(302)=dcl1(6,3)*v(437)+v(567)+v(568)
    v(303)=c1(1)*v(300)+dcl1(5,1)*v(426)-dcl1(4,1)*v(427)+dcl1(1,1)*v(428)+c1(5)*v(557)-c1(4)*v(560)
    v(304)=c1(1)*v(301)+dcl1(5,2)*v(426)-dcl1(4,2)*v(427)+dcl1(1,2)*v(428)+c1(5)*v(558)-c1(4)*v(561)
    v(305)=c1(1)*v(302)+dcl1(5,3)*v(426)-dcl1(4,3)*v(427)+dcl1(1,3)*v(428)+c1(5)*v(559)-c1(4)*v(562)
    v(306)=-(v(271)*v(303))
    v(307)=-(v(271)*v(304))
    v(308)=-(v(271)*v(305))
    v(309)=v(303)*v(569)
    v(604)=-(v(262)*v(309))
    v(311)=v(304)*v(569)
    v(671)=-(v(262)*v(311))
    v(312)=v(305)*v(569)
    v(737)=-(v(262)*v(312))
    v(319)=dcl0(6,1)*v(449)+v(576)+v(577)
    v(320)=dcl0(6,2)*v(449)+v(578)+v(579)
    v(321)=dcl0(6,3)*v(449)+v(580)+v(581)
    v(322)=c0(1)*v(319)+dcl0(5,1)*v(429)-dcl0(4,1)*v(430)+dcl0(1,1)*v(431)+c0(5)*v(570)-c0(4)*v(573)
    v(605)=-(v(271)*v(322))
    v(323)=c0(1)*v(320)+dcl0(5,2)*v(429)-dcl0(4,2)*v(430)+dcl0(1,2)*v(431)+c0(5)*v(571)-c0(4)*v(574)
    v(672)=-(v(271)*v(323))
    v(324)=c0(1)*v(321)+dcl0(5,3)*v(429)-dcl0(4,3)*v(430)+dcl0(1,3)*v(431)+c0(5)*v(572)-c0(4)*v(575)
    v(738)=-(v(271)*v(324))
    v(325)=v(262)*v(306)+v(322)/v(261)
    v(326)=v(262)*v(307)+v(323)/v(261)
    v(327)=v(262)*v(308)+v(324)/v(261)
    v(328)=v(325)*v(648)
    v(368)=(v(328)*v(367)+v(270)*(v(306)*v(362)+v(452)*(dck1(1,3)*v(300)+dcl1(1,1)*v(359)+dcl1(5,1)*v(440)-dcl1(4,1)*v(441)&
         &+v(428)*v(463)-v(427)*v(490)+v(426)*v(499)+dck1(5,3)*v(557)-dck1(4,3)*v(560)+c1(5)*(-(c1(5)*v(472))+c1(6)*v(490)-c1(2&
         &)*v(499)+c1(4)*v(508)+v(582)+v(583)+v(584)+v(585))-c1(4)*(c1(4)*v(481)+c1(3)*v(490)-c1(6)*v(499)-c1(5)*v(508)+v(586)+v&
         &(587)+v(588)+v(589))+c1(1)*(c1(3)*v(472)+c1(2)*v(481)+v(437)*v(508)+v(590)+v(591)+v(592)))+(dck0(1,3)*v(319)+dcl0(1,1&
         &)*v(364)+d2c0(1,3,1)*v(431)+dcl0(5,1)*v(453)-dcl0(4,1)*v(454)-v(430)*v(535)+v(429)*v(544)+dck0(5,3)*v(570)-dck0(4,3)*v&
         &(573)+c0(5)*(-(c0(5)*v(517))+c0(6)*v(535)-c0(2)*v(544)+c0(4)*v(553)+v(593)+v(594)+v(595)+v(596))-c0(4)*(c0(4)*v(526)+c0&
         &(3)*v(535)-c0(6)*v(544)-c0(5)*v(553)+v(597)+v(598)+v(599)+v(600))+c0(1)*(c0(3)*v(517)+c0(2)*v(526)+v(449)*v(553)+v(601)&
         &+v(602)+v(603)))/v(261)+v(363)*v(604)+v(363)*v(605)))/3d0
    v(356)=(v(328)*v(355)+v(270)*(v(306)*v(350)+v(351)*v(604)+v(351)*v(605)+v(452)*(dck1(1,2)*v(300)+dcl1(1,1)*v(347)+dcl1&
         &(5,1)*v(435)-dcl1(4,1)*v(436)+v(428)*v(460)-v(427)*v(487)+v(426)*v(496)+dck1(5,2)*v(557)-dck1(4,2)*v(560)+c1(5)*(-(c1(5&
         &)*v(469))+c1(6)*v(487)-c1(2)*v(496)+c1(4)*v(505)+v(608)+v(609)+v(610)+v(611))-c1(4)*(c1(4)*v(478)+c1(3)*v(487)-c1(6)*v&
         &(496)-c1(5)*v(505)+v(612)+v(613)+v(614)+v(615))+c1(1)*(c1(3)*v(469)+c1(2)*v(478)+v(437)*v(505)+dck1(6,2)*v(606)+v(616)&
         &+v(617)))+(dck0(1,2)*v(319)+dcl0(1,1)*v(352)+d2c0(1,2,1)*v(431)+dcl0(5,1)*v(447)-dcl0(4,1)*v(448)-v(430)*v(532)+v(429&
         &)*v(541)+dck0(5,2)*v(570)-dck0(4,2)*v(573)+c0(5)*(-(c0(5)*v(514))+c0(6)*v(532)-c0(2)*v(541)+c0(4)*v(550)+v(618)+v(619)&
         &+v(620)+v(621))-c0(4)*(c0(4)*v(523)+c0(3)*v(532)-c0(6)*v(541)-c0(5)*v(550)+v(622)+v(623)+v(624)+v(625))+c0(1)*(c0(3)*v&
         &(514)+c0(2)*v(523)+v(449)*v(550)+dck0(6,2)*v(607)+v(626)+v(627)))/v(261)))/3d0
    v(344)=(v(328)*v(343)+v(270)*(v(306)*v(338)+v(339)*v(604)+v(339)*v(605)+v(452)*(dck1(1,1)*v(300)+dcl1(5,1)*v(432)-dcl1&
         &(4,1)*v(433)+dcl1(1,1)*v(434)+v(428)*v(457)-v(427)*v(484)+v(426)*v(493)+dck1(5,1)*v(557)-dck1(4,1)*v(560)+c1(5)*(-(c1(5&
         &)*v(466))+c1(6)*v(484)-c1(2)*v(493)+c1(4)*v(502)+v(628)+v(629)+v(630)+v(631))-c1(4)*(c1(4)*v(475)+c1(3)*v(484)-c1(6)*v&
         &(493)-c1(5)*v(502)+v(632)+v(633)+v(634)+v(635))+c1(1)*(c1(3)*v(466)+c1(2)*v(475)+v(437)*v(502)+dck1(6,1)*v(606)+v(636)&
         &+v(637)))+(dck0(1,1)*v(319)+d2c0(1,1,1)*v(431)+dcl0(5,1)*v(444)-dcl0(4,1)*v(445)+dcl0(1,1)*v(446)-v(430)*v(529)+v(429&
         &)*v(538)+dck0(5,1)*v(570)-dck0(4,1)*v(573)+c0(5)*(-(c0(5)*v(511))+c0(6)*v(529)-c0(2)*v(538)+c0(4)*v(547)+v(638)+v(639)&
         &+v(640)+v(641))-c0(4)*(c0(4)*v(520)+c0(3)*v(529)-c0(6)*v(538)-c0(5)*v(547)+v(642)+v(643)+v(644)+v(645))+c0(1)*(c0(3)*v&
         &(511)+c0(2)*v(520)+v(449)*v(547)+dck0(6,1)*v(607)+v(646)+v(647)))/v(261)))/3d0
    v(330)=v(326)*v(648)
    v(369)=(v(330)*v(367)+v(270)*(v(307)*v(362)+v(452)*(dck1(1,3)*v(301)+dcl1(1,2)*v(359)+dcl1(5,2)*v(440)-dcl1(4,2)*v(441)&
         &+v(428)*v(464)-v(427)*v(491)+v(426)*v(500)+dck1(5,3)*v(558)-dck1(4,3)*v(561)+c1(5)*(-(c1(5)*v(473))+c1(6)*v(491)-c1(2&
         &)*v(500)+c1(4)*v(509)+v(649)+v(650)+v(651)+v(652))-c1(4)*(c1(4)*v(482)+c1(3)*v(491)-c1(6)*v(500)-c1(5)*v(509)+v(653)+v&
         &(654)+v(655)+v(656))+c1(1)*(c1(3)*v(473)+c1(2)*v(482)+v(437)*v(509)+v(657)+v(658)+v(659)))+(dck0(1,3)*v(320)+dcl0(1,2&
         &)*v(364)+d2c0(1,3,2)*v(431)+dcl0(5,2)*v(453)-dcl0(4,2)*v(454)-v(430)*v(536)+v(429)*v(545)+dck0(5,3)*v(571)-dck0(4,3)*v&
         &(574)+c0(5)*(-(c0(5)*v(518))+c0(6)*v(536)-c0(2)*v(545)+c0(4)*v(554)+v(660)+v(661)+v(662)+v(663))-c0(4)*(c0(4)*v(527)+c0&
         &(3)*v(536)-c0(6)*v(545)-c0(5)*v(554)+v(664)+v(665)+v(666)+v(667))+c0(1)*(c0(3)*v(518)+c0(2)*v(527)+v(449)*v(554)+v(668)&
         &+v(669)+v(670)))/v(261)+v(363)*v(671)+v(363)*v(672)))/3d0
    v(357)=(v(330)*v(355)+v(270)*(v(307)*v(350)+v(351)*v(671)+v(351)*v(672)+v(452)*(dck1(1,2)*v(301)+dcl1(1,2)*v(347)+dcl1&
         &(5,2)*v(435)-dcl1(4,2)*v(436)+v(428)*v(461)-v(427)*v(488)+v(426)*v(497)+dck1(5,2)*v(558)-dck1(4,2)*v(561)+c1(5)*(-(c1(5&
         &)*v(470))+c1(6)*v(488)-c1(2)*v(497)+c1(4)*v(506)+v(675)+v(676)+v(677)+v(678))-c1(4)*(c1(4)*v(479)+c1(3)*v(488)-c1(6)*v&
         &(497)-c1(5)*v(506)+v(679)+v(680)+v(681)+v(682))+c1(1)*(c1(3)*v(470)+c1(2)*v(479)+v(437)*v(506)+dck1(6,2)*v(673)+v(683)&
         &+v(684)))+(dck0(1,2)*v(320)+dcl0(1,2)*v(352)+d2c0(1,2,2)*v(431)+dcl0(5,2)*v(447)-dcl0(4,2)*v(448)-v(430)*v(533)+v(429&
         &)*v(542)+dck0(5,2)*v(571)-dck0(4,2)*v(574)+c0(5)*(-(c0(5)*v(515))+c0(6)*v(533)-c0(2)*v(542)+c0(4)*v(551)+v(685)+v(686)&
         &+v(687)+v(688))-c0(4)*(c0(4)*v(524)+c0(3)*v(533)-c0(6)*v(542)-c0(5)*v(551)+v(689)+v(690)+v(691)+v(692))+c0(1)*(c0(3)*v&
         &(515)+c0(2)*v(524)+v(449)*v(551)+dck0(6,2)*v(674)+v(693)+v(694)))/v(261)))/3d0
    v(345)=(v(330)*v(343)+v(270)*(v(307)*v(338)+v(339)*v(671)+v(339)*v(672)+v(452)*(dck1(1,1)*v(301)+dcl1(5,2)*v(432)-dcl1&
         &(4,2)*v(433)+dcl1(1,2)*v(434)+v(428)*v(458)-v(427)*v(485)+v(426)*v(494)+dck1(5,1)*v(558)-dck1(4,1)*v(561)+c1(5)*(-(c1(5&
         &)*v(467))+c1(6)*v(485)-c1(2)*v(494)+c1(4)*v(503)+v(695)+v(696)+v(697)+v(698))-c1(4)*(c1(4)*v(476)+c1(3)*v(485)-c1(6)*v&
         &(494)-c1(5)*v(503)+v(699)+v(700)+v(701)+v(702))+c1(1)*(c1(3)*v(467)+c1(2)*v(476)+v(437)*v(503)+dck1(6,1)*v(673)+v(703)&
         &+v(704)))+(dck0(1,1)*v(320)+d2c0(1,1,2)*v(431)+dcl0(5,2)*v(444)-dcl0(4,2)*v(445)+dcl0(1,2)*v(446)-v(430)*v(530)+v(429&
         &)*v(539)+dck0(5,1)*v(571)-dck0(4,1)*v(574)+c0(5)*(-(c0(5)*v(512))+c0(6)*v(530)-c0(2)*v(539)+c0(4)*v(548)+v(705)+v(706)&
         &+v(707)+v(708))-c0(4)*(c0(4)*v(521)+c0(3)*v(530)-c0(6)*v(539)-c0(5)*v(548)+v(709)+v(710)+v(711)+v(712))+c0(1)*(c0(3)*v&
         &(512)+c0(2)*v(521)+v(449)*v(548)+dck0(6,1)*v(674)+v(713)+v(714)))/v(261)))/3d0
    v(331)=v(327)*v(648)
    v(370)=(v(331)*v(367)+v(270)*(v(308)*v(362)+v(452)*(dck1(1,3)*v(302)+dcl1(1,3)*v(359)+dcl1(5,3)*v(440)-dcl1(4,3)*v(441)&
         &+v(428)*v(465)-v(427)*v(492)+v(426)*v(501)+dck1(5,3)*v(559)-dck1(4,3)*v(562)+c1(5)*(-(c1(5)*v(474))+c1(6)*v(492)-c1(2&
         &)*v(501)+c1(4)*v(510)+v(715)+v(716)+v(717)+v(718))-c1(4)*(c1(4)*v(483)+c1(3)*v(492)-c1(6)*v(501)-c1(5)*v(510)+v(719)+v&
         &(720)+v(721)+v(722))+c1(1)*(c1(3)*v(474)+c1(2)*v(483)+v(437)*v(510)+v(723)+v(724)+v(725)))+(dck0(1,3)*v(321)+dcl0(1,3&
         &)*v(364)+d2c0(1,3,3)*v(431)+dcl0(5,3)*v(453)-dcl0(4,3)*v(454)-v(430)*v(537)+v(429)*v(546)+dck0(5,3)*v(572)-dck0(4,3)*v&
         &(575)+c0(5)*(-(c0(5)*v(519))+c0(6)*v(537)-c0(2)*v(546)+c0(4)*v(555)+v(726)+v(727)+v(728)+v(729))-c0(4)*(c0(4)*v(528)+c0&
         &(3)*v(537)-c0(6)*v(546)-c0(5)*v(555)+v(730)+v(731)+v(732)+v(733))+c0(1)*(c0(3)*v(519)+c0(2)*v(528)+v(449)*v(555)+v(734)&
         &+v(735)+v(736)))/v(261)+v(363)*v(737)+v(363)*v(738)))/3d0
    v(358)=(v(331)*v(355)+v(270)*(v(308)*v(350)+v(351)*v(737)+v(351)*v(738)+v(452)*(dck1(1,2)*v(302)+dcl1(1,3)*v(347)+dcl1&
         &(5,3)*v(435)-dcl1(4,3)*v(436)+v(428)*v(462)-v(427)*v(489)+v(426)*v(498)+dck1(5,2)*v(559)-dck1(4,2)*v(562)+c1(5)*(-(c1(5&
         &)*v(471))+c1(6)*v(489)-c1(2)*v(498)+c1(4)*v(507)+v(741)+v(742)+v(743)+v(744))-c1(4)*(c1(4)*v(480)+c1(3)*v(489)-c1(6)*v&
         &(498)-c1(5)*v(507)+v(745)+v(746)+v(747)+v(748))+c1(1)*(c1(3)*v(471)+c1(2)*v(480)+v(437)*v(507)+dck1(6,2)*v(739)+v(749)&
         &+v(750)))+(dck0(1,2)*v(321)+dcl0(1,3)*v(352)+d2c0(1,2,3)*v(431)+dcl0(5,3)*v(447)-dcl0(4,3)*v(448)-v(430)*v(534)+v(429&
         &)*v(543)+dck0(5,2)*v(572)-dck0(4,2)*v(575)+c0(5)*(-(c0(5)*v(516))+c0(6)*v(534)-c0(2)*v(543)+c0(4)*v(552)+v(751)+v(752)&
         &+v(753)+v(754))-c0(4)*(c0(4)*v(525)+c0(3)*v(534)-c0(6)*v(543)-c0(5)*v(552)+v(755)+v(756)+v(757)+v(758))+c0(1)*(c0(3)*v&
         &(516)+c0(2)*v(525)+v(449)*v(552)+dck0(6,2)*v(740)+v(759)+v(760)))/v(261)))/3d0
    v(346)=(v(331)*v(343)+v(270)*(v(308)*v(338)+v(339)*v(737)+v(339)*v(738)+v(452)*(dck1(1,1)*v(302)+dcl1(5,3)*v(432)-dcl1&
         &(4,3)*v(433)+dcl1(1,3)*v(434)+v(428)*v(459)-v(427)*v(486)+v(426)*v(495)+dck1(5,1)*v(559)-dck1(4,1)*v(562)+c1(5)*(-(c1(5&
         &)*v(468))+c1(6)*v(486)-c1(2)*v(495)+c1(4)*v(504)+v(761)+v(762)+v(763)+v(764))-c1(4)*(c1(4)*v(477)+c1(3)*v(486)-c1(6)*v&
         &(495)-c1(5)*v(504)+v(765)+v(766)+v(767)+v(768))+c1(1)*(c1(3)*v(468)+c1(2)*v(477)+v(437)*v(504)+dck1(6,1)*v(739)+v(769)&
         &+v(770)))+(dck0(1,1)*v(321)+d2c0(1,1,3)*v(431)+dcl0(5,3)*v(444)-dcl0(4,3)*v(445)+dcl0(1,3)*v(446)-v(430)*v(531)+v(429&
         &)*v(540)+dck0(5,1)*v(572)-dck0(4,1)*v(575)+c0(5)*(-(c0(5)*v(513))+c0(6)*v(531)-c0(2)*v(540)+c0(4)*v(549)+v(771)+v(772)&
         &+v(773)+v(774))-c0(4)*(c0(4)*v(522)+c0(3)*v(531)-c0(6)*v(540)-c0(5)*v(549)+v(775)+v(776)+v(777)+v(778))+c0(1)*(c0(3)*v&
         &(513)+c0(2)*v(522)+v(449)*v(549)+dck0(6,1)*v(740)+v(779)+v(780)))/v(261)))/3d0
    v(332)=v(325)*v(556)
    v(333)=v(326)*v(556)
    v(334)=v(327)*v(556)
    d2c(1,1,1)=dcl1(1,1)*v(269)+dck1(1,1)*v(332)+c1(1)*v(344)+v(254)*v(457)
    d2c(1,1,2)=dcl1(1,2)*v(269)+dck1(1,1)*v(333)+c1(1)*v(345)+v(254)*v(458)
    d2c(1,1,3)=dcl1(1,3)*v(269)+dck1(1,1)*v(334)+c1(1)*v(346)+v(254)*v(459)
    d2c(1,2,1)=dcl1(1,1)*v(273)+dck1(1,2)*v(332)+c1(1)*v(356)+v(254)*v(460)
    d2c(1,2,2)=dcl1(1,2)*v(273)+dck1(1,2)*v(333)+c1(1)*v(357)+v(254)*v(461)
    d2c(1,2,3)=dcl1(1,3)*v(273)+dck1(1,2)*v(334)+c1(1)*v(358)+v(254)*v(462)
    d2c(1,3,1)=dcl1(1,1)*v(275)+dck1(1,3)*v(332)+c1(1)*v(368)+v(254)*v(463)
    d2c(1,3,2)=dcl1(1,2)*v(275)+dck1(1,3)*v(333)+c1(1)*v(369)+v(254)*v(464)
    d2c(1,3,3)=dcl1(1,3)*v(275)+dck1(1,3)*v(334)+c1(1)*v(370)+v(254)*v(465)
    d2c(2,1,1)=dcl1(2,1)*v(269)+dck1(2,1)*v(332)+c1(2)*v(344)+v(254)*v(466)
    d2c(2,1,2)=dcl1(2,2)*v(269)+dck1(2,1)*v(333)+c1(2)*v(345)+v(254)*v(467)
    d2c(2,1,3)=dcl1(2,3)*v(269)+dck1(2,1)*v(334)+c1(2)*v(346)+v(254)*v(468)
    d2c(2,2,1)=dcl1(2,1)*v(273)+dck1(2,2)*v(332)+c1(2)*v(356)+v(254)*v(469)
    d2c(2,2,2)=dcl1(2,2)*v(273)+dck1(2,2)*v(333)+c1(2)*v(357)+v(254)*v(470)
    d2c(2,2,3)=dcl1(2,3)*v(273)+dck1(2,2)*v(334)+c1(2)*v(358)+v(254)*v(471)
    d2c(2,3,1)=dcl1(2,1)*v(275)+dck1(2,3)*v(332)+c1(2)*v(368)+v(254)*v(472)
    d2c(2,3,2)=dcl1(2,2)*v(275)+dck1(2,3)*v(333)+c1(2)*v(369)+v(254)*v(473)
    d2c(2,3,3)=dcl1(2,3)*v(275)+dck1(2,3)*v(334)+c1(2)*v(370)+v(254)*v(474)
    d2c(3,1,1)=dcl1(3,1)*v(269)+dck1(3,1)*v(332)+c1(3)*v(344)+v(254)*v(475)
    d2c(3,1,2)=dcl1(3,2)*v(269)+dck1(3,1)*v(333)+c1(3)*v(345)+v(254)*v(476)
    d2c(3,1,3)=dcl1(3,3)*v(269)+dck1(3,1)*v(334)+c1(3)*v(346)+v(254)*v(477)
    d2c(3,2,1)=dcl1(3,1)*v(273)+dck1(3,2)*v(332)+c1(3)*v(356)+v(254)*v(478)
    d2c(3,2,2)=dcl1(3,2)*v(273)+dck1(3,2)*v(333)+c1(3)*v(357)+v(254)*v(479)
    d2c(3,2,3)=dcl1(3,3)*v(273)+dck1(3,2)*v(334)+c1(3)*v(358)+v(254)*v(480)
    d2c(3,3,1)=dcl1(3,1)*v(275)+dck1(3,3)*v(332)+c1(3)*v(368)+v(254)*v(481)
    d2c(3,3,2)=dcl1(3,2)*v(275)+dck1(3,3)*v(333)+c1(3)*v(369)+v(254)*v(482)
    d2c(3,3,3)=dcl1(3,3)*v(275)+dck1(3,3)*v(334)+c1(3)*v(370)+v(254)*v(483)
    d2c(4,1,1)=dcl1(4,1)*v(269)+dck1(4,1)*v(332)+c1(4)*v(344)+v(254)*v(484)
    d2c(4,1,2)=dcl1(4,2)*v(269)+dck1(4,1)*v(333)+c1(4)*v(345)+v(254)*v(485)
    d2c(4,1,3)=dcl1(4,3)*v(269)+dck1(4,1)*v(334)+c1(4)*v(346)+v(254)*v(486)
    d2c(4,2,1)=dcl1(4,1)*v(273)+dck1(4,2)*v(332)+c1(4)*v(356)+v(254)*v(487)
    d2c(4,2,2)=dcl1(4,2)*v(273)+dck1(4,2)*v(333)+c1(4)*v(357)+v(254)*v(488)
    d2c(4,2,3)=dcl1(4,3)*v(273)+dck1(4,2)*v(334)+c1(4)*v(358)+v(254)*v(489)
    d2c(4,3,1)=dcl1(4,1)*v(275)+dck1(4,3)*v(332)+c1(4)*v(368)+v(254)*v(490)
    d2c(4,3,2)=dcl1(4,2)*v(275)+dck1(4,3)*v(333)+c1(4)*v(369)+v(254)*v(491)
    d2c(4,3,3)=dcl1(4,3)*v(275)+dck1(4,3)*v(334)+c1(4)*v(370)+v(254)*v(492)
    d2c(5,1,1)=dcl1(5,1)*v(269)+dck1(5,1)*v(332)+c1(5)*v(344)+v(254)*v(493)
    d2c(5,1,2)=dcl1(5,2)*v(269)+dck1(5,1)*v(333)+c1(5)*v(345)+v(254)*v(494)
    d2c(5,1,3)=dcl1(5,3)*v(269)+dck1(5,1)*v(334)+c1(5)*v(346)+v(254)*v(495)
    d2c(5,2,1)=dcl1(5,1)*v(273)+dck1(5,2)*v(332)+c1(5)*v(356)+v(254)*v(496)
    d2c(5,2,2)=dcl1(5,2)*v(273)+dck1(5,2)*v(333)+c1(5)*v(357)+v(254)*v(497)
    d2c(5,2,3)=dcl1(5,3)*v(273)+dck1(5,2)*v(334)+c1(5)*v(358)+v(254)*v(498)
    d2c(5,3,1)=dcl1(5,1)*v(275)+dck1(5,3)*v(332)+c1(5)*v(368)+v(254)*v(499)
    d2c(5,3,2)=dcl1(5,2)*v(275)+dck1(5,3)*v(333)+c1(5)*v(369)+v(254)*v(500)
    d2c(5,3,3)=dcl1(5,3)*v(275)+dck1(5,3)*v(334)+c1(5)*v(370)+v(254)*v(501)
    d2c(6,1,1)=dcl1(6,1)*v(269)+dck1(6,1)*v(332)+c1(6)*v(344)+v(254)*v(502)
    d2c(6,1,2)=dcl1(6,2)*v(269)+dck1(6,1)*v(333)+c1(6)*v(345)+v(254)*v(503)
    d2c(6,1,3)=dcl1(6,3)*v(269)+dck1(6,1)*v(334)+c1(6)*v(346)+v(254)*v(504)
    d2c(6,2,1)=dcl1(6,1)*v(273)+dck1(6,2)*v(332)+c1(6)*v(356)+v(254)*v(505)
    d2c(6,2,2)=dcl1(6,2)*v(273)+dck1(6,2)*v(333)+c1(6)*v(357)+v(254)*v(506)
    d2c(6,2,3)=dcl1(6,3)*v(273)+dck1(6,2)*v(334)+c1(6)*v(358)+v(254)*v(507)
    d2c(6,3,1)=dcl1(6,1)*v(275)+dck1(6,3)*v(332)+c1(6)*v(368)+v(254)*v(508)
    d2c(6,3,2)=dcl1(6,2)*v(275)+dck1(6,3)*v(333)+c1(6)*v(369)+v(254)*v(509)
    d2c(6,3,3)=dcl1(6,3)*v(275)+dck1(6,3)*v(334)+c1(6)*v(370)+v(254)*v(510)
  END SUBROUTINE mls_combinedfbar2


!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           8 Feb 22 12:29:31  *
!**************************************************************
! User     : Full professional version
! Notebook : mls_combinedfbar2
! Evaluation time                 : 19 s    Mode  : Optimal
! Number of formulae              : 401     Method: Automatic
! Subroutine                      : mls_combinedfbar2 size: 12767
! Total size of Mathematica  code : 12767 subexpressions
! Total size of Fortran code      : 26108 bytes

!******************* S U B R O U T I N E **********************
  SUBROUTINE mls_combinedfbar2old(c,cb,dcleft,dcright,dcbleft,dcbright,dc2,dcb2,dcs2)
    IMPLICIT NONE
    DOUBLE PRECISION v(770),c(6),cb(6),dcleft(6,3),dcright(6,3),dcbleft(6,3),dcbright(6,3),dc2(6,3,3),dcb2(6,3,3)&
         &,dcs2(6,3,3)
    v(765)=dcbleft(2,1)*dcbright(3,3)
    v(764)=dcbleft(3,1)*dcbright(2,3)
    v(763)=-(dcbleft(5,1)*dcbright(6,3))
    v(762)=-(dcbleft(6,1)*dcbright(5,3))
    v(761)=dcbleft(3,1)*dcbright(4,3)
    v(760)=dcbleft(4,1)*dcbright(3,3)
    v(759)=dcbleft(4,1)*dcbright(6,3)
    v(758)=-(dcbleft(2,1)*dcbright(5,3))
    v(757)=dcbleft(6,1)*dcbright(4,3)
    v(756)=-(dcbleft(5,1)*dcbright(2,3))
    v(755)=dcleft(2,1)*dcright(3,3)
    v(754)=dcleft(3,1)*dcright(2,3)
    v(753)=-(dcleft(5,1)*dcright(6,3))
    v(752)=-(dcleft(6,1)*dcright(5,3))
    v(751)=dcleft(3,1)*dcright(4,3)
    v(750)=dcleft(4,1)*dcright(3,3)
    v(749)=dcleft(4,1)*dcright(6,3)
    v(748)=-(dcleft(2,1)*dcright(5,3))
    v(747)=dcleft(6,1)*dcright(4,3)
    v(746)=-(dcleft(5,1)*dcright(2,3))
    v(745)=dcbleft(2,2)*dcbright(3,3)
    v(744)=dcbleft(3,2)*dcbright(2,3)
    v(743)=-(dcbleft(5,2)*dcbright(6,3))
    v(742)=-(dcbleft(6,2)*dcbright(5,3))
    v(741)=dcbleft(3,2)*dcbright(4,3)
    v(740)=dcbleft(4,2)*dcbright(3,3)
    v(739)=dcbleft(4,2)*dcbright(6,3)
    v(738)=-(dcbleft(2,2)*dcbright(5,3))
    v(737)=dcbleft(6,2)*dcbright(4,3)
    v(736)=-(dcbleft(5,2)*dcbright(2,3))
    v(735)=dcleft(2,2)*dcright(3,3)
    v(734)=dcleft(3,2)*dcright(2,3)
    v(733)=-(dcleft(5,2)*dcright(6,3))
    v(732)=-(dcleft(6,2)*dcright(5,3))
    v(731)=dcleft(3,2)*dcright(4,3)
    v(730)=dcleft(4,2)*dcright(3,3)
    v(729)=dcleft(4,2)*dcright(6,3)
    v(728)=-(dcleft(2,2)*dcright(5,3))
    v(727)=dcleft(6,2)*dcright(4,3)
    v(726)=-(dcleft(5,2)*dcright(2,3))
    v(725)=(-2d0)*dcbright(6,3)
    v(724)=(-2d0)*dcright(6,3)
    v(722)=dcbleft(6,3)*v(725)
    v(721)=dcbleft(2,3)*dcbright(3,3)
    v(720)=dcbleft(3,3)*dcbright(2,3)
    v(719)=-(dcbleft(5,3)*dcbright(6,3))
    v(718)=-(dcbleft(6,3)*dcbright(5,3))
    v(717)=dcbleft(3,3)*dcbright(4,3)
    v(716)=dcbleft(4,3)*dcbright(3,3)
    v(715)=dcbleft(4,3)*dcbright(6,3)
    v(714)=-(dcbleft(2,3)*dcbright(5,3))
    v(713)=dcbleft(6,3)*dcbright(4,3)
    v(712)=-(dcbleft(5,3)*dcbright(2,3))
    v(711)=dcleft(6,3)*v(724)
    v(710)=dcleft(2,3)*dcright(3,3)
    v(709)=dcleft(3,3)*dcright(2,3)
    v(708)=-(dcleft(5,3)*dcright(6,3))
    v(707)=-(dcleft(6,3)*dcright(5,3))
    v(706)=dcleft(3,3)*dcright(4,3)
    v(705)=dcleft(4,3)*dcright(3,3)
    v(704)=dcleft(4,3)*dcright(6,3)
    v(703)=-(dcleft(2,3)*dcright(5,3))
    v(702)=dcleft(6,3)*dcright(4,3)
    v(701)=-(dcleft(5,3)*dcright(2,3))
    v(700)=dcbleft(2,1)*dcbright(3,2)
    v(699)=dcbleft(3,1)*dcbright(2,2)
    v(698)=-(dcbleft(5,1)*dcbright(6,2))
    v(697)=-(dcbleft(6,1)*dcbright(5,2))
    v(696)=dcbleft(3,1)*dcbright(4,2)
    v(695)=dcbleft(4,1)*dcbright(3,2)
    v(694)=dcbleft(4,1)*dcbright(6,2)
    v(693)=-(dcbleft(2,1)*dcbright(5,2))
    v(692)=dcbleft(6,1)*dcbright(4,2)
    v(691)=-(dcbleft(5,1)*dcbright(2,2))
    v(690)=dcleft(2,1)*dcright(3,2)
    v(689)=dcleft(3,1)*dcright(2,2)
    v(688)=-(dcleft(5,1)*dcright(6,2))
    v(687)=-(dcleft(6,1)*dcright(5,2))
    v(686)=dcleft(3,1)*dcright(4,2)
    v(685)=dcleft(4,1)*dcright(3,2)
    v(684)=dcleft(4,1)*dcright(6,2)
    v(683)=-(dcleft(2,1)*dcright(5,2))
    v(682)=dcleft(6,1)*dcright(4,2)
    v(681)=-(dcleft(5,1)*dcright(2,2))
    v(680)=dcbleft(2,2)*dcbright(3,2)
    v(679)=dcbleft(3,2)*dcbright(2,2)
    v(678)=-(dcbleft(5,2)*dcbright(6,2))
    v(677)=-(dcbleft(6,2)*dcbright(5,2))
    v(676)=dcbleft(3,2)*dcbright(4,2)
    v(675)=dcbleft(4,2)*dcbright(3,2)
    v(674)=dcbleft(4,2)*dcbright(6,2)
    v(673)=-(dcbleft(2,2)*dcbright(5,2))
    v(672)=dcbleft(6,2)*dcbright(4,2)
    v(671)=-(dcbleft(5,2)*dcbright(2,2))
    v(670)=dcleft(2,2)*dcright(3,2)
    v(669)=dcleft(3,2)*dcright(2,2)
    v(668)=-(dcleft(5,2)*dcright(6,2))
    v(667)=-(dcleft(6,2)*dcright(5,2))
    v(666)=dcleft(3,2)*dcright(4,2)
    v(665)=dcleft(4,2)*dcright(3,2)
    v(664)=dcleft(4,2)*dcright(6,2)
    v(663)=-(dcleft(2,2)*dcright(5,2))
    v(662)=dcleft(6,2)*dcright(4,2)
    v(661)=-(dcleft(5,2)*dcright(2,2))
    v(660)=(-2d0)*dcbright(6,2)
    v(659)=(-2d0)*dcright(6,2)
    v(657)=dcbleft(6,3)*v(660)
    v(656)=dcbleft(2,3)*dcbright(3,2)
    v(655)=dcbleft(3,3)*dcbright(2,2)
    v(654)=-(dcbleft(5,3)*dcbright(6,2))
    v(653)=-(dcbleft(6,3)*dcbright(5,2))
    v(652)=dcbleft(3,3)*dcbright(4,2)
    v(651)=dcbleft(4,3)*dcbright(3,2)
    v(650)=dcbleft(4,3)*dcbright(6,2)
    v(649)=-(dcbleft(2,3)*dcbright(5,2))
    v(648)=dcbleft(6,3)*dcbright(4,2)
    v(647)=-(dcbleft(5,3)*dcbright(2,2))
    v(646)=dcleft(6,3)*v(659)
    v(645)=dcleft(2,3)*dcright(3,2)
    v(644)=dcleft(3,3)*dcright(2,2)
    v(643)=-(dcleft(5,3)*dcright(6,2))
    v(642)=-(dcleft(6,3)*dcright(5,2))
    v(641)=dcleft(3,3)*dcright(4,2)
    v(640)=dcleft(4,3)*dcright(3,2)
    v(639)=dcleft(4,3)*dcright(6,2)
    v(638)=-(dcleft(2,3)*dcright(5,2))
    v(637)=dcleft(6,3)*dcright(4,2)
    v(636)=-(dcleft(5,3)*dcright(2,2))
    v(634)=dcbleft(2,1)*dcbright(3,1)
    v(633)=dcbleft(3,1)*dcbright(2,1)
    v(632)=-(dcbleft(5,1)*dcbright(6,1))
    v(631)=-(dcbleft(6,1)*dcbright(5,1))
    v(630)=dcbleft(3,1)*dcbright(4,1)
    v(629)=dcbleft(4,1)*dcbright(3,1)
    v(628)=dcbleft(4,1)*dcbright(6,1)
    v(627)=-(dcbleft(2,1)*dcbright(5,1))
    v(626)=dcbleft(6,1)*dcbright(4,1)
    v(625)=-(dcbleft(5,1)*dcbright(2,1))
    v(624)=dcleft(2,1)*dcright(3,1)
    v(623)=dcleft(3,1)*dcright(2,1)
    v(622)=-(dcleft(5,1)*dcright(6,1))
    v(621)=-(dcleft(6,1)*dcright(5,1))
    v(620)=dcleft(3,1)*dcright(4,1)
    v(619)=dcleft(4,1)*dcright(3,1)
    v(618)=dcleft(4,1)*dcright(6,1)
    v(617)=-(dcleft(2,1)*dcright(5,1))
    v(616)=dcleft(6,1)*dcright(4,1)
    v(615)=-(dcleft(5,1)*dcright(2,1))
    v(614)=dcbleft(2,2)*dcbright(3,1)
    v(613)=dcbleft(3,2)*dcbright(2,1)
    v(612)=-(dcbleft(5,2)*dcbright(6,1))
    v(611)=-(dcbleft(6,2)*dcbright(5,1))
    v(610)=dcbleft(3,2)*dcbright(4,1)
    v(609)=dcbleft(4,2)*dcbright(3,1)
    v(608)=dcbleft(4,2)*dcbright(6,1)
    v(607)=-(dcbleft(2,2)*dcbright(5,1))
    v(606)=dcbleft(6,2)*dcbright(4,1)
    v(605)=-(dcbleft(5,2)*dcbright(2,1))
    v(604)=dcleft(2,2)*dcright(3,1)
    v(603)=dcleft(3,2)*dcright(2,1)
    v(602)=-(dcleft(5,2)*dcright(6,1))
    v(601)=-(dcleft(6,2)*dcright(5,1))
    v(600)=dcleft(3,2)*dcright(4,1)
    v(599)=dcleft(4,2)*dcright(3,1)
    v(598)=dcleft(4,2)*dcright(6,1)
    v(597)=-(dcleft(2,2)*dcright(5,1))
    v(596)=dcleft(6,2)*dcright(4,1)
    v(595)=-(dcleft(5,2)*dcright(2,1))
    v(594)=(-2d0)*dcbright(6,1)
    v(593)=(-2d0)*dcright(6,1)
    v(591)=dcbleft(6,3)*v(594)
    v(590)=dcbleft(2,3)*dcbright(3,1)
    v(589)=dcbleft(3,3)*dcbright(2,1)
    v(588)=-(dcbleft(5,3)*dcbright(6,1))
    v(587)=-(dcbleft(6,3)*dcbright(5,1))
    v(586)=dcbleft(3,3)*dcbright(4,1)
    v(585)=dcbleft(4,3)*dcbright(3,1)
    v(584)=dcbleft(4,3)*dcbright(6,1)
    v(583)=-(dcbleft(2,3)*dcbright(5,1))
    v(582)=dcbleft(6,3)*dcbright(4,1)
    v(581)=-(dcbleft(5,3)*dcbright(2,1))
    v(580)=dcleft(6,3)*v(593)
    v(579)=dcleft(2,3)*dcright(3,1)
    v(578)=dcleft(3,3)*dcright(2,1)
    v(577)=-(dcleft(5,3)*dcright(6,1))
    v(576)=-(dcleft(6,3)*dcright(5,1))
    v(575)=dcleft(3,3)*dcright(4,1)
    v(574)=dcleft(4,3)*dcright(3,1)
    v(573)=dcleft(4,3)*dcright(6,1)
    v(572)=-(dcleft(2,3)*dcright(5,1))
    v(571)=dcleft(6,3)*dcright(4,1)
    v(570)=-(dcleft(5,3)*dcright(2,1))
    v(569)=cb(2)*dcbright(3,3)
    v(568)=cb(3)*dcbright(2,3)
    v(567)=cb(2)*dcbright(3,2)
    v(566)=cb(3)*dcbright(2,2)
    v(565)=cb(2)*dcbright(3,1)
    v(564)=cb(3)*dcbright(2,1)
    v(563)=cb(4)*dcbright(3,3)+cb(3)*dcbright(4,3)-cb(6)*dcbright(5,3)-cb(5)*dcbright(6,3)
    v(562)=cb(4)*dcbright(3,2)+cb(3)*dcbright(4,2)-cb(6)*dcbright(5,2)-cb(5)*dcbright(6,2)
    v(561)=cb(4)*dcbright(3,1)+cb(3)*dcbright(4,1)-cb(6)*dcbright(5,1)-cb(5)*dcbright(6,1)
    v(560)=-(cb(5)*dcbright(2,3))+cb(6)*dcbright(4,3)-cb(2)*dcbright(5,3)+cb(4)*dcbright(6,3)
    v(559)=-(cb(5)*dcbright(2,2))+cb(6)*dcbright(4,2)-cb(2)*dcbright(5,2)+cb(4)*dcbright(6,2)
    v(558)=-(cb(5)*dcbright(2,1))+cb(6)*dcbright(4,1)-cb(2)*dcbright(5,1)+cb(4)*dcbright(6,1)
    v(557)=c(2)*dcright(3,3)
    v(556)=c(3)*dcright(2,3)
    v(555)=c(2)*dcright(3,2)
    v(554)=c(3)*dcright(2,2)
    v(553)=c(2)*dcright(3,1)
    v(552)=c(3)*dcright(2,1)
    v(551)=c(4)*dcright(3,3)+c(3)*dcright(4,3)-c(6)*dcright(5,3)-c(5)*dcright(6,3)
    v(550)=c(4)*dcright(3,2)+c(3)*dcright(4,2)-c(6)*dcright(5,2)-c(5)*dcright(6,2)
    v(549)=c(4)*dcright(3,1)+c(3)*dcright(4,1)-c(6)*dcright(5,1)-c(5)*dcright(6,1)
    v(548)=-(c(5)*dcright(2,3))+c(6)*dcright(4,3)-c(2)*dcright(5,3)+c(4)*dcright(6,3)
    v(547)=-(c(5)*dcright(2,2))+c(6)*dcright(4,2)-c(2)*dcright(5,2)+c(4)*dcright(6,2)
    v(546)=-(c(5)*dcright(2,1))+c(6)*dcright(4,1)-c(2)*dcright(5,1)+c(4)*dcright(6,1)
    v(544)=dcb2(6,3,3)
    v(543)=dcb2(6,3,2)
    v(542)=dcb2(6,3,1)
    v(541)=dcb2(6,2,3)
    v(540)=dcb2(6,2,2)
    v(539)=dcb2(6,2,1)
    v(538)=dcb2(6,1,3)
    v(537)=dcb2(6,1,2)
    v(536)=dcb2(6,1,1)
    v(535)=dcb2(5,3,3)
    v(534)=dcb2(5,3,2)
    v(533)=dcb2(5,3,1)
    v(532)=dcb2(5,2,3)
    v(531)=dcb2(5,2,2)
    v(530)=dcb2(5,2,1)
    v(529)=dcb2(5,1,3)
    v(528)=dcb2(5,1,2)
    v(527)=dcb2(5,1,1)
    v(526)=dcb2(4,3,3)
    v(525)=dcb2(4,3,2)
    v(524)=dcb2(4,3,1)
    v(523)=dcb2(4,2,3)
    v(522)=dcb2(4,2,2)
    v(521)=dcb2(4,2,1)
    v(520)=dcb2(4,1,3)
    v(519)=dcb2(4,1,2)
    v(518)=dcb2(4,1,1)
    v(517)=dcb2(3,3,3)
    v(516)=dcb2(3,3,2)
    v(515)=dcb2(3,3,1)
    v(514)=dcb2(3,2,3)
    v(513)=dcb2(3,2,2)
    v(512)=dcb2(3,2,1)
    v(511)=dcb2(3,1,3)
    v(510)=dcb2(3,1,2)
    v(509)=dcb2(3,1,1)
    v(508)=dcb2(2,3,3)
    v(507)=dcb2(2,3,2)
    v(506)=dcb2(2,3,1)
    v(505)=dcb2(2,2,3)
    v(504)=dcb2(2,2,2)
    v(503)=dcb2(2,2,1)
    v(502)=dcb2(2,1,3)
    v(501)=dcb2(2,1,2)
    v(500)=dcb2(2,1,1)
    v(499)=dc2(6,3,3)
    v(498)=dc2(6,3,2)
    v(497)=dc2(6,3,1)
    v(496)=dc2(6,2,3)
    v(495)=dc2(6,2,2)
    v(494)=dc2(6,2,1)
    v(493)=dc2(6,1,3)
    v(492)=dc2(6,1,2)
    v(491)=dc2(6,1,1)
    v(490)=dc2(5,3,3)
    v(489)=dc2(5,3,2)
    v(488)=dc2(5,3,1)
    v(487)=dc2(5,2,3)
    v(486)=dc2(5,2,2)
    v(485)=dc2(5,2,1)
    v(484)=dc2(5,1,3)
    v(483)=dc2(5,1,2)
    v(482)=dc2(5,1,1)
    v(481)=dc2(4,3,3)
    v(480)=dc2(4,3,2)
    v(479)=dc2(4,3,1)
    v(478)=dc2(4,2,3)
    v(477)=dc2(4,2,2)
    v(476)=dc2(4,2,1)
    v(475)=dc2(4,1,3)
    v(474)=dc2(4,1,2)
    v(473)=dc2(4,1,1)
    v(472)=dc2(3,3,3)
    v(471)=dc2(3,3,2)
    v(470)=dc2(3,3,1)
    v(469)=dc2(3,2,3)
    v(468)=dc2(3,2,2)
    v(467)=dc2(3,2,1)
    v(466)=dc2(3,1,3)
    v(465)=dc2(3,1,2)
    v(464)=dc2(3,1,1)
    v(463)=dc2(2,3,3)
    v(462)=dc2(2,3,2)
    v(461)=dc2(2,3,1)
    v(460)=dc2(2,2,3)
    v(459)=dc2(2,2,2)
    v(458)=dc2(2,2,1)
    v(457)=dc2(2,1,3)
    v(456)=dc2(2,1,2)
    v(455)=dc2(2,1,1)
    v(454)=dc2(1,3,3)
    v(453)=dc2(1,3,2)
    v(452)=dc2(1,3,1)
    v(451)=dc2(1,2,3)
    v(450)=dc2(1,2,2)
    v(449)=dc2(1,2,1)
    v(448)=dc2(1,1,3)
    v(447)=dc2(1,1,2)
    v(446)=dc2(1,1,1)
    v(445)=cb(2)*dcbleft(3,3)
    v(444)=cb(3)*dcbleft(2,3)
    v(443)=cb(4)*dcbleft(3,3)+cb(3)*dcbleft(4,3)-cb(6)*dcbleft(5,3)-cb(5)*dcbleft(6,3)
    v(442)=-(cb(5)*dcbleft(2,3))+cb(6)*dcbleft(4,3)-cb(2)*dcbleft(5,3)+cb(4)*dcbleft(6,3)
    v(441)=cb(2)*dcbleft(3,2)
    v(440)=cb(3)*dcbleft(2,2)
    v(439)=(-2d0)*cb(6)
    v(438)=cb(4)*dcbleft(3,2)+cb(3)*dcbleft(4,2)-cb(6)*dcbleft(5,2)-cb(5)*dcbleft(6,2)
    v(437)=-(cb(5)*dcbleft(2,2))+cb(6)*dcbleft(4,2)-cb(2)*dcbleft(5,2)+cb(4)*dcbleft(6,2)
    v(436)=cb(3)*dcbleft(2,1)+cb(2)*dcbleft(3,1)+dcbleft(6,1)*v(439)
    v(435)=cb(4)*dcbleft(3,1)+cb(3)*dcbleft(4,1)-cb(6)*dcbleft(5,1)-cb(5)*dcbleft(6,1)
    v(434)=-(cb(5)*dcbleft(2,1))+cb(6)*dcbleft(4,1)-cb(2)*dcbleft(5,1)+cb(4)*dcbleft(6,1)
    v(433)=c(2)*dcleft(3,3)
    v(432)=c(3)*dcleft(2,3)
    v(431)=c(4)*dcleft(3,3)+c(3)*dcleft(4,3)-c(6)*dcleft(5,3)-c(5)*dcleft(6,3)
    v(430)=-(c(5)*dcleft(2,3))+c(6)*dcleft(4,3)-c(2)*dcleft(5,3)+c(4)*dcleft(6,3)
    v(429)=c(2)*dcleft(3,2)
    v(428)=c(3)*dcleft(2,2)
    v(427)=(-2d0)*c(6)
    v(426)=c(4)*dcleft(3,2)+c(3)*dcleft(4,2)-c(6)*dcleft(5,2)-c(5)*dcleft(6,2)
    v(425)=-(c(5)*dcleft(2,2))+c(6)*dcleft(4,2)-c(2)*dcleft(5,2)+c(4)*dcleft(6,2)
    v(424)=c(3)*dcleft(2,1)+c(2)*dcleft(3,1)+dcleft(6,1)*v(427)
    v(423)=c(4)*dcleft(3,1)+c(3)*dcleft(4,1)-c(6)*dcleft(5,1)-c(5)*dcleft(6,1)
    v(422)=-(c(5)*dcleft(2,1))+c(6)*dcleft(4,1)-c(2)*dcleft(5,1)+c(4)*dcleft(6,1)
    v(420)=cb(2)*cb(3)-cb(6)**2
    v(419)=cb(3)*cb(4)-cb(5)*cb(6)
    v(418)=-(cb(2)*cb(5))+cb(4)*cb(6)
    v(417)=c(2)*c(3)-c(6)**2
    v(416)=c(3)*c(4)-c(5)*c(6)
    v(415)=-(c(2)*c(5))+c(4)*c(6)
    v(260)=c(5)*v(415)-c(4)*v(416)+c(1)*v(417)
    v(545)=1d0/(3d0*v(260))
    v(304)=1d0/v(260)**2
    v(635)=(-1d0/3d0)*v(304)
    v(261)=cb(5)*v(418)-cb(4)*v(419)+cb(1)*v(420)
    v(421)=v(261)/3d0
    v(321)=((2d0/3d0)*v(261))/v(260)**3
    v(269)=-(v(304)*v(421))
    v(254)=v(421)/v(260)
    v(334)=dcleft(5,1)*v(415)-dcleft(4,1)*v(416)+dcleft(1,1)*v(417)+c(5)*v(422)-c(4)*v(423)+c(1)*v(424)
    v(338)=dcleft(6,2)*v(427)+v(428)+v(429)
    v(345)=c(1)*v(338)+dcleft(5,2)*v(415)-dcleft(4,2)*v(416)+dcleft(1,2)*v(417)+c(5)*v(425)-c(4)*v(426)
    v(349)=dcleft(6,3)*v(427)+v(432)+v(433)
    v(356)=c(1)*v(349)+dcleft(5,3)*v(415)-dcleft(4,3)*v(416)+dcleft(1,3)*v(417)+c(5)*v(430)-c(4)*v(431)
    v(330)=dcbleft(5,1)*v(418)-dcbleft(4,1)*v(419)+dcbleft(1,1)*v(420)+cb(5)*v(434)-cb(4)*v(435)+cb(1)*v(436)
    v(342)=dcbleft(6,2)*v(439)+v(440)+v(441)
    v(341)=cb(1)*v(342)+dcbleft(5,2)*v(418)-dcbleft(4,2)*v(419)+dcbleft(1,2)*v(420)+cb(5)*v(437)-cb(4)*v(438)
    v(353)=dcbleft(6,3)*v(439)+v(444)+v(445)
    v(352)=cb(1)*v(353)+dcbleft(5,3)*v(418)-dcbleft(4,3)*v(419)+dcbleft(1,3)*v(420)+cb(5)*v(442)-cb(4)*v(443)
    v(268)=v(269)*v(334)+v(330)*v(545)
    v(271)=v(269)*v(345)+v(341)*v(545)
    v(273)=v(269)*v(356)+v(352)*v(545)
    v(298)=dcright(6,1)*v(427)+v(552)+v(553)
    v(299)=dcright(6,2)*v(427)+v(554)+v(555)
    v(300)=dcright(6,3)*v(427)+v(556)+v(557)
    v(301)=c(1)*v(298)+dcright(5,1)*v(415)-dcright(4,1)*v(416)+dcright(1,1)*v(417)+c(5)*v(546)-c(4)*v(549)
    v(302)=c(1)*v(299)+dcright(5,2)*v(415)-dcright(4,2)*v(416)+dcright(1,2)*v(417)+c(5)*v(547)-c(4)*v(550)
    v(303)=c(1)*v(300)+dcright(5,3)*v(415)-dcright(4,3)*v(416)+dcright(1,3)*v(417)+c(5)*v(548)-c(4)*v(551)
    v(305)=-(v(301)*v(304))
    v(592)=v(305)/3d0
    v(306)=-(v(302)*v(304))
    v(658)=v(306)/3d0
    v(307)=-(v(303)*v(304))
    v(723)=v(307)/3d0
    v(314)=dcbright(6,1)*v(439)+v(564)+v(565)
    v(315)=dcbright(6,2)*v(439)+v(566)+v(567)
    v(316)=dcbright(6,3)*v(439)+v(568)+v(569)
    v(317)=cb(1)*v(314)+dcbright(5,1)*v(418)-dcbright(4,1)*v(419)+dcbright(1,1)*v(420)+cb(5)*v(558)-cb(4)*v(561)
    v(318)=cb(1)*v(315)+dcbright(5,2)*v(418)-dcbright(4,2)*v(419)+dcbright(1,2)*v(420)+cb(5)*v(559)-cb(4)*v(562)
    v(319)=cb(1)*v(316)+dcbright(5,3)*v(418)-dcbright(4,3)*v(419)+dcbright(1,3)*v(420)+cb(5)*v(560)-cb(4)*v(563)
    v(320)=v(301)*v(321)+v(317)*v(635)
    v(357)=v(320)*v(356)+v(269)*(dcleft(1,3)*v(298)+dcright(1,1)*v(349)+dcright(5,1)*v(430)-dcright(4,1)*v(431)+v(417)*v&
         &(452)-v(416)*v(479)+v(415)*v(488)+dcleft(5,3)*v(546)-dcleft(4,3)*v(549)+c(5)*(-(c(5)*v(461))+c(6)*v(479)-c(2)*v(488)+c&
         &(4)*v(497)+v(570)+v(571)+v(572)+v(573))-c(4)*(c(4)*v(470)+c(3)*v(479)-c(6)*v(488)-c(5)*v(497)+v(574)+v(575)+v(576)+v&
         &(577))+c(1)*(c(3)*v(461)+c(2)*v(470)+v(427)*v(497)+v(578)+v(579)+v(580)))+v(545)*(dcbleft(1,3)*v(314)+dcbright(1,1)*v&
         &(353)+dcb2(1,3,1)*v(420)+dcbright(5,1)*v(442)-dcbright(4,1)*v(443)-v(419)*v(524)+v(418)*v(533)+dcbleft(5,3)*v(558)&
         &-dcbleft(4,3)*v(561)+cb(5)*(-(cb(5)*v(506))+cb(6)*v(524)-cb(2)*v(533)+cb(4)*v(542)+v(581)+v(582)+v(583)+v(584))-cb(4)*&
         &(cb(4)*v(515)+cb(3)*v(524)-cb(6)*v(533)-cb(5)*v(542)+v(585)+v(586)+v(587)+v(588))+cb(1)*(cb(3)*v(506)+cb(2)*v(515)+v&
         &(439)*v(542)+v(589)+v(590)+v(591)))+v(352)*v(592)
    v(346)=v(320)*v(345)+v(341)*v(592)+v(269)*(dcleft(1,2)*v(298)+dcright(1,1)*v(338)+dcright(5,1)*v(425)-dcright(4,1)*v&
         &(426)+v(417)*v(449)-v(416)*v(476)+v(415)*v(485)+dcleft(5,2)*v(546)-dcleft(4,2)*v(549)+c(5)*(-(c(5)*v(458))+c(6)*v(476)&
         &-c(2)*v(485)+c(4)*v(494)+v(595)+v(596)+v(597)+v(598))-c(4)*(c(4)*v(467)+c(3)*v(476)-c(6)*v(485)-c(5)*v(494)+v(599)+v&
         &(600)+v(601)+v(602))+c(1)*(c(3)*v(458)+c(2)*v(467)+v(427)*v(494)+dcleft(6,2)*v(593)+v(603)+v(604)))+v(545)*(dcbleft(1,2&
         &)*v(314)+dcbright(1,1)*v(342)+dcb2(1,2,1)*v(420)+dcbright(5,1)*v(437)-dcbright(4,1)*v(438)-v(419)*v(521)+v(418)*v(530)&
         &+dcbleft(5,2)*v(558)-dcbleft(4,2)*v(561)+cb(5)*(-(cb(5)*v(503))+cb(6)*v(521)-cb(2)*v(530)+cb(4)*v(539)+v(605)+v(606)+v&
         &(607)+v(608))-cb(4)*(cb(4)*v(512)+cb(3)*v(521)-cb(6)*v(530)-cb(5)*v(539)+v(609)+v(610)+v(611)+v(612))+cb(1)*(cb(3)*v&
         &(503)+cb(2)*v(512)+v(439)*v(539)+dcbleft(6,2)*v(594)+v(613)+v(614)))
    v(335)=v(320)*v(334)+v(330)*v(592)+v(269)*(dcleft(1,1)*v(298)+dcright(5,1)*v(422)-dcright(4,1)*v(423)+dcright(1,1)*v&
         &(424)+v(417)*v(446)-v(416)*v(473)+v(415)*v(482)+dcleft(5,1)*v(546)-dcleft(4,1)*v(549)+c(5)*(-(c(5)*v(455))+c(6)*v(473)&
         &-c(2)*v(482)+c(4)*v(491)+v(615)+v(616)+v(617)+v(618))-c(4)*(c(4)*v(464)+c(3)*v(473)-c(6)*v(482)-c(5)*v(491)+v(619)+v&
         &(620)+v(621)+v(622))+c(1)*(c(3)*v(455)+c(2)*v(464)+v(427)*v(491)+dcleft(6,1)*v(593)+v(623)+v(624)))+v(545)*(dcbleft(1,1&
         &)*v(314)+dcb2(1,1,1)*v(420)+dcbright(5,1)*v(434)-dcbright(4,1)*v(435)+dcbright(1,1)*v(436)-v(419)*v(518)+v(418)*v(527)&
         &+dcbleft(5,1)*v(558)-dcbleft(4,1)*v(561)+cb(5)*(-(cb(5)*v(500))+cb(6)*v(518)-cb(2)*v(527)+cb(4)*v(536)+v(625)+v(626)+v&
         &(627)+v(628))-cb(4)*(cb(4)*v(509)+cb(3)*v(518)-cb(6)*v(527)-cb(5)*v(536)+v(629)+v(630)+v(631)+v(632))+cb(1)*(cb(3)*v&
         &(500)+cb(2)*v(509)+v(439)*v(536)+dcbleft(6,1)*v(594)+v(633)+v(634)))
    v(322)=v(302)*v(321)+v(318)*v(635)
    v(358)=v(322)*v(356)+v(269)*(dcleft(1,3)*v(299)+dcright(1,2)*v(349)+dcright(5,2)*v(430)-dcright(4,2)*v(431)+v(417)*v&
         &(453)-v(416)*v(480)+v(415)*v(489)+dcleft(5,3)*v(547)-dcleft(4,3)*v(550)+c(5)*(-(c(5)*v(462))+c(6)*v(480)-c(2)*v(489)+c&
         &(4)*v(498)+v(636)+v(637)+v(638)+v(639))-c(4)*(c(4)*v(471)+c(3)*v(480)-c(6)*v(489)-c(5)*v(498)+v(640)+v(641)+v(642)+v&
         &(643))+c(1)*(c(3)*v(462)+c(2)*v(471)+v(427)*v(498)+v(644)+v(645)+v(646)))+v(545)*(dcbleft(1,3)*v(315)+dcbright(1,2)*v&
         &(353)+dcb2(1,3,2)*v(420)+dcbright(5,2)*v(442)-dcbright(4,2)*v(443)-v(419)*v(525)+v(418)*v(534)+dcbleft(5,3)*v(559)&
         &-dcbleft(4,3)*v(562)+cb(5)*(-(cb(5)*v(507))+cb(6)*v(525)-cb(2)*v(534)+cb(4)*v(543)+v(647)+v(648)+v(649)+v(650))-cb(4)*&
         &(cb(4)*v(516)+cb(3)*v(525)-cb(6)*v(534)-cb(5)*v(543)+v(651)+v(652)+v(653)+v(654))+cb(1)*(cb(3)*v(507)+cb(2)*v(516)+v&
         &(439)*v(543)+v(655)+v(656)+v(657)))+v(352)*v(658)
    v(347)=v(322)*v(345)+v(341)*v(658)+v(269)*(dcleft(1,2)*v(299)+dcright(1,2)*v(338)+dcright(5,2)*v(425)-dcright(4,2)*v&
         &(426)+v(417)*v(450)-v(416)*v(477)+v(415)*v(486)+dcleft(5,2)*v(547)-dcleft(4,2)*v(550)+c(5)*(-(c(5)*v(459))+c(6)*v(477)&
         &-c(2)*v(486)+c(4)*v(495)+v(661)+v(662)+v(663)+v(664))-c(4)*(c(4)*v(468)+c(3)*v(477)-c(6)*v(486)-c(5)*v(495)+v(665)+v&
         &(666)+v(667)+v(668))+c(1)*(c(3)*v(459)+c(2)*v(468)+v(427)*v(495)+dcleft(6,2)*v(659)+v(669)+v(670)))+v(545)*(dcbleft(1,2&
         &)*v(315)+dcbright(1,2)*v(342)+dcb2(1,2,2)*v(420)+dcbright(5,2)*v(437)-dcbright(4,2)*v(438)-v(419)*v(522)+v(418)*v(531)&
         &+dcbleft(5,2)*v(559)-dcbleft(4,2)*v(562)+cb(5)*(-(cb(5)*v(504))+cb(6)*v(522)-cb(2)*v(531)+cb(4)*v(540)+v(671)+v(672)+v&
         &(673)+v(674))-cb(4)*(cb(4)*v(513)+cb(3)*v(522)-cb(6)*v(531)-cb(5)*v(540)+v(675)+v(676)+v(677)+v(678))+cb(1)*(cb(3)*v&
         &(504)+cb(2)*v(513)+v(439)*v(540)+dcbleft(6,2)*v(660)+v(679)+v(680)))
    v(336)=v(322)*v(334)+v(330)*v(658)+v(269)*(dcleft(1,1)*v(299)+dcright(5,2)*v(422)-dcright(4,2)*v(423)+dcright(1,2)*v&
         &(424)+v(417)*v(447)-v(416)*v(474)+v(415)*v(483)+dcleft(5,1)*v(547)-dcleft(4,1)*v(550)+c(5)*(-(c(5)*v(456))+c(6)*v(474)&
         &-c(2)*v(483)+c(4)*v(492)+v(681)+v(682)+v(683)+v(684))-c(4)*(c(4)*v(465)+c(3)*v(474)-c(6)*v(483)-c(5)*v(492)+v(685)+v&
         &(686)+v(687)+v(688))+c(1)*(c(3)*v(456)+c(2)*v(465)+v(427)*v(492)+dcleft(6,1)*v(659)+v(689)+v(690)))+v(545)*(dcbleft(1,1&
         &)*v(315)+dcb2(1,1,2)*v(420)+dcbright(5,2)*v(434)-dcbright(4,2)*v(435)+dcbright(1,2)*v(436)-v(419)*v(519)+v(418)*v(528)&
         &+dcbleft(5,1)*v(559)-dcbleft(4,1)*v(562)+cb(5)*(-(cb(5)*v(501))+cb(6)*v(519)-cb(2)*v(528)+cb(4)*v(537)+v(691)+v(692)+v&
         &(693)+v(694))-cb(4)*(cb(4)*v(510)+cb(3)*v(519)-cb(6)*v(528)-cb(5)*v(537)+v(695)+v(696)+v(697)+v(698))+cb(1)*(cb(3)*v&
         &(501)+cb(2)*v(510)+v(439)*v(537)+dcbleft(6,1)*v(660)+v(699)+v(700)))
    v(323)=v(303)*v(321)+v(319)*v(635)
    v(359)=v(323)*v(356)+v(269)*(dcleft(1,3)*v(300)+dcright(1,3)*v(349)+dcright(5,3)*v(430)-dcright(4,3)*v(431)+v(417)*v&
         &(454)-v(416)*v(481)+v(415)*v(490)+dcleft(5,3)*v(548)-dcleft(4,3)*v(551)+c(5)*(-(c(5)*v(463))+c(6)*v(481)-c(2)*v(490)+c&
         &(4)*v(499)+v(701)+v(702)+v(703)+v(704))-c(4)*(c(4)*v(472)+c(3)*v(481)-c(6)*v(490)-c(5)*v(499)+v(705)+v(706)+v(707)+v&
         &(708))+c(1)*(c(3)*v(463)+c(2)*v(472)+v(427)*v(499)+v(709)+v(710)+v(711)))+v(545)*(dcbleft(1,3)*v(316)+dcbright(1,3)*v&
         &(353)+dcb2(1,3,3)*v(420)+dcbright(5,3)*v(442)-dcbright(4,3)*v(443)-v(419)*v(526)+v(418)*v(535)+dcbleft(5,3)*v(560)&
         &-dcbleft(4,3)*v(563)+cb(5)*(-(cb(5)*v(508))+cb(6)*v(526)-cb(2)*v(535)+cb(4)*v(544)+v(712)+v(713)+v(714)+v(715))-cb(4)*&
         &(cb(4)*v(517)+cb(3)*v(526)-cb(6)*v(535)-cb(5)*v(544)+v(716)+v(717)+v(718)+v(719))+cb(1)*(cb(3)*v(508)+cb(2)*v(517)+v&
         &(439)*v(544)+v(720)+v(721)+v(722)))+v(352)*v(723)
    v(348)=v(323)*v(345)+v(341)*v(723)+v(269)*(dcleft(1,2)*v(300)+dcright(1,3)*v(338)+dcright(5,3)*v(425)-dcright(4,3)*v&
         &(426)+v(417)*v(451)-v(416)*v(478)+v(415)*v(487)+dcleft(5,2)*v(548)-dcleft(4,2)*v(551)+c(5)*(-(c(5)*v(460))+c(6)*v(478)&
         &-c(2)*v(487)+c(4)*v(496)+v(726)+v(727)+v(728)+v(729))-c(4)*(c(4)*v(469)+c(3)*v(478)-c(6)*v(487)-c(5)*v(496)+v(730)+v&
         &(731)+v(732)+v(733))+c(1)*(c(3)*v(460)+c(2)*v(469)+v(427)*v(496)+dcleft(6,2)*v(724)+v(734)+v(735)))+v(545)*(dcbleft(1,2&
         &)*v(316)+dcbright(1,3)*v(342)+dcb2(1,2,3)*v(420)+dcbright(5,3)*v(437)-dcbright(4,3)*v(438)-v(419)*v(523)+v(418)*v(532)&
         &+dcbleft(5,2)*v(560)-dcbleft(4,2)*v(563)+cb(5)*(-(cb(5)*v(505))+cb(6)*v(523)-cb(2)*v(532)+cb(4)*v(541)+v(736)+v(737)+v&
         &(738)+v(739))-cb(4)*(cb(4)*v(514)+cb(3)*v(523)-cb(6)*v(532)-cb(5)*v(541)+v(740)+v(741)+v(742)+v(743))+cb(1)*(cb(3)*v&
         &(505)+cb(2)*v(514)+v(439)*v(541)+dcbleft(6,2)*v(725)+v(744)+v(745)))
    v(337)=v(323)*v(334)+v(330)*v(723)+v(269)*(dcleft(1,1)*v(300)+dcright(5,3)*v(422)-dcright(4,3)*v(423)+dcright(1,3)*v&
         &(424)+v(417)*v(448)-v(416)*v(475)+v(415)*v(484)+dcleft(5,1)*v(548)-dcleft(4,1)*v(551)+c(5)*(-(c(5)*v(457))+c(6)*v(475)&
         &-c(2)*v(484)+c(4)*v(493)+v(746)+v(747)+v(748)+v(749))-c(4)*(c(4)*v(466)+c(3)*v(475)-c(6)*v(484)-c(5)*v(493)+v(750)+v&
         &(751)+v(752)+v(753))+c(1)*(c(3)*v(457)+c(2)*v(466)+v(427)*v(493)+dcleft(6,1)*v(724)+v(754)+v(755)))+v(545)*(dcbleft(1,1&
         &)*v(316)+dcb2(1,1,3)*v(420)+dcbright(5,3)*v(434)-dcbright(4,3)*v(435)+dcbright(1,3)*v(436)-v(419)*v(520)+v(418)*v(529)&
         &+dcbleft(5,1)*v(560)-dcbleft(4,1)*v(563)+cb(5)*(-(cb(5)*v(502))+cb(6)*v(520)-cb(2)*v(529)+cb(4)*v(538)+v(756)+v(757)+v&
         &(758)+v(759))-cb(4)*(cb(4)*v(511)+cb(3)*v(520)-cb(6)*v(529)-cb(5)*v(538)+v(760)+v(761)+v(762)+v(763))+cb(1)*(cb(3)*v&
         &(502)+cb(2)*v(511)+v(439)*v(538)+dcbleft(6,1)*v(725)+v(764)+v(765)))
    v(324)=(v(261)*v(305)+v(317)/v(260))/3d0
    v(325)=(v(261)*v(306)+v(318)/v(260))/3d0
    v(326)=(v(261)*v(307)+v(319)/v(260))/3d0
    dcs2(1,1,1)=dcright(1,1)*v(268)+dcleft(1,1)*v(324)+c(1)*v(335)+v(254)*v(446)
    dcs2(1,1,2)=dcright(1,2)*v(268)+dcleft(1,1)*v(325)+c(1)*v(336)+v(254)*v(447)
    dcs2(1,1,3)=dcright(1,3)*v(268)+dcleft(1,1)*v(326)+c(1)*v(337)+v(254)*v(448)
    dcs2(1,2,1)=dcright(1,1)*v(271)+dcleft(1,2)*v(324)+c(1)*v(346)+v(254)*v(449)
    dcs2(1,2,2)=dcright(1,2)*v(271)+dcleft(1,2)*v(325)+c(1)*v(347)+v(254)*v(450)
    dcs2(1,2,3)=dcright(1,3)*v(271)+dcleft(1,2)*v(326)+c(1)*v(348)+v(254)*v(451)
    dcs2(1,3,1)=dcright(1,1)*v(273)+dcleft(1,3)*v(324)+c(1)*v(357)+v(254)*v(452)
    dcs2(1,3,2)=dcright(1,2)*v(273)+dcleft(1,3)*v(325)+c(1)*v(358)+v(254)*v(453)
    dcs2(1,3,3)=dcright(1,3)*v(273)+dcleft(1,3)*v(326)+c(1)*v(359)+v(254)*v(454)
    dcs2(2,1,1)=dcright(2,1)*v(268)+dcleft(2,1)*v(324)+c(2)*v(335)+v(254)*v(455)
    dcs2(2,1,2)=dcright(2,2)*v(268)+dcleft(2,1)*v(325)+c(2)*v(336)+v(254)*v(456)
    dcs2(2,1,3)=dcright(2,3)*v(268)+dcleft(2,1)*v(326)+c(2)*v(337)+v(254)*v(457)
    dcs2(2,2,1)=dcright(2,1)*v(271)+dcleft(2,2)*v(324)+c(2)*v(346)+v(254)*v(458)
    dcs2(2,2,2)=dcright(2,2)*v(271)+dcleft(2,2)*v(325)+c(2)*v(347)+v(254)*v(459)
    dcs2(2,2,3)=dcright(2,3)*v(271)+dcleft(2,2)*v(326)+c(2)*v(348)+v(254)*v(460)
    dcs2(2,3,1)=dcright(2,1)*v(273)+dcleft(2,3)*v(324)+c(2)*v(357)+v(254)*v(461)
    dcs2(2,3,2)=dcright(2,2)*v(273)+dcleft(2,3)*v(325)+c(2)*v(358)+v(254)*v(462)
    dcs2(2,3,3)=dcright(2,3)*v(273)+dcleft(2,3)*v(326)+c(2)*v(359)+v(254)*v(463)
    dcs2(3,1,1)=dcright(3,1)*v(268)+dcleft(3,1)*v(324)+c(3)*v(335)+v(254)*v(464)
    dcs2(3,1,2)=dcright(3,2)*v(268)+dcleft(3,1)*v(325)+c(3)*v(336)+v(254)*v(465)
    dcs2(3,1,3)=dcright(3,3)*v(268)+dcleft(3,1)*v(326)+c(3)*v(337)+v(254)*v(466)
    dcs2(3,2,1)=dcright(3,1)*v(271)+dcleft(3,2)*v(324)+c(3)*v(346)+v(254)*v(467)
    dcs2(3,2,2)=dcright(3,2)*v(271)+dcleft(3,2)*v(325)+c(3)*v(347)+v(254)*v(468)
    dcs2(3,2,3)=dcright(3,3)*v(271)+dcleft(3,2)*v(326)+c(3)*v(348)+v(254)*v(469)
    dcs2(3,3,1)=dcright(3,1)*v(273)+dcleft(3,3)*v(324)+c(3)*v(357)+v(254)*v(470)
    dcs2(3,3,2)=dcright(3,2)*v(273)+dcleft(3,3)*v(325)+c(3)*v(358)+v(254)*v(471)
    dcs2(3,3,3)=dcright(3,3)*v(273)+dcleft(3,3)*v(326)+c(3)*v(359)+v(254)*v(472)
    dcs2(4,1,1)=dcright(4,1)*v(268)+dcleft(4,1)*v(324)+c(4)*v(335)+v(254)*v(473)
    dcs2(4,1,2)=dcright(4,2)*v(268)+dcleft(4,1)*v(325)+c(4)*v(336)+v(254)*v(474)
    dcs2(4,1,3)=dcright(4,3)*v(268)+dcleft(4,1)*v(326)+c(4)*v(337)+v(254)*v(475)
    dcs2(4,2,1)=dcright(4,1)*v(271)+dcleft(4,2)*v(324)+c(4)*v(346)+v(254)*v(476)
    dcs2(4,2,2)=dcright(4,2)*v(271)+dcleft(4,2)*v(325)+c(4)*v(347)+v(254)*v(477)
    dcs2(4,2,3)=dcright(4,3)*v(271)+dcleft(4,2)*v(326)+c(4)*v(348)+v(254)*v(478)
    dcs2(4,3,1)=dcright(4,1)*v(273)+dcleft(4,3)*v(324)+c(4)*v(357)+v(254)*v(479)
    dcs2(4,3,2)=dcright(4,2)*v(273)+dcleft(4,3)*v(325)+c(4)*v(358)+v(254)*v(480)
    dcs2(4,3,3)=dcright(4,3)*v(273)+dcleft(4,3)*v(326)+c(4)*v(359)+v(254)*v(481)
    dcs2(5,1,1)=dcright(5,1)*v(268)+dcleft(5,1)*v(324)+c(5)*v(335)+v(254)*v(482)
    dcs2(5,1,2)=dcright(5,2)*v(268)+dcleft(5,1)*v(325)+c(5)*v(336)+v(254)*v(483)
    dcs2(5,1,3)=dcright(5,3)*v(268)+dcleft(5,1)*v(326)+c(5)*v(337)+v(254)*v(484)
    dcs2(5,2,1)=dcright(5,1)*v(271)+dcleft(5,2)*v(324)+c(5)*v(346)+v(254)*v(485)
    dcs2(5,2,2)=dcright(5,2)*v(271)+dcleft(5,2)*v(325)+c(5)*v(347)+v(254)*v(486)
    dcs2(5,2,3)=dcright(5,3)*v(271)+dcleft(5,2)*v(326)+c(5)*v(348)+v(254)*v(487)
    dcs2(5,3,1)=dcright(5,1)*v(273)+dcleft(5,3)*v(324)+c(5)*v(357)+v(254)*v(488)
    dcs2(5,3,2)=dcright(5,2)*v(273)+dcleft(5,3)*v(325)+c(5)*v(358)+v(254)*v(489)
    dcs2(5,3,3)=dcright(5,3)*v(273)+dcleft(5,3)*v(326)+c(5)*v(359)+v(254)*v(490)
    dcs2(6,1,1)=dcright(6,1)*v(268)+dcleft(6,1)*v(324)+c(6)*v(335)+v(254)*v(491)
    dcs2(6,1,2)=dcright(6,2)*v(268)+dcleft(6,1)*v(325)+c(6)*v(336)+v(254)*v(492)
    dcs2(6,1,3)=dcright(6,3)*v(268)+dcleft(6,1)*v(326)+c(6)*v(337)+v(254)*v(493)
    dcs2(6,2,1)=dcright(6,1)*v(271)+dcleft(6,2)*v(324)+c(6)*v(346)+v(254)*v(494)
    dcs2(6,2,2)=dcright(6,2)*v(271)+dcleft(6,2)*v(325)+c(6)*v(347)+v(254)*v(495)
    dcs2(6,2,3)=dcright(6,3)*v(271)+dcleft(6,2)*v(326)+c(6)*v(348)+v(254)*v(496)
    dcs2(6,3,1)=dcright(6,1)*v(273)+dcleft(6,3)*v(324)+c(6)*v(357)+v(254)*v(497)
    dcs2(6,3,2)=dcright(6,2)*v(273)+dcleft(6,3)*v(325)+c(6)*v(358)+v(254)*v(498)
    dcs2(6,3,3)=dcright(6,3)*v(273)+dcleft(6,3)*v(326)+c(6)*v(359)+v(254)*v(499)
  END SUBROUTINE mls_combinedfbar2old

!----------------------
!*** END MLS
!*** START RBF MESHLESS
!----------------------

!---------------------------
!*** RADIUS AND DERIVATIVES
!*** 
!*** DETERMINED
!*** CHK0,1
!---------------------------
  subroutine rbf_radius(ndi,x,xi,r,drdx,d2rdx2)
    IMPLICIT REAL(8) (a-h,o-z)
    REAL(8),PARAMETER::small=1.0d-10
    real(8),dimension(ndi)::x,xi
    real(8),dimension(ndi)::drdx
    real(8),dimension(ndi,ndi)::d2rdx2
    r=rnorm2(ndi,x-xi)
    drdx=0.0d00
    d2rdx2=0.0d00
    if(r.gt.small)then
       do id=1,ndi
          drdx(id)=(x(id)-xi(id))/r
          do jd=1,ndi
             d2rdx2(id,jd)=deltak(id,jd)/r-(1.0d00/(r**3.0d00))*(x(id)-xi(id))*(x(jd)-xi(jd))
          end do
       end do
    end if
  end subroutine rbf_radius

!--------------------------
!*** radial basis function
!*** basis function
!*** chk0
!--------------------------
  subroutine rbf_c2(rmax,r,a,dadr,d2adr2)
    implicit real(8) (a-h,o-z)
    real(8)::rmax,r,a,dadr,d2adr2
    a=0.0d00
    dadr=0.0d00
    d2adr2=0.0d00
    IF(r.LE.0.9999d00*rmax)THEN
       a=((1.0d00-(r/rmax))**4)*(1.0d00+4.0d00*(r/rmax))
       dadr=(20.0d00*r*(r-rmax)**3.0d00)/(rmax**5.0d00)
       d2adr2=20.0d00*((r-rmax)**2.0d00)*(4.0d00*r-rmax)/(rmax**5.0d00)
    end if
  end subroutine rbf_c2

!------------------------------------
!*** radial and derivatives
!*** xi is the "node" coordinate
!*** x is the gauss point coordinate
!*** chk0
!------------------------------------
  subroutine rbf_one(rmax,ndi,x,xi,a,dadx,d2adx2)
    implicit real(8) (a-h,o-z)
    real(8)::rmax,r,a,dadr,d2adr2
    real(8),dimension(ndi)::x,xi
    real(8),dimension(ndi)::dadx,drdx
    real(8),dimension(ndi,ndi)::d2adx2,d2rdx2
    call rbf_radius(ndi,x,xi,r,drdx,d2rdx2)
    call rbf_c2(rmax,r,a,dadr,d2adr2)
    do id=1,ndi
       dadx(id)=dadr*drdx(id)
       do jd=1,ndi
          d2adx2(id,jd)=d2adr2*drdx(id)*drdx(jd)+dadr*d2rdx2(id,jd)
       end do
    end do
  end subroutine rbf_one

!-------------------------------------------
!*** obtain shape functions and derivatives
!*** from a given list of nodes
!*** at one point with coordinates x
!*** chk0
!-------------------------------------------
  subroutine rbf_shapefunctions(imeshless,rmax,ndi,n,xn,x,xc,ff,dff,dff2)
    implicit real(8) (a-h,o-z)
    integer,parameter::mpol=20
    real(8)::rmax
    real(8),dimension(n,n)::amatrix,invmatrix,gmatrix
    real(8),dimension(n,mpol)::bmatrix
    real(8),dimension(mpol,n)::hmatrix
    real(8),dimension(ndi,ndi,n)::dff2
    real(8),dimension(ndi,n)::dff
    real(8),dimension(n)::ff,a
    real(8),dimension(ndi,n)::xn
    real(8),dimension(ndi)::x,xc
    real(8),dimension(mpol)::b
    real(8),dimension(ndi,mpol)::dbdx
    real(8),dimension(ndi,ndi,mpol)::dbdx2
    real(8),dimension(ndi,ndi,mpol)::d2bdx2
    real(8),dimension(mpol,mpol)::b1b,invb1b
    real(8),dimension(ndi,n)::dadx
    real(8),dimension(ndi,ndi,n)::d2adx2
    real(8),dimension(mpol)::polyn
    real(8),dimension(ndi,mpol)::dpolyn
    real(8),dimension(ndi,ndi,mpol)::dpolyn2
!------------------
!*** forms amatrix
!------------------
    do jn=1,n
       do in=1,n
!-----------------------------------
!*** dadx and d2adx2 are trash here
!-----------------------------------
          call rbf_one(rmax,ndi,xn(1:ndi,in),xn(1:ndi,jn),amatrix(in,jn),dadx(1:ndi,in),d2adx2(1:ndi,1:ndi,in))
       end do
    end do
!-------------------------
!*** invmatrix=amatrix^-1
!-------------------------
    call invmat(n,deta,amatrix,invmatrix)
!------------------
!*** forms bmatrix
!------------------
    do in=1,n
       mp=0
       select case(imeshless)
       case(1)
          call mls_polynlinbase(ndi,1.0d00,xn(1:ndi,in),xc(1:ndi),polyn(1:mpol),dpolyn(1:ndi,1:mpol),dpolyn2(1:ndi,1:ndi,1:mpol),mp)
       case(2)
          call mls_polynqbase(ndi,1.0d00,xn(1:ndi,in),xc(1:ndi),polyn(1:mpol),dpolyn(1:ndi,1:mpol),dpolyn2(1:ndi,1:ndi,1:mpol),mp)
       case(3)
          call mls_polyncubbase(ndi,1.0d00,xn(1:ndi,in),xc(1:ndi),polyn(1:mpol),dpolyn(1:ndi,1:mpol),dpolyn2(1:ndi,1:ndi,1:mpol),mp)
       end select
       do ip=1,mp
          bmatrix(in,ip)=polyn(ip)
       end do
    end do
!-----------------------
!*** forms big matrices
!-----------------------
    do jd=1,mp
       do id=1,mp
          b1b(id,jd)=0.0d00
          do in=1,n
             do jn=1,n
                b1b(id,jd)=b1b(id,jd)+bmatrix(in,id)*invmatrix(in,jn)*bmatrix(jn,jd)
             end do
          end do
       end do
    end do
!----------------
!*** invert b1b
!*** into invb1b
!----------------
    call invmat(mp,det,b1b(1:mp,1:mp),invb1b(1:mp,1:mp))
!------------
!*** hmatrix
!------------
    do in=1,n
       do id=1,mp
          hmatrix(id,in)=0.0d00          
          do jn=1,n
             do jd=1,mp
                hmatrix(id,in)=hmatrix(id,in)+invb1b(id,jd)*bmatrix(jn,jd)*invmatrix(jn,in)
             end do
          end do
       end do
    end do
!------------
!*** gmatrix
!------------
    do jn=1,n
       do in=1,n
          gmatrix(in,jn)=invmatrix(in,jn)
          do id=1,mp
             do kn=1,n
                gmatrix(in,jn)=gmatrix(in,jn)-invmatrix(in,kn)*bmatrix(kn,id)*hmatrix(id,jn)
             end do
          end do
       end do
    end do
!--------------------------
!*** forms vectors a and b
!--------------------------
    do in=1,n
       call rbf_one(rmax,ndi,x(1:ndi),xn(1:ndi,in),a(in),dadx(1:ndi,in),d2adx2(1:ndi,1:ndi,in))
       select case(imeshless)
       case(1)
          call mls_polynlinbase(ndi,1.0d00,x(1:ndi),xc(1:ndi),b,dbdx(1:ndi,1:mp),dbdx2(1:ndi,1:ndi,1:mp),mp)
       case(2)
          call mls_polynqbase(ndi,1.0d00,x(1:ndi),xc(1:ndi),b,dbdx(1:ndi,1:mp),dbdx2(1:ndi,1:ndi,1:mp),mp)
       case(3)
          call mls_polyncubbase(ndi,1.0d00,x(1:ndi),xc(1:ndi),b,dbdx(1:ndi,1:mp),dbdx2(1:ndi,1:ndi,1:mp),mp)
       end select
    end do
!------------------------------------------------------
!*** now determines rb shape functions and derivatives
!------------------------------------------------------
    if(imeshless.eq.0)then
       do in=1,n
          ff(in)=dotprod(n,a(1:n),invmatrix(1:n,in))
          do id=1,ndi
             dff(id,in)=dotprod(n,dadx(id,1:n),invmatrix(1:n,in))
             do jd=1,ndi
                dff2(id,jd,in)=dotprod(n,d2adx2(id,jd,1:n),invmatrix(1:n,in))
             enddo
          enddo
       enddo
    else
       do in=1,n
          ff(in)=dotprod(n,a(1:n),gmatrix(1:n,in))+dotprod(mp,b(1:mp),hmatrix(1:mp,in))
          DO id=1,ndi
             dff(id,in)=dotprod(n,dadx(id,1:n),gmatrix(1:n,in))+dotprod(mp,dbdx(id,1:mp),hmatrix(1:mp,in))
             do jd=1,ndi
                dff2(id,jd,in)=dotprod(n,d2adx2(id,jd,1:n),gmatrix(1:n,in))+dotprod(mp,dbdx2(id,jd,1:mp),hmatrix(1:mp,in))
             enddo
          enddo
       enddo
    end if
  end subroutine rbf_shapefunctions

!-------------------------------------------------
!*** from an element support estimate, determines
!*** gauss point support radius
!-------------------------------------------------
  subroutine rbf_gausslevelshapefunctions(imeshless,rcentroid,ndi,xg,xc,ntot,xtot,ff,dff,dff2)
    implicit real(8) (a-h,o-z)
    integer::imeshless
    real(8)::rcentroid,rgauss
    real(8),dimension(ndi)::xg,xc
    real(8),dimension(ntot)::ff
    real(8),dimension(ndi,ntot)::dff
    real(8),dimension(ndi,ndi,ntot)::dff2
    real(8),dimension(ndi,ntot)::xtot
    rgauss=rcentroid!-rnorm2(ndi,xg-xc)    
    call rbf_testingallrequired(imeshless,rgauss,ndi,ntot,xtot,xg,xc,ff,dff,dff2)    
  end subroutine rbf_gausslevelshapefunctions

!--------------------------------------------------------------------
!*** for representation purposes only:
!*** get all nodes, makes a selection
!*** and returns all quantities for all nodes, not only the selection
!--------------------------------------------------------------------
  SUBROUTINE rbf_testingallrequired(imeshless,rmax,ndi,ntot,xtot,x,xc,ff,dff,dff2)
    implicit real(8) (a-h,o-z)
    integer::imeshless
    integer,parameter::mpol=20
    integer,dimension(:),allocatable::listn
    real(8)::rmax
    real(8),dimension(ntot)::ff,ffloc
    real(8),dimension(ndi,ntot)::dff,dffloc
    real(8),dimension(ndi,ndi,ntot)::dff2,dff2loc
    real(8),dimension(ndi,ntot)::xtot,xn
    real(8),dimension(ndi)::x,xc
    real(8),dimension(mpol)::polyn
    real(8),dimension(ndi,mpol)::dpolyn
    real(8),dimension(ndi,ndi,mpol)::dpolyn2
    real(8),dimension(mpol,ntot)::u2
!-----------------------
!*** gets closest nodes
!-----------------------
    call mls_getsnodes(rmax,x,ndi,ntot,xtot,n,listn)
    do i=1,n
       do id=1,ndi
          xn(id,i)=xtot(id,listn(i))
       end do
    end do
!-----------------------------------
!*** check for repeated coordinates
!-----------------------------------
!!$    do in=1,n
!!$       do jn=in+1,n
!!$          if(rnorm2(ndi,xn(1:ndi,in)-xn(1:ndi,jn)).le.1.0d-20)then
!!$             write(*,*)"node in and jn",in,jn," are the same"
!!$             write(*,*)"in,x=",in,xn(1:ndi,in)
!!$             write(*,*)"jn,x=",jn,xn(1:ndi,jn)
!!$             read(*,*)
!!$          end if
!!$       end do
!!$    end do
!-------------------------------
!*** shape function derivatives
!-------------------------------
    if(imeshless.ge.0)then
       call rbf_shapefunctions(imeshless,rmax,ndi,n,xn,x,xc,ffloc,dffloc,dff2loc)
    else
!---------  
!**** mls
!---------
       tol=1.0d-2
       call mls_u2atacoordinate(imeshless,rmax,tol,ndi,n,xc,xn,x,u2(1:mpol,1:n))
       select case(imeshless)
       case(-1)
          call mls_polynlinbase(ndi,rmax,x(1:ndi),xc(1:ndi),polyn(1:mpol),dpolyn(1:ndi,1:mpol),dpolyn2(1:ndi,1:ndi,1:mpol),m)
       case(-2)
          call mls_polynqbase(ndi,rmax,x(1:ndi),xc(1:ndi),polyn(1:mpol),dpolyn(1:ndi,1:mpol),dpolyn2(1:ndi,1:ndi,1:mpol),m)
       case(-3)
          CALL mls_polyncubbase(ndi,rmax,x(1:ndi),xc(1:ndi),polyn(1:mpol),dpolyn(1:ndi,1:mpol),dpolyn2(1:ndi,1:ndi,1:mpol),m)
       end select
       call mls_sfder(ndi,n,m,polyn(1:m),dpolyn(1:ndi,1:m),dpolyn2(1:ndi,1:ndi,1:m),u2(1:m,1:n),ffloc,dffloc,dff2loc)
    end if
    ff=0.0d00
    dff=0.0d00
    dff2=0.0d00
    do i=1,n
       ff(listn(i))=ffloc(i)
       dff(1:ndi,listn(i))=dffloc(1:ndi,i)
       dff2(1:ndi,1:ndi,listn(i))=dff2loc(1:ndi,1:ndi,i)
    end do
    deallocate(listn)
  END SUBROUTINE rbf_testingallrequired
!--------------------------
!*** END RBF, END MESHLESS
!--------------------------

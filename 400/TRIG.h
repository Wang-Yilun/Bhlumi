*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                     Pseudo-CLASS  TRIG                                   //
*//                                                                          //
*//   Purpose: Emulate LEP event selections for LEP luminosity detectors     //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      DOUBLE PRECISION      m_pi
      PARAMETER( m_pi=3.1415926535897932d0)

      INTEGER              NPHI,NTHE,Npad,Nseg   ! segmentation and cluster size
      DOUBLE PRECISION     TMIND,TMAXD           ! Fidutial angular range
      DOUBLE PRECISION     th1w,th2w, th1n,th2n  ! Wide and narrow ang acceptance
      DOUBLE PRECISION     Acopl,Acoli           ! Acoplanarity and acolinearity cuts
      DOUBLE PRECISION     EcutMin,EcutMax,EcutSum,EcutProd ! Energy cut params

      COMMON /c_TRIG/
     $ TMIND,TMAXD,NPHI,NTHE,Npad,Nseg,    ! Fidutial angular range
     $ th1w,th2w, th1n,th2n,               ! Wide and narrow ang. acceptance
     $ Acopl,Acoli,                        ! Acoplanarity and acolinearity cuts
     $ EcutMin,EcutMax,EcutSum,EcutProd    ! Energy cut params

***      SAVE /c_TRIG/   ! not needed

*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of Pseudo-CLASS  TRIG                           //
*//////////////////////////////////////////////////////////////////////////////

      implicit real*8(a-h, o-z)
	integer SEED,MX,MY,MZ,NA,LIST,NSP1,NSP2
      integer NM,NMZ,MAPSZE,NCELL,MAPS,HEAD
      integer molstart,molend,molrange
      double precision r,axes,q,CL,transstep,rotstep,volstep
        double precision HN,HNI,VOLN,a8,dia,dh,dminsq,dmax2
      double precision dp,dpsq
	double precision rcutsq11,rcut,dgap,vertdh
	double precision roversq11
	double precision press,temp,beta,cellix,sqrt3by2
      double precision celliy,celliz,daxis
      integer itype
      PARAMETER (sqrt3by2=3.0d0**0.5d0/2.0d0)
      PARAMETER (a6 = 2.0d0,rcutsq11=3.0d0*a6**2,
     >           roversq11=a6**2,roverh11 = 0.50d0*a6) !is roverh the L1

       PARAMETER(NA=1728)   ! for cubes + spheres (12x12x12)
       PARAMETER( NM=1600,NMZ=26*NM) ! PARAMETER( NM=5832,NMZ=26*NM) 

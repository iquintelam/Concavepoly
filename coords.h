        common /coords/ r(3,NA),axes(3,3,NA),q(3,NA),
     >                 itype(NA),NSP1,NSP2
	common /boxsize/ CL(3),cellix,celliy,celliz,MX,MY,MZ
	common /movsizes/ transstep(3),rotstep(3),volstep(3,3)
	common /var/ press,temp,beta,SEED,LIST(NA)
	common /map/ MAPS(NMZ),HEAD(NM),NCELL,MAPSZE
        common /tensor/ HN(3,3),HNI(3,3),VOLN
 	common /sle/ molstart,molend,molrange
        common /vectorij/qij(3),tt(3)
        common /forcefield/ene12,rcutsq12ene,rcutsq12e,cos2d2,deltamax
        common /limits/rcut,rcuti,rcut2(2,2),rover2(2,2),iboxshape
        common /spheres/dia,rcutsq22,rcutsq12,roversq22,roversq12
        common /alchemi/deltad,rcut2new(2,2,2),rover2new(2,2,2),
     >       rcutsq12enew(2),roversqh22new(0:2),roversqh22
        common /identation/dh,dhsq,dminsq,dmax2,dminsqnew(0:2),
     >       dhsqnew(0:2),dhnew(0:2),daxis(8,3),dp,dpsq,
     >       dpn(0:2),dpsqn(0:2)

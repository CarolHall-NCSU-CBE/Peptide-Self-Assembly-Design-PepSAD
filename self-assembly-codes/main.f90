module constant

    implicit none
! The Data type "groupdetails" is used to define an array to store the information on the residues and atoms.
	type groupdetails
		integer				:: cnum1, cnum2, cnum3
		character*4			:: atype1(20), atype2(60), atype3(20)
		character*4			:: gtype
		real				:: coo1(20,3), coo2(60,3), coo3(20,3)
	end type

	type energyparameters
		integer				:: iac1(20), iac2(60), iac3(20), atomid1(20), atomid2(60), atomid3(20)
		real				:: charge1(20), epsion1(20), r1(20), rborn1(20), fs1(20), dielecons1(20)
		real				:: charge2(60), epsion2(60), r2(60), rborn2(60), fs2(60), dielecons2(60)
		real				:: charge3(20), epsion3(20), r3(20), rborn3(20), fs3(20), dielecons3(20)
	end type
	
	type lib4aminoacid
		integer				:: cnum1, cnum2, cnum3
		integer				:: dihedralangle(34,4)		
		character*4			:: atype1(20), atype2(60), atype3(20)
		character*4			:: gtype
		real				:: coo1(20,3), coo2(60,3), coo3(20,3)
		integer				:: grade, rotanum
	end type
	
	type index4sidechain
		integer				:: class_No, member_No
	end type

	type conformer4sidechain
		real				:: member(15,3)
	end type

	type dihedralparameters
		integer*8 			:: iph(36), jph(36), kph(36), lph(36), multiply(36)
		real*8    			:: pk(36,4), pn(36,4), phase(36,4)
	end type

	type peptideassembly
		integer			    :: num4peptides, peptideID(6)
	end type

	integer, parameter		:: gnum=14
	integer, parameter		:: repeated_unit=8
	integer, parameter		:: num4category=2
	integer, parameter		:: atom_num=gnum*60
	integer, parameter		:: num4pdbsave=10
	integer					:: nstep_start, nstep_terminal, idum, recalcu_switch
	integer					:: fpho,fpol,fchg,foth

	real, parameter			:: scmfswitch1=0.6
!	real, parameter			:: scmfswitch2=0.8
	real, parameter			:: dihedral_weighting_factor=0.20
	real, parameter			:: propensity_weighting_factor=3.0
	real, parameter			:: vdw14_coeff=2.0  !1-4 VDW scale factor
	real, parameter			:: ele14_coeff=1.2  !1-4 ELE scale factor

	real					:: ekt
	real					:: energy_min(num4pdbsave)

	type RES4chain
		integer				:: num
		integer				:: IDs(gnum)
	end type

	type(lib4aminoacid)		:: aa_lib(20)
	
	type(peptideassembly)	:: selfassembly(num4category)
	character*20			:: filename

end module constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module randomgenerator

	use constant

	contains
! Generate the random number
	subroutine ran_gen(ran2,flag_ran)
	implicit none
    integer			:: im1,im2,imm2,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv,flag_ran,start_flag
    real			:: ran2,am,eps,rnmx
    parameter       (im1=2147483563,im2=2147483399,am=1./im1,imm2=im1-1,ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211, &
					ir2=3791,ntab=32,ndiv=1+imm2/ntab,eps=1.2e-7,rnmx=1.-eps)
	integer			:: idum2,j,k,iv(ntab),iy
    save			iv,iy,idum2
	save			start_flag
    data            idum2/123456789/, iv/ntab*0/, iy/0/, start_flag/0/

	if(recalcu_switch.ne.0.and.start_flag.eq.0) then
		start_flag=1
		open(20, file="randomnumber.txt", status="old")
			read(20,*)
			read(20,*) idum
			read(20,*)
			read(20,*) idum2
			read(20,*)
			do j=1, ntab
				read(20,*) iv(j)
			enddo
			read(20,*)
			read(20,*) iy
		close(20)		
	endif

	if(flag_ran==0) then
		if (idum.le.0) then
			idum=max(-idum,1)
			idum2=idum
			do 100 j=ntab+8,1,-1
				k=idum/iq1
				idum=ia1*(idum-k*iq1)-k*ir1
				if (idum.lt.0) idum=idum+im1
				if (j.le.ntab) iv(j)=idum
100			continue
			iy=iv(1)
		endif
		k=idum/iq1  
		idum=ia1*(idum-k*iq1)-k*ir1	  
		if (idum.lt.0) idum=idum+im1
		k=idum2/iq2 
		idum2=ia2*(idum2-k*iq2)-k*ir2
		if (idum2.lt.0) idum2=idum2+im2  
		j=1+iy/ndiv 
		iy=iv(j)-idum2	  
		iv(j)=idum  
		if(iy.lt.1) iy=iy+imm2  
		ran2=min(am*iy,rnmx)
	else
		open(20, file="randomnumber.txt", status="replace")
			write(20,*) "idum="
			write(20,*) idum
			write(20,*) "idum2="
			write(20,*) idum2
			write(20,*) "iv(ntab=32)="
			do j=1, ntab
				write(20,*) iv(j)
			enddo
			write(20,*) "iy="
			write(20,*) iy
		close(20)
	endif		

    return
    end subroutine ran_gen

end module randomgenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module input

	use constant

	contains
	subroutine inputfile
	implicit none
	integer						:: i,j
	
	open(10, file="input.txt", status="old")
		read(10,*) filename, recalcu_switch
		read(10,*) nstep_start, nstep_terminal
		read(10,*) idum, ekt
		read(10,*) fpho,fpol,fchg,foth
		do i=1, num4category
			read(10,*) selfassembly(i)%num4peptides
			read(10,*) (selfassembly(i)%peptideID(j),j=1,selfassembly(i)%num4peptides)
		enddo
	close(10)

	return
	end subroutine inputfile

end module input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdbfile

	use	constant

	contains
! Read an initial file.
	subroutine readpdb(group)
	implicit none
	integer						:: anum, status, chainID, ic
	real						:: x, y, z
	character*4					:: char, atype, name
	type(groupdetails)			:: group(repeated_unit,gnum)

	group%cnum1=0
	group%cnum2=0
	group%cnum3=0
	chainID=1

	open(10, file=filename)
	do while(.true.)
		read(10, *, iostat=status) char, anum, atype, name, ic, x, y, z
		if(status.ne.0) exit
		if(ic.gt.(chainID*gnum)) chainID=chainID+1
		ic=ic-(chainID-1)*gnum
		group(chainID,ic)%gtype=name
		if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
		   .or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then  
			group(chainID,ic)%cnum1=group(chainID,ic)%cnum1+1
			group(chainID,ic)%atype1(group(chainID,ic)%cnum1)=atype
			group(chainID,ic)%coo1(group(chainID,ic)%cnum1,1)=x
			group(chainID,ic)%coo1(group(chainID,ic)%cnum1,2)=y
			group(chainID,ic)%coo1(group(chainID,ic)%cnum1,3)=z
		elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
			group(chainID,ic)%cnum3=group(chainID,ic)%cnum3+1
			group(chainID,ic)%atype3(group(chainID,ic)%cnum3)=atype
			group(chainID,ic)%coo3(group(chainID,ic)%cnum3,1)=x
			group(chainID,ic)%coo3(group(chainID,ic)%cnum3,2)=y
			group(chainID,ic)%coo3(group(chainID,ic)%cnum3,3)=z
		else
			group(chainID,ic)%cnum2=group(chainID,ic)%cnum2+1
			group(chainID,ic)%atype2(group(chainID,ic)%cnum2)=atype
			group(chainID,ic)%coo2(group(chainID,ic)%cnum2,1)=x
			group(chainID,ic)%coo2(group(chainID,ic)%cnum2,2)=y
			group(chainID,ic)%coo2(group(chainID,ic)%cnum2,3)=z
		endif	
	enddo
	close(10)
	
	return
	end subroutine readpdb


	subroutine generatepdb(step, attempt, group)
	implicit none
	integer							:: step, attempt
	integer							:: i, j, k, atomnum
	character*5						:: stepchar, attemptchar
	type(groupdetails)				:: group(repeated_unit,gnum)

	write(stepchar,"(i5)") step
	write(attemptchar,"(i4)") attempt
	open(10,file='pdbfiles/'//trim(adjustl(stepchar))//'_'//trim(adjustl(attemptchar))//'.pdb', access="append")
		atomnum=1
		do i=1, repeated_unit
			do j=1, gnum
				do k=1, group(i,j)%cnum1
					if(len_trim(group(i,j)%atype1(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype1(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype1(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum2
					if(len_trim(group(i,j)%atype2(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype2(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype2(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum3
					if(len_trim(group(i,j)%atype3(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype3(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype3(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					endif
					atomnum=atomnum+1
				enddo
			enddo
		enddo
1		format(a4,i7,a6,a4,i5,f12.3,2f8.3,2f6.2,a16)
2		format(a4,i7,a5,a1,a4,i5,f12.3,2f8.3,2f6.2,a16)
	close(10)

	return
	end subroutine generatepdb

end module pdbfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mathfunction

	use constant

	contains
	subroutine normalvector(rsta, rmid, rend, r_nor)
	implicit none
	real						:: rsta(3), rmid(3), rend(3), r_nor(3)
	real						:: a(3), b(3)

	a=rsta-rmid
	b=rend-rmid

	r_nor(1)=a(2)*b(3)-a(3)*b(2)
	r_nor(2)=a(3)*b(1)-a(1)*b(3)
	r_nor(3)=a(1)*b(2)-a(2)*b(1)

	return
	end subroutine normalvector


	subroutine vectorrotation(rsta, rend, m)
	implicit none
	real						:: rsta(3), rend(3), m(3,3)
	real						:: r_cropro(3), a(3), a1(3,3), a2(3,3), a3(3,3)
	real						:: absrsta, absrend, r_dotpro, cos, sin

	absrsta=sqrt(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
	absrend=sqrt(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))

	r_dotpro=dot_product(rsta, rend)

	r_cropro(1)=rsta(2)*rend(3)-rsta(3)*rend(2)
	r_cropro(2)=rsta(3)*rend(1)-rsta(1)*rend(3)
	r_cropro(3)=rsta(1)*rend(2)-rsta(2)*rend(1)

	cos=r_dotpro/(absrsta*absrend)
	sin=sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))/(absrsta*absrend)

	a(1)=r_cropro(1)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
	a(2)=r_cropro(2)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
	a(3)=r_cropro(3)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))

	a1(1,1)=a(1)*a(1)
	a1(1,2)=a(1)*a(2)
	a1(1,3)=a(1)*a(3)
	a1(2,1)=a(2)*a(1)
	a1(2,2)=a(2)*a(2)
	a1(2,3)=a(2)*a(3)
	a1(3,1)=a(3)*a(1)
	a1(3,2)=a(3)*a(2)
	a1(3,3)=a(3)*a(3)

	a2=-a1
	a2(1,1)=1+a2(1,1)
	a2(2,2)=1+a2(2,2)
	a2(3,3)=1+a2(3,3)
	a2=cos*a2

	a3(1,1)=0.0
	a3(1,2)=-a(3)
	a3(1,3)=a(2)
	a3(2,1)=a(3)
	a3(2,2)=0.0
	a3(2,3)=-a(1)
	a3(3,1)=-a(2)
	a3(3,2)=a(1)
	a3(3,3)=0
	a3=sin*a3

	m=a1+a2+a3
	m=transpose(m)

	return	
	end subroutine vectorrotation

	
	subroutine axisrotation(a, cos, sin, m)
	implicit none
	real					:: cos, sin
	real					:: a(3), a1(3,3), a2(3,3), a3(3,3)
	real					:: m(3,3)

	a1(1,1)=a(1)*a(1)
	a1(1,2)=a(1)*a(2)
	a1(1,3)=a(1)*a(3)
	a1(2,1)=a(2)*a(1)
	a1(2,2)=a(2)*a(2)
	a1(2,3)=a(2)*a(3)
	a1(3,1)=a(3)*a(1)
	a1(3,2)=a(3)*a(2)
	a1(3,3)=a(3)*a(3)

	a2=-a1
	a2(1,1)=1+a2(1,1)
	a2(2,2)=1+a2(2,2)
	a2(3,3)=1+a2(3,3)
	a2=cos*a2

	a3(1,1)=0.0
	a3(1,2)=-a(3)
	a3(1,3)=a(2)
	a3(2,1)=a(3)
	a3(2,2)=0.0
	a3(2,3)=-a(1)
	a3(3,1)=-a(2)
	a3(3,2)=a(1)
	a3(3,3)=0
	a3=sin*a3

	m=a1+a2+a3
	m=transpose(m)

	return	
	end subroutine axisrotation

	
	subroutine phipsiomg_angle(p1, p2, p3, p4, angle)
	implicit none
	real						:: p1(3), p2(3), p3(3), p4(3)
	real						:: angle, angle_T1, angle_T2
	real						:: rsta(3), rend(3), rmid(3)
	real						:: r_1(3), r_2(3)
	real						:: absrsta, absrend, r_dotpro, value

	call normalvector(p1, p2, p3, rend)
	call normalvector(p2, p3, p4, rsta)

	absrsta=sqrt(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
	absrend=sqrt(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))

	r_dotpro=dot_product(rsta, rend)
	value=r_dotpro/(absrsta*absrend)
	if(value.lt.-1.00) then
		value=-1.00
	elseif(value.gt.1.00) then	
		value=1.00
	endif
	angle_T1=acosd(value)

	if(abs(180.0-angle_T1).le.(0.5)) then
		angle=180.0
	elseif(abs(angle_T1).le.(0.4)) then
		angle=0.0
	else
		rmid=0.0
		r_2=p3-p2
		call normalvector(rsta, rmid, rend, r_1)

		absrsta=sqrt(r_1(1)*r_1(1)+r_1(2)*r_1(2)+r_1(3)*r_1(3))
		absrend=sqrt(r_2(1)*r_2(1)+r_2(2)*r_2(2)+r_2(3)*r_2(3))

		r_dotpro=dot_product(r_1, r_2)
		if(abs(r_dotpro/(absrsta*absrend)-1.0).le.(0.1)) then
			angle_T2=0.00
		elseif(abs(r_dotpro/(absrsta*absrend)+1.0).le.(0.1)) then
			angle_T2=180
		else
			value=r_dotpro/(absrsta*absrend)
			if(value.lt.-1.00) then
				value=-1.00
			elseif(value.gt.1.00) then	
				value=1.00
			endif	
			angle_T2=acosd(value)
		endif
		if(angle_T2.gt.90) then
			angle=-angle_T1
		else
			angle=angle_T1
		endif
	endif

	return
	end subroutine phipsiomg_angle

	
	subroutine variance_covariance_matrix(matrix,obs,n,det)
	Implicit none

	integer		     					:: obs, n, i, j, k
	real								:: matrix(34,4),a(n,n),avg(n),matrix_dif(obs,n)
	real								:: det

	do i=1,obs
		do j=1, n
			if(matrix(i,j).gt.180) then
					matrix(i,j)=matrix(i,j)-360.0
			elseif(matrix(i,j).lt.-180) then
					matrix(i,j)=matrix(i,j)+360.0
			endif
		enddo
	enddo
	matrix=matrix*0.0174533

	avg=0.0
	do j=1, n
		do i=1, obs
			avg(j)=avg(j)+matrix(i,j)
		enddo
	enddo
	avg=avg/obs

	do j=1, n
		do i=1, obs
			matrix_dif(i,j)=matrix(i,j)-avg(j)
		enddo
	enddo

	a=0.0
	do j=1, n
		do k=1, j
			do i=1, obs
				a(j,k)=a(j,k)+matrix_dif(i,j)*matrix_dif(i,k)
			enddo
			a(j,k)=a(j,k)/(obs-1)
			a(k,j)=a(j,k)
		enddo
	enddo

	call eigenfunction(n,a,det)

	return
	end subroutine variance_covariance_matrix


	subroutine eigenfunction(n,a,det)
	implicit none
	
	integer		:: n,i,j,p,q
	real		:: amax,temp,zemp,coo,sii,co,si,app,aqq,apq,api,aqi
	real		:: rip,riq,det
	real		:: a(n,n),r(n,n)
	character	:: name*12

	do i=1,n
		do j=1,n
			r(i,j)=0.0
		enddo
		r(i,i)=1.0
	enddo

10	amax=abs(a(2,1))
	p=2
	q=1
	do i=2,n
		do j=1,i-1
			if(abs(a(i,j)).gt.amax) then
				amax=abs(a(i,j))
	 			p=i
				q=j
			endif
		enddo
	enddo

	if(amax.le.1.0e-7) then
		goto 20
	endif

	temp=2*a(p,q)/(a(p,p)-a(q,q)+1.0e-30)
	zemp=(a(p,p)-a(q,q))/(2*a(p,q))
	
	if(abs(temp).lt.1.0) then
		coo=(1+temp**2)**(-0.5)
		sii=temp*(1+temp**2)**(-0.5)
	else
		coo=abs(zemp)*(1+zemp**2)**(-0.5)
		sii=sign(1.0,zemp)*(1+zemp**2)**(-0.5)
	endif
	
	co=sqrt(0.5*(1+coo))
	si=sii/(2.0*co)

	do i=1,n
		rip=r(i,p)*co+r(i,q)*si
		riq=-r(i,p)*si+r(i,q)*co
		r(i,p)=rip
		r(i,q)=riq
	enddo

	app=a(p,p)*co**2+a(q,q)*si**2+2*a(p,q)*co*si
	aqq=a(p,p)*si**2+a(q,q)*co**2-2*a(p,q)*co*si
	apq=0.5*(a(q,q)-a(p,p))*sii+a(p,q)*coo
	a(p,p)=app
	a(q,q)=aqq
	a(p,q)=apq
	a(q,p)=a(p,q)
	
	do i=1,n
		if(i.eq.p.or.i.eq.q) then
		else
			api=a(p,i)*co+a(q,i)*si
			aqi=-a(p,i)*si+a(q,i)*co
			a(p,i)=api
			a(q,i)=aqi
			a(i,p)=a(p,i)
			a(i,q)=a(q,i)
		endif
	enddo	
    goto 10

20	continue
	det=a(1,1)
	do i=2,n
		det=det*a(i,i)
	enddo

	return
	end	subroutine eigenfunction

end module mathfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module database

	use constant
	use randomgenerator
	use mathfunction

	contains
	subroutine rotamerlib
	implicit none
	integer							:: status, grade, rotanum, anum, num, i, j, k
	real							:: x, y, z
	character*4						:: char, atype, name

	open(10, file="lib/rotamer", status="old")
		do while(.true.)
			read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum		
			if(status.ne.0) exit
			i=0
			if(name=="GLY") then
				i=1
			elseif(name=="LEU") then
				i=2
			elseif(name=="VAL") then
				i=3
			elseif(name=="ILE") then
				i=4
			elseif(name=="MET") then
				i=5
			elseif(name=="PHE") then
				i=6
			elseif(name=="TYR") then
				i=7
			elseif(name=="TRP") then
				i=8
			elseif(name=="ARG") then
				i=9
			elseif(name=="LYS") then
				i=10
			elseif(name=="SER") then
				i=11
			elseif(name=="THR") then
				i=12
			elseif(name=="ASN") then
				i=13
			elseif(name=="GLN") then
				i=14
			elseif(name=="HIE") then
				i=15
			elseif(name=="PRO") then
				i=16
			elseif(name=="CYS") then
				i=17
			elseif(name=="ALA") then
				i=18
			elseif(name=="GLU") then
				i=19
			elseif(name=="ASP") then
				i=20
			endif
				
			if(i.ne.0) then
				aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
				if(grade.ne.0) then
					do j=1, rotanum
						read(10,"(<grade>(i5))") (aa_lib(i)%dihedralangle(j,k), k=1, grade)
					enddo
				endif
			endif
		enddo
	close(10)

	aa_lib%cnum1=0
	aa_lib%cnum2=0
	aa_lib%cnum3=0

	do i=1, 20
		open (10, file='lib/RotamerLibrary/'//trim(aa_lib(i)%gtype), status="old")	
		do while(.true.)
			read(10, *, iostat=status) char, anum, atype, name, char, num, x, y, z
			if(status.ne.0) exit
			if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
			   .or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then  
				aa_lib(i)%cnum1=aa_lib(i)%cnum1+1
				aa_lib(i)%atype1(aa_lib(i)%cnum1)=atype
				aa_lib(i)%coo1(aa_lib(i)%cnum1,1)=x
				aa_lib(i)%coo1(aa_lib(i)%cnum1,2)=y
				aa_lib(i)%coo1(aa_lib(i)%cnum1,3)=z
			elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
				aa_lib(i)%cnum3=aa_lib(i)%cnum3+1
				aa_lib(i)%atype3(aa_lib(i)%cnum3)=atype
				aa_lib(i)%coo3(aa_lib(i)%cnum3,1)=x
				aa_lib(i)%coo3(aa_lib(i)%cnum3,2)=y
				aa_lib(i)%coo3(aa_lib(i)%cnum3,3)=z
			else
				aa_lib(i)%cnum2=aa_lib(i)%cnum2+1
				aa_lib(i)%atype2(aa_lib(i)%cnum2)=atype
				aa_lib(i)%coo2(aa_lib(i)%cnum2,1)=x
				aa_lib(i)%coo2(aa_lib(i)%cnum2,2)=y
				aa_lib(i)%coo2(aa_lib(i)%cnum2,3)=z
			endif
		enddo
		close(10)
	enddo
	
	return
	end subroutine rotamerlib


	subroutine findrotamer(ic, group, name_original, rotanum, aa_group, ip)
	implicit none
	integer								:: categoryID, status, ic, rotanum, i, ii, j, k, l, n, ip
	integer								:: grade, grade_num(6), monitor(6)
	real								:: nr(3), car(3), cr(3), r_norpep(3)
	real								:: aa_nr(3), aa_car(3), aa_cr(3), r_norrot(3)
	real								:: r_nca(3), aa_r_nca(3), r_trans(3)
	real								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	real								:: delta_chi, cos_angle, sin_angle	
	real								:: temp1(20,3), temp2(60,3), temp3(20,3)
	character*4							:: name_original

	type(groupdetails)					:: group(repeated_unit,gnum), aa_group(repeated_unit,40)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6)

	if(name_original=="GLY".or.name_original=="NGLY".or.name_original=="CGLY") then
		ip=1
	elseif(name_original=="LEU".or.name_original=="NLEU".or.name_original=="CLEU") then
		ip=2
	elseif(name_original=="VAL".or.name_original=="NVAL".or.name_original=="CVAL") then
		ip=3
	elseif(name_original=="ILE".or.name_original=="NILE".or.name_original=="CILE") then
		ip=4
	elseif(name_original=="MET".or.name_original=="NMET".or.name_original=="CMET") then
		ip=5
	elseif(name_original=="PHE".or.name_original=="NPHE".or.name_original=="CPHE") then
		ip=6
	elseif(name_original=="TYR".or.name_original=="NTYR".or.name_original=="CTYR".or. &
	       name_original=="TYX".or.name_original=="NTYX".or.name_original=="CTYX") then
		ip=7
	elseif(name_original=="TRP".or.name_original=="NTRP".or.name_original=="CTRP") then
		ip=8
	elseif(name_original=="ARG".or.name_original=="NARG".or.name_original=="CARG".or. &
	       name_original=="ARN".or.name_original=="NARN".or.name_original=="CARN") then
		ip=9
	elseif(name_original=="LYS".or.name_original=="NLYS".or.name_original=="CLYS".or. &
	       name_original=="LYN".or.name_original=="NLYN".or.name_original=="CLYN") then
		ip=10
	elseif(name_original=="SER".or.name_original=="NSER".or.name_original=="CSER") then
		ip=11
	elseif(name_original=="THR".or.name_original=="NTHR".or.name_original=="CTHR") then
		ip=12
	elseif(name_original=="ASN".or.name_original=="NASN".or.name_original=="CASN") then
		ip=13
	elseif(name_original=="GLN".or.name_original=="NGLN".or.name_original=="CGLN") then
		ip=14
	elseif(name_original=="HIE".or.name_original=="NHIE".or.name_original=="CHIE".or. &
	       name_original=="HIP".or.name_original=="NHIP".or.name_original=="CHIP") then
		ip=15
	elseif(name_original=="PRO".or.name_original=="NPRO".or.name_original=="CPRO") then
		ip=16
	elseif(name_original=="CYS".or.name_original=="NCYS".or.name_original=="CCYS".or. &
	       name_original=="CYT".or.name_original=="NCYT".or.name_original=="CCYT") then
		ip=17
	elseif(name_original=="ALA".or.name_original=="NALA".or.name_original=="CALA") then
		ip=18
	elseif(name_original=="GLU".or.name_original=="NGLU".or.name_original=="CGLU".or. &
	       name_original=="GLH".or.name_original=="NGLH".or.name_original=="CGLH") then
		ip=19
	elseif(name_original=="ASP".or.name_original=="NASP".or.name_original=="CASP".or. &
	       name_original=="ASH".or.name_original=="NASH".or.name_original=="CASH") then
		ip=20
	endif
	
	do ii=1, repeated_unit	
		rotanum=aa_lib(ip)%rotanum
		aa_group(ii,1)%cnum1=aa_lib(ip)%cnum1; aa_group(ii,1)%cnum2=aa_lib(ip)%cnum2; aa_group(ii,1)%cnum3=aa_lib(ip)%cnum3
		aa_group(ii,1)%gtype=name_original
		aa_group(ii,1)%atype1=aa_lib(ip)%atype1; aa_group(ii,1)%atype2=aa_lib(ip)%atype2; aa_group(ii,1)%atype3=aa_lib(ip)%atype3
		aa_group(ii,1)%coo1=aa_lib(ip)%coo1; aa_group(ii,1)%coo2=aa_lib(ip)%coo2; aa_group(ii,1)%coo3=aa_lib(ip)%coo3
	
		do k=1, group(ii,ic)%cnum1
			if(group(ii,ic)%atype1(k)=="N") then
				nr(1)=group(ii,ic)%coo1(k,1)
				nr(2)=group(ii,ic)%coo1(k,2)
				nr(3)=group(ii,ic)%coo1(k,3)
			elseif(group(ii,ic)%atype1(k)=="CA") then
				car(1)=group(ii,ic)%coo1(k,1)
				car(2)=group(ii,ic)%coo1(k,2)
				car(3)=group(ii,ic)%coo1(k,3)
			endif
		enddo
		do k=1, group(ii,ic)%cnum3
			if(group(ii,ic)%atype3(k)=="C") then
				cr(1)=group(ii,ic)%coo3(k,1)
				cr(2)=group(ii,ic)%coo3(k,2)
				cr(3)=group(ii,ic)%coo3(k,3)
			endif
		enddo

		call normalvector(nr, car, cr, r_norpep)

		r_nca(1)=nr(1)-car(1)
		r_nca(2)=nr(2)-car(2)
		r_nca(3)=nr(3)-car(3)

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="N") then
				aa_nr(1)=aa_group(ii,1)%coo1(k,1)
				aa_nr(2)=aa_group(ii,1)%coo1(k,2)
				aa_nr(3)=aa_group(ii,1)%coo1(k,3)
			elseif(aa_group(ii,1)%atype1(k)=="CA") then
				aa_car(1)=aa_group(ii,1)%coo1(k,1)
				aa_car(2)=aa_group(ii,1)%coo1(k,2)
				aa_car(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo
		do k=1, aa_group(ii,1)%cnum3
			if(aa_group(ii,1)%atype3(k)=="C") then
				aa_cr(1)=aa_group(ii,1)%coo3(k,1)
				aa_cr(2)=aa_group(ii,1)%coo3(k,2)
				aa_cr(3)=aa_group(ii,1)%coo3(k,3)
			endif
		enddo	
		
		call normalvector(aa_nr, aa_car, aa_cr, r_norrot)

		call vectorrotation(r_norrot, r_norpep, m)
		
		temp1=matmul(aa_group(ii,1)%coo1, m)
		aa_group(ii,1)%coo1=temp1

		temp2=matmul(aa_group(ii,1)%coo2, m)
		aa_group(ii,1)%coo2=temp2

		temp3=matmul(aa_group(ii,1)%coo3, m)
		aa_group(ii,1)%coo3=temp3

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="N") then
				aa_nr(1)=aa_group(ii,1)%coo1(k,1)
				aa_nr(2)=aa_group(ii,1)%coo1(k,2)
				aa_nr(3)=aa_group(ii,1)%coo1(k,3)
			elseif(aa_group(ii,1)%atype1(k)=="CA") then
				aa_car(1)=aa_group(ii,1)%coo1(k,1)
				aa_car(2)=aa_group(ii,1)%coo1(k,2)
				aa_car(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo

		aa_r_nca(1)=aa_nr(1)-aa_car(1)
		aa_r_nca(2)=aa_nr(2)-aa_car(2)
		aa_r_nca(3)=aa_nr(3)-aa_car(3)

		call vectorrotation(aa_r_nca, r_nca, m)

		temp1=matmul(aa_group(ii,1)%coo1, m)
		aa_group(ii,1)%coo1=temp1

		temp2=matmul(aa_group(ii,1)%coo2, m)
		aa_group(ii,1)%coo2=temp2

		temp3=matmul(aa_group(ii,1)%coo3, m)
		aa_group(ii,1)%coo3=temp3

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="CA") then
				aa_car(1)=aa_group(ii,1)%coo1(k,1)
				aa_car(2)=aa_group(ii,1)%coo1(k,2)
				aa_car(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo

		r_trans(1)=car(1)-aa_car(1)
		r_trans(2)=car(2)-aa_car(2)
		r_trans(3)=car(3)-aa_car(3)

		do k=1, aa_group(ii,1)%cnum1
			aa_group(ii,1)%coo1(k,1)=anint((aa_group(ii,1)%coo1(k,1)+r_trans(1))*1000)/1000
			aa_group(ii,1)%coo1(k,2)=anint((aa_group(ii,1)%coo1(k,2)+r_trans(2))*1000)/1000
			aa_group(ii,1)%coo1(k,3)=anint((aa_group(ii,1)%coo1(k,3)+r_trans(3))*1000)/1000
		enddo
		do k=1, aa_group(ii,1)%cnum2
			aa_group(ii,1)%coo2(k,1)=anint((aa_group(ii,1)%coo2(k,1)+r_trans(1))*1000)/1000
			aa_group(ii,1)%coo2(k,2)=anint((aa_group(ii,1)%coo2(k,2)+r_trans(2))*1000)/1000
			aa_group(ii,1)%coo2(k,3)=anint((aa_group(ii,1)%coo2(k,3)+r_trans(3))*1000)/1000
		enddo
		do k=1, aa_group(ii,1)%cnum3
			aa_group(ii,1)%coo3(k,1)=anint((aa_group(ii,1)%coo3(k,1)+r_trans(1))*1000)/1000
			aa_group(ii,1)%coo3(k,2)=anint((aa_group(ii,1)%coo3(k,2)+r_trans(2))*1000)/1000
			aa_group(ii,1)%coo3(k,3)=anint((aa_group(ii,1)%coo3(k,3)+r_trans(3))*1000)/1000
		enddo	
	
		if(rotanum.le.1) goto 10

		grade_num=0
		if(aa_group(ii,1)%gtype=="VAL".or.aa_group(ii,1)%gtype=="NVAL".or.aa_group(ii,1)%gtype=="CVAL") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="LEU".or.aa_group(ii,1)%gtype=="NLEU".or.aa_group(ii,1)%gtype=="CLEU") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else	
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo		
		elseif(aa_group(ii,1)%gtype=="ILE".or.aa_group(ii,1)%gtype=="NILE".or.aa_group(ii,1)%gtype=="CILE") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2	
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3) 
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB".or.aa_group(ii,1)%atype2(k)=="CG2".or.aa_group(ii,1)%atype2(k)=="HG21".or.aa_group(ii,1)%atype2(k)=="HG22".or. &
					aa_group(ii,1)%atype2(k)=="HG23".or.aa_group(ii,1)%atype2(k)=="CG1") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG1") monitor(2)=grade_num(2)
				else	
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="PHE".or.aa_group(ii,1)%gtype=="NPHE".or.aa_group(ii,1)%gtype=="CPHE") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="TRP".or.aa_group(ii,1)%gtype=="NTRP".or.aa_group(ii,1)%gtype=="CTRP") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else	
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="TYR".or.aa_group(ii,1)%gtype=="NTYR".or.aa_group(ii,1)%gtype=="CTYR".or. &
			   aa_group(ii,1)%gtype=="TYX".or.aa_group(ii,1)%gtype=="NTYX".or.aa_group(ii,1)%gtype=="CTYX") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo						
		elseif(aa_group(ii,1)%gtype=="SER".or.aa_group(ii,1)%gtype=="NSER".or.aa_group(ii,1)%gtype=="CSER") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo		
		elseif(aa_group(ii,1)%gtype=="THR".or.aa_group(ii,1)%gtype=="NTHR".or.aa_group(ii,1)%gtype=="CTHR") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="CYS".or.aa_group(ii,1)%gtype=="NCYS".or.aa_group(ii,1)%gtype=="CCYS".or. &
			   aa_group(ii,1)%gtype=="CYT".or.aa_group(ii,1)%gtype=="NCYT".or.aa_group(ii,1)%gtype=="CCYT") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="MET".or.aa_group(ii,1)%gtype=="NMET".or.aa_group(ii,1)%gtype=="CMET") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="SD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="SD") monitor(3)=grade_num(3)
				else
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="ASN".or.aa_group(ii,1)%gtype=="NASN".or.aa_group(ii,1)%gtype=="CASN") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="GLN".or.aa_group(ii,1)%gtype=="NGLN".or.aa_group(ii,1)%gtype=="CGLN") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				else	
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="ASP".or.aa_group(ii,1)%gtype=="NASP".or.aa_group(ii,1)%gtype=="CASP".or. &
			   aa_group(ii,1)%gtype=="ASH".or.aa_group(ii,1)%gtype=="NASH".or.aa_group(ii,1)%gtype=="CASH") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="GLU".or.aa_group(ii,1)%gtype=="NGLU".or.aa_group(ii,1)%gtype=="CGLU".or. &
			   aa_group(ii,1)%gtype=="GLH".or.aa_group(ii,1)%gtype=="NGLH".or.aa_group(ii,1)%gtype=="CGLH") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				else
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="HIE".or.aa_group(ii,1)%gtype=="NHIE".or.aa_group(ii,1)%gtype=="CHIE".or. &
			   aa_group(ii,1)%gtype=="HIP".or.aa_group(ii,1)%gtype=="NHIP".or.aa_group(ii,1)%gtype=="CHIP") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="LYS".or.aa_group(ii,1)%gtype=="NLYS".or.aa_group(ii,1)%gtype=="CLYS".or. &
			   aa_group(ii,1)%gtype=="LYN".or.aa_group(ii,1)%gtype=="NLYN".or.aa_group(ii,1)%gtype=="CLYN") then
			grade=4
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				elseif(aa_group(ii,1)%atype2(k)=="HD2".or.aa_group(ii,1)%atype2(k)=="HD3".or.aa_group(ii,1)%atype2(k)=="CE") then
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
					if(aa_group(ii,1)%atype2(k)=="CE") monitor(4)=grade_num(4)
				else
					grade_num(5)=grade_num(5)+1
					Iclass(5)%member(grade_num(5),1)=aa_group(ii,1)%coo2(k,1); Iclass(5)%member(grade_num(5),2)=aa_group(ii,1)%coo2(k,2); Iclass(5)%member(grade_num(5),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=5; index(k)%member_No=grade_num(5)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="ARG".or.aa_group(ii,1)%gtype=="NARG".or.aa_group(ii,1)%gtype=="CARG".or. &
			   aa_group(ii,1)%gtype=="ARN".or.aa_group(ii,1)%gtype=="NARN".or.aa_group(ii,1)%gtype=="CARN") then
			grade=4
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				elseif(aa_group(ii,1)%atype2(k)=="HD2".or.aa_group(ii,1)%atype2(k)=="HD3".or.aa_group(ii,1)%atype2(k)=="NE") then
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
					if(aa_group(ii,1)%atype2(k)=="NE") monitor(4)=grade_num(4)
				else
					grade_num(5)=grade_num(5)+1
					Iclass(5)%member(grade_num(5),1)=aa_group(ii,1)%coo2(k,1); Iclass(5)%member(grade_num(5),2)=aa_group(ii,1)%coo2(k,2); Iclass(5)%member(grade_num(5),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=5; index(k)%member_No=grade_num(5)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="PRO".or.aa_group(ii,1)%gtype=="NPRO".or.aa_group(ii,1)%gtype=="CPRO") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				else
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo
		endif
			
		if(grade.ne.aa_lib(ip)%grade) then
			open(10, file="error.txt", access="append")
				write(10,*) aa_lib(ip)%gtype
				write(10,*) "grade=", grade
				write(10,*) "aa_lib(",ip,")%grade=", aa_lib(ip)%grade
				write(10,*) "They are not equal with each other!"
			close(10)
			stop
		endif

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="CA") then
				CA(1)=aa_group(ii,1)%coo1(k,1); CA(2)=aa_group(ii,1)%coo1(k,2); CA(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo

		do n=2, rotanum
			aa_group(ii,n)%cnum1=aa_group(ii,1)%cnum1; aa_group(ii,n)%cnum2=aa_group(ii,1)%cnum2; aa_group(ii,n)%cnum3=aa_group(ii,1)%cnum3
			aa_group(ii,n)%gtype=name_original
			aa_group(ii,n)%atype1=aa_group(ii,1)%atype1; aa_group(ii,n)%atype2=aa_group(ii,1)%atype2; aa_group(ii,n)%atype3=aa_group(ii,1)%atype3
			aa_group(ii,n)%coo1=aa_group(ii,1)%coo1; aa_group(ii,n)%coo2=aa_group(ii,1)%coo2; aa_group(ii,n)%coo3=aa_group(ii,1)%coo3

			Tclass=Iclass	
			do j=1, grade
				delta_chi=real(aa_lib(ip)%dihedralangle(n,j)-aa_lib(ip)%dihedralangle(1,j))
				cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)			
				if(j==1) then
					rotaxis_x=Tclass(j)%member(monitor(j),1)-CA(1)
					rotaxis_y=Tclass(j)%member(monitor(j),2)-CA(2)
					rotaxis_z=Tclass(j)%member(monitor(j),3)-CA(3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				else
					rotaxis_x=Tclass(j)%member(monitor(j),1)-Tclass(j-1)%member(monitor(j-1),1)
					rotaxis_y=Tclass(j)%member(monitor(j),2)-Tclass(j-1)%member(monitor(j-1),2)
					rotaxis_z=Tclass(j)%member(monitor(j),3)-Tclass(j-1)%member(monitor(j-1),3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				endif

				call axisrotation(rotaxis, cos_angle, sin_angle, m)
				
				do l=(j+1), (grade+1)
					do k=1, grade_num(l)
						Tclass(l)%member(k,1)=Tclass(l)%member(k,1)-Tclass(j)%member(monitor(j),1)
						Tclass(l)%member(k,2)=Tclass(l)%member(k,2)-Tclass(j)%member(monitor(j),2)
						Tclass(l)%member(k,3)=Tclass(l)%member(k,3)-Tclass(j)%member(monitor(j),3)
					enddo
					
					Tmember=matmul(Tclass(l)%member, m)
					Tclass(l)%member=Tmember
					
					do k=1, grade_num(l)
						Tclass(l)%member(k,1)=anint((Tclass(l)%member(k,1)+Tclass(j)%member(monitor(j),1))*1000)/1000
						Tclass(l)%member(k,2)=anint((Tclass(l)%member(k,2)+Tclass(j)%member(monitor(j),2))*1000)/1000				
						Tclass(l)%member(k,3)=anint((Tclass(l)%member(k,3)+Tclass(j)%member(monitor(j),3))*1000)/1000
					enddo
				enddo
			enddo
				
			do k=1, aa_group(ii,n)%cnum2
				aa_group(ii,n)%coo2(k,1)=Tclass(index(k)%class_No)%member(index(k)%member_No,1)
				aa_group(ii,n)%coo2(k,2)=Tclass(index(k)%class_No)%member(index(k)%member_No,2)
				aa_group(ii,n)%coo2(k,3)=Tclass(index(k)%class_No)%member(index(k)%member_No,3)
			enddo
		enddo
10	continue
	enddo

	return
	end subroutine findrotamer

	
	subroutine energy_parameter(group, group_para)
	implicit none
	integer							:: i, j, k, status, atomid
	real							:: charge, epsion, r, rborn, fs, dielecons
	character*4						:: lbres, igraph
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)

	type(energyparameters), dimension(:,:), allocatable &
									:: Tgroup_para

	allocate(Tgroup_para(repeated_unit,gnum))
	
	do i=1, repeated_unit
		do j=1, gnum
			open(10, file='lib/ForceField/'//trim(group(i,j)%gtype), status="old")
			read(10, *)
			do while(.true.)
				read(10, 20, iostat=status) lbres, igraph, charge, epsion, r, rborn, fs, dielecons, atomid
				if(status.ne.0) goto 30
				do k=1, group(i,j)%cnum1
					if(group(i,j)%atype1(k)==igraph) then
						Tgroup_para(i,j)%charge1(k)=charge
						Tgroup_para(i,j)%epsion1(k)=epsion
						Tgroup_para(i,j)%r1(k)=r
						Tgroup_para(i,j)%rborn1(k)=rborn
						Tgroup_para(i,j)%fs1(k)=fs
						Tgroup_para(i,j)%dielecons1(k)=dielecons
						Tgroup_para(i,j)%atomid1(k)=atomid
						goto 40
					endif
				enddo
				do k=1, group(i,j)%cnum2
					if(group(i,j)%atype2(k)==igraph) then
						Tgroup_para(i,j)%charge2(k)=charge
						Tgroup_para(i,j)%epsion2(k)=epsion
						Tgroup_para(i,j)%r2(k)=r
						Tgroup_para(i,j)%rborn2(k)=rborn
						Tgroup_para(i,j)%fs2(k)=fs
						Tgroup_para(i,j)%dielecons2(k)=dielecons
						Tgroup_para(i,j)%atomid2(k)=atomid
						goto 40
					endif
				enddo
				do k=1, group(i,j)%cnum3
					if(group(i,j)%atype3(k)==igraph) then
						Tgroup_para(i,j)%charge3(k)=charge
						Tgroup_para(i,j)%epsion3(k)=epsion
						Tgroup_para(i,j)%r3(k)=r
						Tgroup_para(i,j)%rborn3(k)=rborn
						Tgroup_para(i,j)%fs3(k)=fs
						Tgroup_para(i,j)%dielecons3(k)=dielecons
						Tgroup_para(i,j)%atomid3(k)=atomid
						goto 40
					endif
				enddo
40				continue
			enddo
30			continue
			close(10)

20	format(2a4, 6e16.8, i8)
		enddo
	enddo

	do i=1, repeated_unit
		do j=1, gnum
			do k=1, group(i,j)%cnum1
				if(Tgroup_para(i,j)%dielecons1(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, "and", Tgroup_para(i,j)%dielecons1(k), "has wrong force field parameter in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
			do k=1, group(i,j)%cnum2
				if(Tgroup_para(i,j)%dielecons2(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, "and", Tgroup_para(i,j)%dielecons2(k), "has wrong force field parameter in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
			do k=1, group(i,j)%cnum3
				if(Tgroup_para(i,j)%dielecons3(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, "and", Tgroup_para(i,j)%dielecons3(k), "has wrong force field parameter in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
		enddo
	enddo	
	group_para=Tgroup_para
	deallocate(Tgroup_para)
	
	return
	end subroutine energy_parameter

	
	subroutine atom_links(group, numex, inb, numex4, inb4)
	implicit none
	type atomlink      ! The Data type "atomlink" is used to store the neighboring atoms for each atom.
		integer				:: linknum
		integer				:: linkindex(4)
	end type

	integer							:: categoryID, i, ii, ic, j, k, status, i1, j1, i2, j2, atomid, natom
	integer							:: id, linknum, linkindex(4)
	integer							:: ipres, numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)

	type(groupdetails)				:: group(repeated_unit,gnum)
	type(atomlink)					:: atom(repeated_unit*atom_num)									

	natom=0
	do categoryID=1, num4category
		do i=1, selfassembly(categoryID)%num4peptides
			ii=selfassembly(categoryID)%peptideID(i)
			do ic=1, gnum
				ipres=natom
				open(10, file='lib/Atomlink/'//trim(group(ii,ic)%gtype), status="old")
				do while(.true.)
					read(10, 20, iostat=status) id, linknum, (linkindex(k), k=1, linknum)
					if(status.ne.0) goto 30
					atomid=ipres+id
					atom(atomid)%linknum=linknum
					do k=1, linknum
						atom(atomid)%linkindex(k)=ipres+linkindex(k)
					enddo
				enddo			
30				continue
				close(10)
				natom=atomid
			enddo
		enddo
	enddo
20  format(i6, i7, 4i3)	

	do i1=1, natom
		numex(i1)=0
		do j1=1, atom(i1)%linknum
			numex(i1)=numex(i1)+1
			inb(i1,numex(i1))=atom(i1)%linkindex(j1)
		enddo

		do j1=1, atom(i1)%linknum
			i2=atom(i1)%linkindex(j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 40
				do k=1, atom(i1)%linknum
					if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
				enddo
				numex(i1)=numex(i1)+1
				inb(i1,numex(i1))=atom(i2)%linkindex(j2)
40				continue
			enddo
		enddo
	enddo

	do i1=1, natom
		numex4(i1)=0
		do j1=1, numex(i1)
			i2=inb(i1,j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 50
				do k=1, numex(i1)
					if(atom(i2)%linkindex(j2).eq.inb(i1,k)) goto 50
				enddo
				numex4(i1)=numex4(i1)+1
				inb4(i1,numex4(i1))=atom(i2)%linkindex(j2)
50			continue
			enddo
		enddo
	enddo

	return
	end subroutine atom_links

	
	subroutine atom_links4sidechain(chainID, ic, group, S_numex, S_inb, S_numex4, S_inb4)
	implicit none
	type atomlink
		integer				:: linknum
		integer				:: linkindex(4)
	end type
	
	integer							:: chainID, i, ii, j, k, status, i1, j1, i2, j2, atomid, natom
	integer							:: ic, id, linknum, linkindex(4)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(atomlink)					:: atom(50)

	ii=chainID
	open(10, file='lib/Atomlink/'//trim(group(ii,ic)%gtype), status="old")
	do while(.true.)
		read(10, 20, iostat=status) id, linknum, (linkindex(j), j=1, linknum)
		if(status.ne.0) goto 30
		atomid=id
		atom(atomid)%linknum=linknum
		do j=1, linknum
			atom(atomid)%linkindex(j)=linkindex(j)
		enddo
	enddo
30	continue
	close(10)
	natom=atomid
20  format(i6, i7, 4i3)

	do i1=1, natom
		S_numex(i1)=0
		do j1=1, atom(i1)%linknum
			S_numex(i1)=S_numex(i1)+1
			S_inb(i1,S_numex(i1))=atom(i1)%linkindex(j1)
		enddo

		do j1=1, atom(i1)%linknum
			i2=atom(i1)%linkindex(j1)
			if(i2.lt.1.or.i2.gt.natom) atom(i2)%linknum=0
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 40
				do k=1, atom(i1)%linknum
					if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
				enddo
				S_numex(i1)=S_numex(i1)+1
				S_inb(i1,S_numex(i1))=atom(i2)%linkindex(j2)
40				continue
			enddo
		enddo
	enddo

	do i1=1, natom
		S_numex4(i1)=0
		do j1=1, S_numex(i1)
			i2=S_inb(i1,j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 50
				do k=1, S_numex(i1)
					if(atom(i2)%linkindex(j2).eq.S_inb(i1,k)) goto 50
				enddo
				S_numex4(i1)=S_numex4(i1)+1
				S_inb4(i1,S_numex4(i1))=atom(i2)%linkindex(j2)
50			continue
			enddo
		enddo
	enddo

	return
	end subroutine atom_links4sidechain


	subroutine mc_choose_aminoacid(ic, group, aminoacid_name)
	implicit none
	integer						:: ii, ic, ip, ip1
	real						:: ran2
	character*4					:: aminoacid_name
	type(groupdetails)			:: group(repeated_unit,gnum)

	ii=1
	if(group(ii,ic)%gtype=="ALA".or.group(ii,ic)%gtype=="LEU".or.group(ii,ic)%gtype=="VAL".or.group(ii,ic)%gtype=="ILE".or.group(ii,ic)%gtype=="MET".or. &
	   group(ii,ic)%gtype=="PHE".or.group(ii,ic)%gtype=="TYR".or.group(ii,ic)%gtype=="TRP".or.group(ii,ic)%gtype=="GLY") then
		ip=10
	elseif(group(ii,ic)%gtype=="ASN".or.group(ii,ic)%gtype=="GLN".or.group(ii,ic)%gtype=="SER".or.group(ii,ic)%gtype=="THR".or.group(ii,ic)%gtype=="HIE") then
		ip=20
	elseif(group(ii,ic)%gtype=="ARG".or.group(ii,ic)%gtype=="LYS".or.group(ii,ic)%gtype=="GLU".or.group(ii,ic)%gtype=="ASP") then
		ip=30
	elseif(group(ii,ic)%gtype=="PRO".or.group(ii,ic)%gtype=="CYS") then
		ip=40
	elseif(group(ii,ic)%gtype=="NALA".or.group(ii,ic)%gtype=="NLEU".or.group(ii,ic)%gtype=="NVAL".or.group(ii,ic)%gtype=="NILE".or.group(ii,ic)%gtype=="NMET".or. &
		   group(ii,ic)%gtype=="NPHE".or.group(ii,ic)%gtype=="NTYR".or.group(ii,ic)%gtype=="NTRP".or.group(ii,ic)%gtype=="NGLY") then
		ip=11
	elseif(group(ii,ic)%gtype=="NASN".or.group(ii,ic)%gtype=="NGLN".or.group(ii,ic)%gtype=="NSER".or.group(ii,ic)%gtype=="NTHR".or.group(ii,ic)%gtype=="NHIE") then
		ip=21
	elseif(group(ii,ic)%gtype=="NARG".or.group(ii,ic)%gtype=="NLYS".or.group(ii,ic)%gtype=="NGLU".or.group(ii,ic)%gtype=="NASP") then
		ip=31
	elseif(group(ii,ic)%gtype=="NPRO".or.group(ii,ic)%gtype=="NCYS") then
		ip=41
	elseif(group(ii,ic)%gtype=="CALA".or.group(ii,ic)%gtype=="CLEU".or.group(ii,ic)%gtype=="CVAL".or.group(ii,ic)%gtype=="CILE".or.group(ii,ic)%gtype=="CMET".or. &
		   group(ii,ic)%gtype=="CPHE".or.group(ii,ic)%gtype=="CTYR".or.group(ii,ic)%gtype=="CTRP".or.group(ii,ic)%gtype=="CGLY") then
		ip=12
	elseif(group(ii,ic)%gtype=="CASN".or.group(ii,ic)%gtype=="CGLN".or.group(ii,ic)%gtype=="CSER".or.group(ii,ic)%gtype=="CTHR".or.group(ii,ic)%gtype=="CHIE") then
		ip=22
	elseif(group(ii,ic)%gtype=="CARG".or.group(ii,ic)%gtype=="CLYS".or.group(ii,ic)%gtype=="CGLU".or.group(ii,ic)%gtype=="CASP") then
		ip=32
	elseif(group(ii,ic)%gtype=="CPRO".or.group(ii,ic)%gtype=="CCYS") then
		ip=42
	endif

40	continue
	if(ip.eq.10) then
		call ran_gen(ran2,0)
		ip1=int(ran2*9-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="ALA"
		elseif(ip1.eq.2) then
			aminoacid_name="LEU"
		elseif(ip1.eq.3) then
			aminoacid_name="VAL"
		elseif(ip1.eq.4) then
			aminoacid_name="ILE"
		elseif(ip1.eq.5) then
			aminoacid_name="MET"
		elseif(ip1.eq.6) then
			aminoacid_name="PHE"
		elseif(ip1.eq.7) then
			aminoacid_name="TYR"
		elseif(ip1.eq.8) then
			aminoacid_name="TRP"
		elseif(ip1.eq.9) then
			aminoacid_name="GLY"
		endif
	elseif(ip.eq.20) then
		call ran_gen(ran2,0)
		ip1=int(ran2*5-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="SER"
		elseif(ip1.eq.2) then
			aminoacid_name="THR"
		elseif(ip1.eq.3) then
			aminoacid_name="HIE"
		elseif(ip1.eq.4) then
			aminoacid_name="ASN"
		elseif(ip1.eq.5) then
			aminoacid_name="GLN"
		endif
	elseif(ip.eq.30) then
		call ran_gen(ran2,0)
		ip1=int(ran2*4-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="ARG"
		elseif(ip1.eq.2) then
			aminoacid_name="LYS"
		elseif(ip1.eq.3) then
			aminoacid_name="GLU"
		elseif(ip1.eq.4) then
			aminoacid_name="ASP"
		endif
	elseif(ip.eq.40) then
		call ran_gen(ran2,0)
		ip1=int(ran2*2-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="PRO"
		elseif(ip1.eq.2) then
			aminoacid_name="CYS"
		endif
	elseif(ip.eq.11) then
		call ran_gen(ran2,0)
		ip1=int(ran2*9-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="NALA"
		elseif(ip1.eq.2) then
			aminoacid_name="NLEU"
		elseif(ip1.eq.3) then
			aminoacid_name="NVAL"
		elseif(ip1.eq.4) then
			aminoacid_name="NILE"
		elseif(ip1.eq.5) then
			aminoacid_name="NMET"
		elseif(ip1.eq.6) then
			aminoacid_name="NPHE"
		elseif(ip1.eq.7) then
			aminoacid_name="NTYR"
		elseif(ip1.eq.8) then
			aminoacid_name="NTRP"
		elseif(ip1.eq.9) then
			aminoacid_name="NGLY"
		endif
	elseif(ip.eq.21) then
		call ran_gen(ran2,0)
		ip1=int(ran2*5-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="NSER"
		elseif(ip1.eq.2) then
			aminoacid_name="NTHR"
		elseif(ip1.eq.3) then
			aminoacid_name="NHIE"
		elseif(ip1.eq.4) then
			aminoacid_name="NASN"
		elseif(ip1.eq.5) then
			aminoacid_name="NGLN"
		endif
	elseif(ip.eq.31) then
		call ran_gen(ran2,0)
		ip1=int(ran2*4-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="NARG"
		elseif(ip1.eq.2) then
			aminoacid_name="NLYS"
		elseif(ip1.eq.3) then
			aminoacid_name="NGLU"
		elseif(ip1.eq.4) then
			aminoacid_name="NASP"
		endif
	elseif(ip.eq.41) then
		call ran_gen(ran2,0)
		ip1=int(ran2*2-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="NPRO"
		elseif(ip1.eq.2) then
			aminoacid_name="NCYS"
		endif
	elseif(ip.eq.12) then
		call ran_gen(ran2,0)
		ip1=int(ran2*9-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="CALA"
		elseif(ip1.eq.2) then
			aminoacid_name="CLEU"
		elseif(ip1.eq.3) then
			aminoacid_name="CVAL"
		elseif(ip1.eq.4) then
			aminoacid_name="CILE"
		elseif(ip1.eq.5) then
			aminoacid_name="CMET"
		elseif(ip1.eq.6) then
			aminoacid_name="CPHE"
		elseif(ip1.eq.7) then
			aminoacid_name="CTYR"
		elseif(ip1.eq.8) then
			aminoacid_name="CTRP"
		elseif(ip1.eq.9) then
			aminoacid_name="CGLY"
		endif
	elseif(ip.eq.22) then
		call ran_gen(ran2,0)
		ip1=int(ran2*5-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="CSER"
		elseif(ip1.eq.2) then
			aminoacid_name="CTHR"
		elseif(ip1.eq.3) then
			aminoacid_name="CHIE"
		elseif(ip1.eq.4) then
			aminoacid_name="CASN"
		elseif(ip1.eq.5) then
			aminoacid_name="CGLN"			
		endif
	elseif(ip.eq.32) then
		call ran_gen(ran2,0)
		ip1=int(ran2*4-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="CARG"
		elseif(ip1.eq.2) then
			aminoacid_name="CLYS"
		elseif(ip1.eq.3) then
			aminoacid_name="CGLU"
		elseif(ip1.eq.4) then
			aminoacid_name="CASP"		
		endif
	elseif(ip.eq.42) then
		call ran_gen(ran2,0)
		ip1=int(ran2*2-1.0e-3)+1
		if(ip1.eq.1) then
			aminoacid_name="CPRO"
		elseif(ip1.eq.2) then
			aminoacid_name="CCYS"
		endif
	endif

	return
	end subroutine mc_choose_aminoacid	

	
	subroutine groupinfo(name, group_name, flag)
	implicit none
	integer					:: i, flag
	character*4				:: name, group_name(3)

	if(name=="GLY".or.name=="NGLY".or.name=="CGLY") then
		group_name(1)="GLY"
		group_name(2)="NGLY"
		group_name(3)="CGLY"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LEU".or.name=="NLEU".or.name=="CLEU") then
		group_name(1)="LEU"
		group_name(2)="NLEU"
		group_name(3)="CLEU"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="VAL".or.name=="NVAL".or.name=="CVAL") then
		group_name(1)="VAL"
		group_name(2)="NVAL"
		group_name(3)="CVAL"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ILE".or.name=="NILE".or.name=="CILE") then
		group_name(1)="ILE"
		group_name(2)="NILE"
		group_name(3)="CILE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="MET".or.name=="NMET".or.name=="CMET") then
		group_name(1)="MET"
		group_name(2)="NMET"
		group_name(3)="CMET"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="PHE".or.name=="NPHE".or.name=="CPHE") then
		group_name(1)="PHE"
		group_name(2)="NPHE"
		group_name(3)="CPHE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TYR".or.name=="NTYR".or.name=="CTYR") then
		group_name(1)="TYR"
		group_name(2)="NTYR"
		group_name(3)="CTYR"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TYX".or.name=="NTYX".or.name=="CTYX") then
		group_name(1)="TYX"
		group_name(2)="NTYX"
		group_name(3)="CTYX"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TRP".or.name=="NTRP".or.name=="CTRP") then
		group_name(1)="TRP"
		group_name(2)="NTRP"
		group_name(3)="CTRP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ARG".or.name=="NARG".or.name=="CARG") then
		group_name(1)="ARG"
		group_name(2)="NARG"
		group_name(3)="CARG"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ARN".or.name=="NARN".or.name=="CARN") then
		group_name(1)="ARN"
		group_name(2)="NARN"
		group_name(3)="CARN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LYN".or.name=="NLYN".or.name=="CLYN") then
		group_name(1)="LYN"
		group_name(2)="NLYN"
		group_name(3)="CLYN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LYS".or.name=="NLYS".or.name=="CLYS") then
		group_name(1)="LYS"
		group_name(2)="NLYS"
		group_name(3)="CLYS"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="SER".or.name=="NSER".or.name=="CSER") then
		group_name(1)="SER"
		group_name(2)="NSER"
		group_name(3)="CSER"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="THR".or.name=="NTHR".or.name=="CTHR") then
		group_name(1)="THR"
		group_name(2)="NTHR"
		group_name(3)="CTHR"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASN".or.name=="NASN".or.name=="CASN") then
		group_name(1)="ASN"
		group_name(2)="NASN"
		group_name(3)="CASN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLN".or.name=="NGLN".or.name=="CGLN") then
		group_name(1)="GLN"
		group_name(2)="NGLN"
		group_name(3)="CGLN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="HIE".or.name=="NHIE".or.name=="CHIE") then
		group_name(1)="HIE"
		group_name(2)="NHIE"
		group_name(3)="CHIE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="HIP".or.name=="NHIP".or.name=="CHIP") then
		group_name(1)="HIP"
		group_name(2)="NHIP"
		group_name(3)="CHIP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
		group_name(1)="PRO"
		group_name(2)="NPRO"
		group_name(3)="CPRO"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="CYS".or.name=="NCYS".or.name=="CCYS") then
		group_name(1)="CYS"
		group_name(2)="NCYS"
		group_name(3)="CCYS"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="CYT".or.name=="NCYT".or.name=="CCYT") then
		group_name(1)="CYT"
		group_name(2)="NCYT"
		group_name(3)="CCYT"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ALA".or.name=="NALA".or.name=="CALA") then
		group_name(1)="ALA"
		group_name(2)="NALA"
		group_name(3)="CALA"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLH".or.name=="NGLH".or.name=="CGLH") then
		group_name(1)="GLH"
		group_name(2)="NGLH"
		group_name(3)="CGLH"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLU".or.name=="NGLU".or.name=="CGLU") then
		group_name(1)="GLU"
		group_name(2)="NGLU"
		group_name(3)="CGLU"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASH".or.name=="NASH".or.name=="CASH") then
		group_name(1)="ASH"
		group_name(2)="NASH"
		group_name(3)="CASH"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASP".or.name=="NASP".or.name=="CASP") then
		group_name(1)="ASP"
		group_name(2)="NASP"
		group_name(3)="CASP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	endif
5	continue

	return
	end subroutine groupinfo

	
	subroutine scmf_choose_aminoacid(ip, aminoacid_number, aminoacid_name)
	implicit none
	integer					:: aminoacid_number
	integer					:: ip, ip1, ip2, i
	real					:: ran2
	character*4				:: char, aminoacid_name(10)

	if(ip==1) then
		aminoacid_number=9
		aminoacid_name(1)="ALA"
		aminoacid_name(2)="LEU"
		aminoacid_name(3)="VAL"
		aminoacid_name(4)="ILE"
		aminoacid_name(5)="MET"
		aminoacid_name(6)="PHE"
		aminoacid_name(7)="TYR"
		aminoacid_name(8)="TRP"
		aminoacid_name(9)="GLY"	
	elseif(ip==2) then
		aminoacid_number=5
		aminoacid_name(1)="ASN"
		aminoacid_name(2)="GLN"
		aminoacid_name(3)="SER"
		aminoacid_name(4)="THR"
		aminoacid_name(5)="HIE"
	elseif(ip==3) then
		aminoacid_number=4
		aminoacid_name(1)="ARG"
		aminoacid_name(2)="LYS"
		aminoacid_name(3)="GLU"
		aminoacid_name(4)="ASP"
	elseif(ip==4) then
		aminoacid_number=2
		aminoacid_name(1)="PRO"
		aminoacid_name(2)="CYS"
	endif

	do i=1, (aminoacid_number-1)
		call ran_gen(ran2,0)
		ip1=int(ran2*aminoacid_number-1.0e-3)+1
		if(ip1.gt.aminoacid_number) ip1=aminoacid_number
		call ran_gen(ran2,0)
		ip2=int(ran2*aminoacid_number-1.0e-3)+1
		if(ip2.gt.aminoacid_number) ip2=aminoacid_number

		char=aminoacid_name(ip1)
		aminoacid_name(ip2)=aminoacid_name(ip1)
		aminoacid_name(ip1)=char
	enddo

	return
	end subroutine scmf_choose_aminoacid
	

	subroutine sidechain_category(chainID, ic, group, Iclass, grade, grade_num, index, monitor)
	implicit none
	integer								:: grade, grade_num(6), monitor(6), chainID, i, ii, ic
	type(groupdetails)					:: group(repeated_unit,gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6)

	ii=chainID
	grade_num=0
	if(group(ii,ic)%gtype=="VAL".or.group(ii,ic)%gtype=="NVAL".or.group(ii,ic)%gtype=="CVAL") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1)
				Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2)
				Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo				
	elseif(group(ii,ic)%gtype=="LEU".or.group(ii,ic)%gtype=="NLEU".or.group(ii,ic)%gtype=="CLEU") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo		
	elseif(group(ii,ic)%gtype=="ILE".or.group(ii,ic)%gtype=="NILE".or.group(ii,ic)%gtype=="CILE") then
		grade=2
		do i=1, group(ii,ic)%cnum2	
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3) 
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB".or.group(ii,ic)%atype2(i)=="CG2".or.group(ii,ic)%atype2(i)=="HG21".or.group(ii,ic)%atype2(i)=="HG22".or. &
			       group(ii,ic)%atype2(i)=="HG23".or.group(ii,ic)%atype2(i)=="CG1") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG1") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo		
	elseif(group(ii,ic)%gtype=="PHE".or.group(ii,ic)%gtype=="NPHE".or.group(ii,ic)%gtype=="CPHE") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="TRP".or.group(ii,ic)%gtype=="NTRP".or.group(ii,ic)%gtype=="CTRP") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="TYR".or.group(ii,ic)%gtype=="NTYR".or.group(ii,ic)%gtype=="CTYR".or.  &
	       group(ii,ic)%gtype=="TYX".or.group(ii,ic)%gtype=="NTYX".or.group(ii,ic)%gtype=="CTYX") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="SER".or.group(ii,ic)%gtype=="NSER".or.group(ii,ic)%gtype=="CSER") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="THR".or.group(ii,ic)%gtype=="NTHR".or.group(ii,ic)%gtype=="CTHR") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="CYS".or.group(ii,ic)%gtype=="NCYS".or.group(ii,ic)%gtype=="CCYS".or.  &
	       group(ii,ic)%gtype=="CYT".or.group(ii,ic)%gtype=="NCYT".or.group(ii,ic)%gtype=="CCYT") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="MET".or.group(ii,ic)%gtype=="NMET".or.group(ii,ic)%gtype=="CMET") then
		grade=3
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="SD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="SD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="ASN".or.group(ii,ic)%gtype=="NASN".or.group(ii,ic)%gtype=="CASN") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="GLN".or.group(ii,ic)%gtype=="NGLN".or.group(ii,ic)%gtype=="CGLN") then
		grade=3
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else	
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="ASP".or.group(ii,ic)%gtype=="NASP".or.group(ii,ic)%gtype=="CASP".or.  &
	       group(ii,ic)%gtype=="ASH".or.group(ii,ic)%gtype=="NASH".or.group(ii,ic)%gtype=="CASH") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="GLU".or.group(ii,ic)%gtype=="NGLU".or.group(ii,ic)%gtype=="CGLU".or.  &
	       group(ii,ic)%gtype=="GLH".or.group(ii,ic)%gtype=="NGLH".or.group(ii,ic)%gtype=="CGLH") then
		grade=3
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="HIE".or.group(ii,ic)%gtype=="NHIE".or.group(ii,ic)%gtype=="CHIE".or.  &
	       group(ii,ic)%gtype=="HIP".or.group(ii,ic)%gtype=="NHIP".or.group(ii,ic)%gtype=="CHIP") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="LYS".or.group(ii,ic)%gtype=="NLYS".or.group(ii,ic)%gtype=="CLYS".or.  &
	       group(ii,ic)%gtype=="LYN".or.group(ii,ic)%gtype=="NLYN".or.group(ii,ic)%gtype=="CLYN") then
		grade=4
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ii,ic)%atype2(i)=="HD2".or.group(ii,ic)%atype2(i)=="HD3".or.group(ii,ic)%atype2(i)=="CE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ii,ic)%atype2(i)=="CE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ii,ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ii,ic)%coo2(i,2); Iclass(5)%member(grade_num(5),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="ARG".or.group(ii,ic)%gtype=="NARG".or.group(ii,ic)%gtype=="CARG".or.  &
	       group(ii,ic)%gtype=="ARN".or.group(ii,ic)%gtype=="NARN".or.group(ii,ic)%gtype=="CARN") then
		grade=4
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ii,ic)%atype2(i)=="HD2".or.group(ii,ic)%atype2(i)=="HD3".or.group(ii,ic)%atype2(i)=="NE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ii,ic)%atype2(i)=="NE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ii,ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ii,ic)%coo2(i,2); Iclass(5)%member(grade_num(5),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	endif

	return
	end subroutine sidechain_category	

	
	subroutine dihedralangle_reading(gtype, dihedral_num, dihedral)
	implicit none
	integer								:: dihedral_num, i, j
	character*4							:: gtype
	type(dihedralparameters)			:: dihedral	

	open(10, file='lib/DihedralAngle/'//trim(gtype), status="old")
		read(10, "(i8)") dihedral_num
		do i=1, dihedral_num
			read(10,"(5i8)") dihedral%iph(i), dihedral%jph(i), dihedral%kph(i), dihedral%lph(i), dihedral%multiply(i)
			do j=1, dihedral%multiply(i)
				read(10,"(3e16.8)") dihedral%pk(i,j), dihedral%pn(i,j), dihedral%phase(i,j)
			enddo
		enddo
	close(10)
	
	return
	end subroutine dihedralangle_reading	


!	subroutine groupRES_1(ic, group, g_residue)
!	implicit none
!	integer							:: ic, chainID, i

!	type(groupdetails)				:: group(repeated_unit,gnum)
!	type(RES4chain)					:: g_residue	

!	chainID=1
!	g_residue%num=0
!	if(group(chainID,ic)%gtype=="ALA".or.group(chainID,ic)%gtype=="LEU".or.group(chainID,ic)%gtype=="VAL".or.group(chainID,ic)%gtype=="ILE".or.group(chainID,ic)%gtype=="MET".or. &
!	   group(chainID,ic)%gtype=="PHE".or.group(chainID,ic)%gtype=="TYR".or.group(chainID,ic)%gtype=="TRP".or.group(chainID,ic)%gtype=="GLY".or.group(chainID,ic)%gtype=="NALA".or.group(chainID,ic)%gtype=="NLEU".or. &
!	   group(chainID,ic)%gtype=="NVAL".or.group(chainID,ic)%gtype=="NILE".or.group(chainID,ic)%gtype=="NMET".or.group(chainID,ic)%gtype=="NPHE".or.group(chainID,ic)%gtype=="NTYR".or. &
!	   group(chainID,ic)%gtype=="NTRP".or.group(chainID,ic)%gtype=="NGLY".or.group(chainID,ic)%gtype=="CALA".or.group(chainID,ic)%gtype=="CLEU".or.group(chainID,ic)%gtype=="CVAL".or.group(chainID,ic)%gtype=="CILE".or. &
!	   group(chainID,ic)%gtype=="CMET".or.group(chainID,ic)%gtype=="CPHE".or.group(chainID,ic)%gtype=="CTYR".or.group(chainID,ic)%gtype=="CTRP".or.group(chainID,ic)%gtype=="CGLY") then
!	   	do i=1,gnum
!			if(group(chainID,i)%gtype=="ALA".or.group(chainID,i)%gtype=="LEU".or.group(chainID,i)%gtype=="VAL".or.group(chainID,i)%gtype=="ILE".or.group(chainID,i)%gtype=="MET".or. &
!			   group(chainID,i)%gtype=="PHE".or.group(chainID,i)%gtype=="TYR".or.group(chainID,i)%gtype=="TRP".or.group(chainID,i)%gtype=="GLY".or.group(chainID,i)%gtype=="NALA".or.group(chainID,i)%gtype=="NLEU".or. &
!			   group(chainID,i)%gtype=="NVAL".or.group(chainID,i)%gtype=="NILE".or.group(chainID,i)%gtype=="NMET".or.group(chainID,i)%gtype=="NPHE".or.group(chainID,i)%gtype=="NTYR".or. &
!			   group(chainID,i)%gtype=="NTRP".or.group(chainID,i)%gtype=="NGLY".or.group(chainID,i)%gtype=="CALA".or.group(chainID,i)%gtype=="CLEU".or.group(chainID,i)%gtype=="CVAL".or.group(chainID,i)%gtype=="CILE".or. &
!			   group(chainID,i)%gtype=="CMET".or.group(chainID,i)%gtype=="CPHE".or.group(chainID,i)%gtype=="CTYR".or.group(chainID,i)%gtype=="CTRP".or.group(chainID,i)%gtype=="CGLY") then
!				g_residue%num=g_residue%num+1
!				g_residue%IDs(g_residue%num)=i
!			endif
!		enddo
!	elseif(group(chainID,ic)%gtype=="GLU".or.group(chainID,ic)%gtype=="ASP".or.group(chainID,ic)%gtype=="ARG".or.group(chainID,ic)%gtype=="LYS".or.group(chainID,ic)%gtype=="ASN".or. &
!		   group(chainID,ic)%gtype=="GLN".or.group(chainID,ic)%gtype=="SER".or.group(chainID,ic)%gtype=="THR".or.group(chainID,ic)%gtype=="HIE".or.group(chainID,ic)%gtype=="NGLU".or. &
!		   group(chainID,ic)%gtype=="NASP".or.group(chainID,ic)%gtype=="NARG".or.group(chainID,ic)%gtype=="NLYS".or.group(chainID,ic)%gtype=="NASN".or.group(chainID,ic)%gtype=="NGLN".or. &
!		   group(chainID,ic)%gtype=="NSER".or.group(chainID,ic)%gtype=="NTHR".or.group(chainID,ic)%gtype=="NHIE".or.group(chainID,ic)%gtype=="CGLU".or.group(chainID,ic)%gtype=="CASP".or. &
!		   group(chainID,ic)%gtype=="CARG".or.group(chainID,ic)%gtype=="CLYS".or.group(chainID,ic)%gtype=="CASN".or.group(chainID,ic)%gtype=="CGLN".or.group(chainID,ic)%gtype=="CSER".or. &
!		   group(chainID,ic)%gtype=="CTHR".or.group(chainID,ic)%gtype=="CHIE") then
!		do i=1,gnum
!			if(group(chainID,i)%gtype=="GLU".or.group(chainID,i)%gtype=="ASP".or.group(chainID,i)%gtype=="ARG".or.group(chainID,i)%gtype=="LYS".or.group(chainID,i)%gtype=="ASN".or. &
!			   group(chainID,i)%gtype=="GLN".or.group(chainID,i)%gtype=="SER".or.group(chainID,i)%gtype=="THR".or.group(chainID,i)%gtype=="HIE".or.group(chainID,i)%gtype=="NGLU".or. &
!			   group(chainID,i)%gtype=="NASP".or.group(chainID,i)%gtype=="NARG".or.group(chainID,i)%gtype=="NLYS".or.group(chainID,i)%gtype=="NASN".or.group(chainID,i)%gtype=="NGLN".or. &
!			   group(chainID,i)%gtype=="NSER".or.group(chainID,i)%gtype=="NTHR".or.group(chainID,i)%gtype=="NHIE".or.group(chainID,i)%gtype=="CGLU".or.group(chainID,i)%gtype=="CASP".or. &
!			   group(chainID,i)%gtype=="CARG".or.group(chainID,i)%gtype=="CLYS".or.group(chainID,i)%gtype=="CASN".or.group(chainID,i)%gtype=="CGLN".or.group(chainID,i)%gtype=="CSER".or. &
!			   group(chainID,i)%gtype=="CTHR".or.group(chainID,i)%gtype=="CHIE") then
!				g_residue%num=g_residue%num+1
!				g_residue%IDs(g_residue%num)=i
!			endif
!		enddo
!	endif
	
!	return
!	end subroutine groupRES_1

	
end module database

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module transplant
	
	use constant
	use mathfunction

	contains
	subroutine residue_replace(chainID, ic, group, ip, aa_group, temp_group)
	implicit none
	integer								:: chainID, ic, ip, ii, j, k, flag
	type(groupdetails)					:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), aa_group(repeated_unit,40)

	temp_group=group

	ii=chainID
	if(temp_group(ii,ic)%gtype=="PRO".or.temp_group(ii,ic)%gtype=="NPRO".or.temp_group(ii,ic)%gtype=="CPRO") then
		if(aa_group(ii,ip)%gtype=="PRO".or.aa_group(ii,ip)%gtype=="NPRO".or.aa_group(ii,ip)%gtype=="CPRO") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		elseif(aa_group(ii,ip)%gtype=="GLY".or.aa_group(ii,ip)%gtype=="NGLY".or.aa_group(ii,ip)%gtype=="CGLY") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, temp_group(ii,ic)%cnum1
				if(temp_group(ii,ic)%atype1(j)=="H2".or.temp_group(ii,ic)%atype1(j)=="H3") then
					flag=1
				endif
			enddo

			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="H") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="N") then
							if(flag==1) then
								temp_group(ii,ic)%atype1(k+1)="H1"
							else
								temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							endif
							goto 10
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				elseif(aa_group(ii,ip)%atype1(j)=="HA2") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="CA") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 10
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo	
				elseif(aa_group(ii,ip)%atype1(j)=="HA3") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="HA2") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 10
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				endif					
10				continue
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
		else
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, temp_group(ii,ic)%cnum1
				if(temp_group(ii,ic)%atype1(j)=="H2".or.temp_group(ii,ic)%atype1(j)=="H3") then
					flag=1
				endif
			enddo

			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="H") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="N") then
							if(flag==1) then
								temp_group(ii,ic)%atype1(k+1)="H1"
							else
								temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							endif
							goto 20			
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				endif
			enddo
20			continue
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		endif	
	elseif(temp_group(ii,ic)%gtype=="GLY".or.temp_group(ii,ic)%gtype=="NGLY".or.temp_group(ii,ic)%gtype=="CGLY") then
		if(aa_group(ii,ip)%gtype=="PRO".or.aa_group(ii,ip)%gtype=="NPRO".or.aa_group(ii,ip)%gtype=="CPRO") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="H1".or.temp_group(ii,ic)%atype1(j)=="H".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1

			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="HA2") then
					temp_group(ii,ic)%atype1(j)="HA"
				elseif(temp_group(ii,ic)%atype1(j)=="HA3".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		elseif(aa_group(ii,ip)%gtype=="GLY".or.aa_group(ii,ip)%gtype=="NGLY".or.aa_group(ii,ip)%gtype=="CGLY") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype			
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		else
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype			
			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="HA2") then
					temp_group(ii,ic)%atype1(j)="HA"
				elseif(temp_group(ii,ic)%atype1(j)=="HA3".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1											
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		endif			
	else
		if(aa_group(ii,ip)%gtype=="PRO".or.aa_group(ii,ip)%gtype=="NPRO".or.aa_group(ii,ip)%gtype=="CPRO") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="H1".or.temp_group(ii,ic)%atype1(j)=="H".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo			
		elseif(aa_group(ii,ip)%gtype=="GLY".or.aa_group(ii,ip)%gtype=="NGLY".or.aa_group(ii,ip)%gtype=="CGLY") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="HA2") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="CA") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 30
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo										
				elseif(aa_group(ii,ip)%atype1(j)=="HA3") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="HA2") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 30
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				endif					
30				continue
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
		else
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype		
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		endif
	endif

	return
	end subroutine residue_replace

	
	subroutine backup4sidechain(flag, chainID, ic, group, aa_backup)
	implicit none
	integer							:: flag, chainID, ic, ii, ik
	type(groupdetails)				:: group(repeated_unit,gnum), aa_backup

	ii=chainID
	if(flag==0) then
		aa_backup%gtype=group(ii,ic)%gtype

		aa_backup%cnum1=group(ii,ic)%cnum1
		do ik=1, group(ii,ic)%cnum1
			aa_backup%atype1(ik)=group(ii,ic)%atype1(ik)
			aa_backup%coo1(ik,1)=group(ii,ic)%coo1(ik,1)
			aa_backup%coo1(ik,2)=group(ii,ic)%coo1(ik,2)
			aa_backup%coo1(ik,3)=group(ii,ic)%coo1(ik,3)
		enddo

		aa_backup%cnum2=group(ii,ic)%cnum2
		do ik=1, group(ii,ic)%cnum2
			aa_backup%atype2(ik)=group(ii,ic)%atype2(ik)
			aa_backup%coo2(ik,1)=group(ii,ic)%coo2(ik,1)
			aa_backup%coo2(ik,2)=group(ii,ic)%coo2(ik,2)
			aa_backup%coo2(ik,3)=group(ii,ic)%coo2(ik,3)
		enddo

		aa_backup%cnum3=group(ii,ic)%cnum3
		do ik=1, group(ii,ic)%cnum3
			aa_backup%atype3(ik)=group(ii,ic)%atype3(ik)
			aa_backup%coo3(ik,1)=group(ii,ic)%coo3(ik,1)
			aa_backup%coo3(ik,2)=group(ii,ic)%coo3(ik,2)
			aa_backup%coo3(ik,3)=group(ii,ic)%coo3(ik,3)
		enddo

	elseif(flag==1) then
		group(ii,ic)%gtype=aa_backup%gtype

		group(ii,ic)%cnum1=aa_backup%cnum1
		do ik=1, aa_backup%cnum1
			group(ii,ic)%atype1(ik)=aa_backup%atype1(ik)
			group(ii,ic)%coo1(ik,1)=aa_backup%coo1(ik,1)
			group(ii,ic)%coo1(ik,2)=aa_backup%coo1(ik,2)
			group(ii,ic)%coo1(ik,3)=aa_backup%coo1(ik,3)
		enddo

		group(ii,ic)%cnum2=aa_backup%cnum2
		do ik=1, aa_backup%cnum2
			group(ii,ic)%atype2(ik)=aa_backup%atype2(ik)
			group(ii,ic)%coo2(ik,1)=aa_backup%coo2(ik,1)
			group(ii,ic)%coo2(ik,2)=aa_backup%coo2(ik,2)
			group(ii,ic)%coo2(ik,3)=aa_backup%coo2(ik,3)
		enddo

		group(ii,ic)%cnum3=aa_backup%cnum3
		do ik=1, aa_backup%cnum3
			group(ii,ic)%atype3(ik)=aa_backup%atype3(ik)
			group(ii,ic)%coo3(ik,1)=aa_backup%coo3(ik,1)
			group(ii,ic)%coo3(ik,2)=aa_backup%coo3(ik,2)
			group(ii,ic)%coo3(ik,3)=aa_backup%coo3(ik,3)
		enddo
	endif
	
	return
	end subroutine backup4sidechain

	
	subroutine check_transplant(chainID, ic, group, feedback)
	implicit none
	integer							:: chainID, ii, ic, ik, jj, jc, jk, feedback
	real							:: rij
	type(groupdetails)				:: group(repeated_unit,gnum)

	feedback=1
	ii=chainID
	do ik=1, group(ii,ic)%cnum2
		do jj=1, repeated_unit
			do jc=1, gnum
				if (jj.eq.ii.and.jc.eq.ic) goto 20
				do jk=1, group(jj,jc)%cnum1
					if(jj==ii.and.jc==(ic+1).and.group(ii,ic)%atype2(ik)=="CB".and.group(jj,jc)%atype1(jk)=="N") goto 50
					if(group(ii,ic)%gtype=="PRO".or.group(ii,ic)%gtype=="NPRO".or.group(ii,ic)%gtype=="CPRO") then
						if(jj==ii.and.jc==(ic-1).and.group(ii,ic)%atype2(ik)=="CD".and.group(jj,jc)%atype1(jk)=="CA") goto 50
					endif
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
						(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
					rij=sqrt(rij)
					if(rij<1.55) then
						feedback=0
						goto 10
					endif
50					continue
				enddo
				do jk=1, group(jj,jc)%cnum2
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
						(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
					rij=sqrt(rij)
					if(rij<1.55) then
						feedback=0
						goto 10
					endif
				enddo
				do jk=1, group(jj,jc)%cnum3
					if(jj==ii.and.jc==(ic-1).and.group(ii,ic)%atype2(ik)=="CB".and.group(jj,jc)%atype3(jk)=="C") goto 60
					if(group(ii,ic)%gtype=="PRO".or.group(ii,ic)%gtype=="NPRO".or.group(ii,ic)%gtype=="CPRO") then
						if(jj==ii.and.jc==(ic-1).and.group(ii,ic)%atype2(ik)=="CD") then
							goto 60
						elseif(jj==ii.and.jc==(ic-1).and.(group(ii,ic)%atype2(ik)=="CG".or.group(ii,ic)%atype2(ik)=="HD2".or.group(ii,ic)%atype2(ik)=="HD3")) then
							if(group(jj,jc)%atype3(jk)=="C") goto 60
						endif
					endif
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
						(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
					rij=sqrt(rij)
					if(rij<1.55) then
						feedback=0
						goto 10
					endif
60		    		continue
				enddo
20				continue
			enddo
		enddo
	enddo
10	continue

	return
	end subroutine check_transplant


	subroutine torsionangle4sidechain(group, chainID, ic, grade, dihedral, Dihedral4entropy)
	implicit none
	integer								:: chainID, ic, i, j
	integer								:: natom, grade, ip, jp, kp, lp
	real								:: p1(3), p2(3), p3(3), p4(3), angle
	real							    :: mdcrd(60,3)
	real								:: Dihedral4entropy(4)
	type(groupdetails)					:: group(repeated_unit,gnum)
	type(dihedralparameters)		    :: dihedral

	natom=0
	do i=1, group(chainID,ic)%cnum1
		natom=natom+1
		mdcrd(natom,1)=group(chainID,ic)%coo1(i,1)
		mdcrd(natom,2)=group(chainID,ic)%coo1(i,2)
		mdcrd(natom,3)=group(chainID,ic)%coo1(i,3)
	enddo
	do i=1, group(chainID,ic)%cnum2
		natom=natom+1
		mdcrd(natom,1)=group(chainID,ic)%coo2(i,1)
		mdcrd(natom,2)=group(chainID,ic)%coo2(i,2)
		mdcrd(natom,3)=group(chainID,ic)%coo2(i,3)
	enddo
	do i=1, group(chainID,ic)%cnum3
		natom=natom+1
		mdcrd(natom,1)=group(chainID,ic)%coo3(i,1)
		mdcrd(natom,2)=group(chainID,ic)%coo3(i,2)
		mdcrd(natom,3)=group(chainID,ic)%coo3(i,3)
	enddo

	do j=1, grade
		ip=dihedral%iph(j); jp=dihedral%jph(j); kp=dihedral%kph(j); lp=dihedral%lph(j)
		p1(1)=mdcrd(ip,1); p1(2)=mdcrd(ip,2); p1(3)=mdcrd(ip,3)
		p2(1)=mdcrd(jp,1); p2(2)=mdcrd(jp,2); p2(3)=mdcrd(jp,3)
		p3(1)=mdcrd(kp,1); p3(2)=mdcrd(kp,2); p3(3)=mdcrd(kp,3)
		p4(1)=mdcrd(lp,1); p4(2)=mdcrd(lp,2); p4(3)=mdcrd(lp,3)	
		call phipsiomg_angle(p1,p2,p3,p4,angle)
		Dihedral4entropy(j)=angle
	enddo

	return
	end subroutine torsionangle4sidechain

end module transplant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module energy_calculation

	use constant
	use mathfunction
	use database

	contains
	subroutine vdwenergy(chainID, ic, group, group_para, energy)
	implicit none
	integer							:: chainID, ii, ic, ik, jj, jc, jk
	real							:: energy, rij, epsion, r0, acoeff, bcoeff, vdw
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)

	energy=0.0
	ii=chainID
	do ik=1, group(ii,ic)%cnum2
		do jj=1, repeated_unit
			do jc=1, gnum
				if (jj==ii.and.(jc==(ic-1).or.jc==ic.or.jc==(ic+1))) goto 10

				do jk=1, group(jj,jc)%cnum1
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
					    (group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
					if(rij>100.0) goto 20
					epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
					r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
20					continue
				enddo
				do jk=1, group(jj,jc)%cnum2
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
					    (group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
					if(rij>100.0) goto 30
					epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
					r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
30					continue
				enddo
				do jk=1, group(jj,jc)%cnum3
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
					    (group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
					if(rij>100.0) goto 40
					epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
					r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
40					continue
				enddo
10				continue
			enddo
		enddo
	enddo

	return
	end subroutine vdwenergy


	subroutine sidechain_energy(stage, chainID, ic, group, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy)	
	implicit none
	integer								:: chainID, j, k, l, ii, ic, ik, i_id, jj, jc, jk, j_id, flag, stage
	integer								:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	integer								:: natom, dihedral_num, ip, jp, kp, lp
	real								:: energy, rij, epsion, r0, acoeff, bcoeff, dielecons4solute, vdw, ele
	real								:: p1(3), p2(3), p3(3), p4(3), angle, dihedral_energy, potential
	real							    :: mdcrd(60,3)
	type(groupdetails)					:: group(repeated_unit,gnum)
	type(energyparameters)				:: group_para(repeated_unit,gnum)
	type(dihedralparameters)		    :: dihedral

	energy=0.0
	ii=chainID
	natom=0
	do ik=1, group(ii,ic)%cnum1
		natom=natom+1
		mdcrd(natom,1)=group(ii,ic)%coo1(ik,1)
		mdcrd(natom,2)=group(ii,ic)%coo1(ik,2)
		mdcrd(natom,3)=group(ii,ic)%coo1(ik,3)
	enddo
	do ik=1, group(ii,ic)%cnum2
		natom=natom+1
		mdcrd(natom,1)=group(ii,ic)%coo2(ik,1)
		mdcrd(natom,2)=group(ii,ic)%coo2(ik,2)
		mdcrd(natom,3)=group(ii,ic)%coo2(ik,3)
	enddo
	do ik=1, group(ii,ic)%cnum3
		natom=natom+1
		mdcrd(natom,1)=group(ii,ic)%coo3(ik,1)
		mdcrd(natom,2)=group(ii,ic)%coo3(ik,2)
		mdcrd(natom,3)=group(ii,ic)%coo3(ik,3)
	enddo

	dihedral_energy=0.0
	do j=1, dihedral_num
		ip=dihedral%iph(j); jp=dihedral%jph(j); kp=dihedral%kph(j); lp=dihedral%lph(j)
		p1(1)=mdcrd(ip,1); p1(2)=mdcrd(ip,2); p1(3)=mdcrd(ip,3)
		p2(1)=mdcrd(jp,1); p2(2)=mdcrd(jp,2); p2(3)=mdcrd(jp,3)
		p3(1)=mdcrd(kp,1); p3(2)=mdcrd(kp,2); p3(3)=mdcrd(kp,3)
		p4(1)=mdcrd(lp,1); p4(2)=mdcrd(lp,2); p4(3)=mdcrd(lp,3)	
		call phipsiomg_angle(p1,p2,p3,p4,angle)
		do k=1, dihedral%multiply(j)
			potential=dihedral%pk(j,k)*(1+cosd(dihedral%pn(j,k)*angle-dihedral%phase(j,k)))
			dihedral_energy=dihedral_energy+potential
		enddo
	enddo

	do ik=1, group(ii,ic)%cnum2
		if(group(ii,ic)%atype2(ik)=="CB") goto 10
		do jj=1,repeated_unit
			do jc=1, gnum
				if(stage==0) then
					if(jj.eq.ii.and.jc.eq.ic) then
						i_id=group_para(ii,ic)%atomid2(ik)
						do jk=1, group(jj,jc)%cnum1
							j_id=group_para(jj,jc)%atomid1(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 20
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 30
								endif
							enddo
30							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(flag==1) then
								vdw=vdw/vdw14_coeff
							endif
							energy=energy+vdw
20							continue
						enddo				
						do jk=1, group(jj,jc)%cnum2
							j_id=group_para(jj,jc)%atomid2(jk)
							if(i_id.eq.j_id) goto 40					
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 40
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 50
								endif
							enddo
50							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(flag==1) then
								vdw=vdw/vdw14_coeff
							endif
							energy=energy+vdw				
40							continue
						enddo
						do jk=1, group(jj,jc)%cnum3
							j_id=group_para(jj,jc)%atomid3(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 60
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 70
								endif
							enddo
70							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(flag==1) then
								vdw=vdw/vdw14_coeff
							endif
							energy=energy+vdw
60							continue
						enddo
					else
						do jk=1, group(jj,jc)%cnum1
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							energy=energy+vdw
						enddo
						do jk=1, group(jj,jc)%cnum2
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							energy=energy+vdw					
						enddo
						do jk=1, group(jj,jc)%cnum3
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							energy=energy+vdw
						enddo
					endif

				elseif(stage==1) then
					if(jj.eq.ii.and.jc.eq.ic) then
						i_id=group_para(ii,ic)%atomid2(ik)
						do jk=1, group(jj,jc)%cnum1
							j_id=group_para(jj,jc)%atomid1(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 80
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 90
								endif
							enddo
90							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons1(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons1(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge1(jk))/(dielecons4solute*sqrt(rij))
							if(flag==1) then
								vdw=vdw/vdw14_coeff
								ele=ele/ele14_coeff
							endif
							energy=energy+vdw+ele
80							continue
						enddo				
						do jk=1, group(jj,jc)%cnum2
							j_id=group_para(jj,jc)%atomid2(jk)
							if(i_id.eq.j_id) goto 100					
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 100
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 110
								endif
							enddo
110							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons2(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons2(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge2(jk))/(dielecons4solute*sqrt(rij))
							if(flag==1) then
								vdw=vdw/vdw14_coeff
								ele=ele/ele14_coeff
							endif
							energy=energy+vdw+ele
100							continue
						enddo
						do jk=1, group(jj,jc)%cnum3
							j_id=group_para(jj,jc)%atomid3(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 120
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 130
								endif
							enddo
130							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons3(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons3(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge3(jk))/(dielecons4solute*sqrt(rij))
							if(flag==1) then
								vdw=vdw/vdw14_coeff
								ele=ele/ele14_coeff
							endif
							energy=energy+vdw+ele
120							continue
						enddo
					else
						do jk=1, group(jj,jc)%cnum1
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons1(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons1(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge1(jk))/(dielecons4solute*sqrt(rij))
							energy=energy+vdw+ele
						enddo
						do jk=1, group(jj,jc)%cnum2
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons2(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons2(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge2(jk))/(dielecons4solute*sqrt(rij))
							energy=energy+vdw+ele
						enddo
						do jk=1, group(jj,jc)%cnum3
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons3(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons3(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge3(jk))/(dielecons4solute*sqrt(rij))
							energy=energy+vdw+ele
						enddo			
					endif
				endif
			enddo
		enddo
10		continue
	enddo
	energy=energy+dihedral_weighting_factor*dihedral_energy
	
	return
	end  subroutine sidechain_energy


	subroutine bindingenergy_noentropy(group, group_para, numex, inb, numex4, inb4, score, binding_energy, score4hydration, Pagg)
	implicit none
	integer							:: ligan_recep_sta(num4category), ligan_recep_end(num4category)	
	integer							:: natom, ipres, atomid(repeated_unit*atom_num)
	integer    						:: numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)
	real							:: mdcrd(repeated_unit*atom_num,3), charge(repeated_unit*atom_num), epsion(repeated_unit*atom_num)
	real							:: r(repeated_unit*atom_num), rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)
	character*4						:: lbres(repeated_unit*atom_num)

	integer							:: categoryID, i, ii, ic, ik, j, k, flag, num4peptides
	real							:: rij, vdw, ele, sgb
	real							:: comp_vdw, comp_ele, comp_sgb, comp_snp
	real							:: ligan_recep_vdw(num4category), ligan_recep_ele(num4category), ligan_recep_sgb(num4category), ligan_recep_snp(num4category)
	real							:: binding_vdw, binding_ele, binding_sgb, binding_snp
	real							:: score, binding_energy, score4hydration, Pagg
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)
	real, dimension(:), allocatable :: alpha

	call hydrationscale(group, score4hydration)
    call aggregation_propensity(group, Pagg)

	natom=0
	flag=0
	num4peptides=0
	do categoryID=1, num4category
		num4peptides=num4peptides+selfassembly(categoryID)%num4peptides
		do i=1, selfassembly(categoryID)%num4peptides
			ii=selfassembly(categoryID)%peptideID(i)
			do ic=1, gnum			
				if(i==1.and.ic==1.and.flag==0) then
					ligan_recep_sta(categoryID)=natom+1
					flag=1
				endif

				ipres=natom
				do ik=1, group(ii,ic)%cnum1
					natom=natom+1
					charge(natom)=group_para(ii,ic)%charge1(ik)
					epsion(natom)=group_para(ii,ic)%epsion1(ik)
					r(natom)=group_para(ii,ic)%r1(ik)
					rborn(natom)=group_para(ii,ic)%rborn1(ik)
					fs(natom)=group_para(ii,ic)%fs1(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons1(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid1(ik)
					mdcrd(natom,1)=group(ii,ic)%coo1(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo1(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo1(ik,3)
					lbres(natom)=group(ii,ic)%gtype	
				enddo
				
				do ik=1, group(ii,ic)%cnum2
					natom=natom+1
					charge(natom)=group_para(ii,ic)%charge2(ik)
					epsion(natom)=group_para(ii,ic)%epsion2(ik)
					r(natom)=group_para(ii,ic)%r2(ik)
					rborn(natom)=group_para(ii,ic)%rborn2(ik)
					fs(natom)=group_para(ii,ic)%fs2(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons2(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid2(ik)
					mdcrd(natom,1)=group(ii,ic)%coo2(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo2(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo2(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo
	
				do ik=1, group(ii,ic)%cnum3
					natom=natom+1
					charge(natom)=group_para(ii,ic)%charge3(ik)
					epsion(natom)=group_para(ii,ic)%epsion3(ik)
					r(natom)=group_para(ii,ic)%r3(ik)
					rborn(natom)=group_para(ii,ic)%rborn3(ik)
					fs(natom)=group_para(ii,ic)%fs3(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons3(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid3(ik)
					mdcrd(natom,1)=group(ii,ic)%coo3(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo3(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo3(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo
				
				if(i==selfassembly(categoryID)%num4peptides.and.ic==gnum.and.flag==1) then
					ligan_recep_end(categoryID)=natom
					flag=0
				endif
				
			enddo
		enddo
	enddo

	allocate (alpha(natom))
	do i=1, natom
		call effgbradi(1,natom,i,rborn,fs,mdcrd,alpha(i))
	enddo

	comp_vdw=0.0; comp_ele=0.0;	comp_sgb=0.0; comp_snp=0.0
	do i=1, natom
		do j=i, natom
			vdw=0.0; ele=0.0; sgb=0.0
			rij=(mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
			if(j.ne.i) then
				do k=1,numex(atomid(i))
					if(atomid(j).eq.inb(atomid(i),k)) goto 40
				enddo
				flag=0
				do k=1,numex4(atomid(i))
					if(atomid(j).eq.inb4(atomid(i),k)) then
						flag=1
						goto 45
					endif
				enddo
45				continue
				call vdwcontri(i,j,rij,epsion,r,vdw)
				call elecontri(i,j,rij,charge,dielecons,ele)
				if(flag==1) then
					vdw=vdw/vdw14_coeff
					ele=ele/ele14_coeff
				endif
			endif
40			continue
			call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
			if(j.eq.i) sgb=sgb/2.0
			comp_vdw=comp_vdw+vdw
			comp_ele=comp_ele+ele
			comp_sgb=comp_sgb+sgb
		enddo
	enddo
	deallocate (alpha)

	do categoryID=1, num4category	
		allocate (alpha(natom))
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			call effgbradi(ligan_recep_sta(categoryID),ligan_recep_end(categoryID),i,rborn,fs,mdcrd,alpha(i))
		enddo

		ligan_recep_vdw(categoryID)=0.0; ligan_recep_ele(categoryID)=0.0; ligan_recep_sgb(categoryID)=0.0; ligan_recep_snp(categoryID)=0.0
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			do j=i, ligan_recep_end(categoryID)
				vdw=0.0; ele=0.0; sgb=0.0
				rij= (mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
				if(j.ne.i) then
					do k=1,numex(atomid(i))
						if(atomid(j).eq.inb(atomid(i),k)) goto 50
					enddo
					flag=0
					do k=1,numex4(atomid(i))
						if(atomid(j).eq.inb4(atomid(i),k)) then
							flag=1
							goto 55
						endif
					enddo
55					continue
					call vdwcontri(i,j,rij,epsion,r,vdw)
					call elecontri(i,j,rij,charge,dielecons,ele)
					if(flag==1) then
						vdw=vdw/vdw14_coeff
						ele=ele/ele14_coeff
					endif						
				endif
50				continue
				call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
				if(j.eq.i) sgb=sgb/2.0		
				ligan_recep_vdw(categoryID)=ligan_recep_vdw(categoryID)+vdw
				ligan_recep_ele(categoryID)=ligan_recep_ele(categoryID)+ele	
 				ligan_recep_sgb(categoryID)=ligan_recep_sgb(categoryID)+sgb
			enddo
		enddo
		deallocate (alpha)
		
	enddo
	
	binding_vdw=0; binding_ele=0; binding_sgb=0; binding_snp=0
	do categoryID=1, num4category
		binding_vdw=binding_vdw+ligan_recep_vdw(categoryID)
		binding_ele=binding_ele+ligan_recep_ele(categoryID)
		binding_sgb=binding_sgb+ligan_recep_sgb(categoryID)
!		binding_snp=binding_snp+ligan_recep_snp(categoryID)
	enddo

    binding_energy=((comp_vdw+comp_ele+comp_sgb+comp_snp)-(binding_vdw+binding_ele+binding_sgb+binding_snp))/real(num4peptides)
!	score=binding_energy+propensity_weighting_factor*(score4hydration+Pagg)
	score=binding_energy-propensity_weighting_factor*(score4hydration+Pagg)

	return
	end subroutine bindingenergy_noentropy

	
	subroutine entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
	implicit none
	integer							:: rotanum, obs, grade
	real							:: matrix(34,4)
	real							:: det, entropy
	character*4						:: aminoacid_name

	if(aminoacid_name=="GLY".or.aminoacid_name=="NGLY".or.aminoacid_name=="CGLY".or.aminoacid_name=="PRO".or.aminoacid_name=="NPRO".or.aminoacid_name=="CPRO".or.  &
	   aminoacid_name=="CYS".or.aminoacid_name=="NCYS".or.aminoacid_name=="CCYS".or.aminoacid_name=="ALA".or.aminoacid_name=="NALA".or.aminoacid_name=="CALA".or.  &
	   aminoacid_name=="VAL".or.aminoacid_name=="NVAL".or.aminoacid_name=="CVAL".or.aminoacid_name=="SER".or.aminoacid_name=="NSER".or.aminoacid_name=="CSER".or.  &
	   aminoacid_name=="THR".or.aminoacid_name=="NTHR".or.aminoacid_name=="CTHR") then
		entropy=0.0
	elseif(aminoacid_name=="PHE".or.aminoacid_name=="NPHE".or.aminoacid_name=="CPHE".or.aminoacid_name=="TYR".or.aminoacid_name=="NTYR".or.aminoacid_name=="CTYR") then
		if(obs.le.2) then
			entropy=-2.260
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.260
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="ARG".or.aminoacid_name=="NARG".or.aminoacid_name=="CARG") then
		if(obs.le.2) then
			entropy=-4.863
		else	
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-4.863
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="LYS".or.aminoacid_name=="NLYS".or.aminoacid_name=="CLYS") then
		if(obs.le.2) then
			entropy=-4.621
		else	
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-4.621
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="ASN".or.aminoacid_name=="NASN".or.aminoacid_name=="CASN") then
		if(obs.le.2) then
			entropy=-2.056
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.056
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif	
	elseif(aminoacid_name=="ASP".or.aminoacid_name=="NASP".or.aminoacid_name=="CASP") then
		if(obs.le.2) then
			entropy=-1.790
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-1.790
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="GLN".or.aminoacid_name=="NGLN".or.aminoacid_name=="CGLN") then
		if(obs.le.2) then
			entropy=-3.205
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-3.205
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif	
	elseif(aminoacid_name=="GLU".or.aminoacid_name=="NGLU".or.aminoacid_name=="CGLU") then
		if(obs.le.2) then
			entropy=-2.541
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.541
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="HIE".or.aminoacid_name=="NHIE".or.aminoacid_name=="CHIE") then
		if(obs.le.2) then
			entropy=-2.411
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.411
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="ILE".or.aminoacid_name=="NILE".or.aminoacid_name=="CILE") then
		if(obs.le.2) then
			entropy=-2.239
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.239
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif	
	elseif(aminoacid_name=="LEU".or.aminoacid_name=="NLEU".or.aminoacid_name=="CLEU") then
		if(obs.le.2) then
			entropy=-1.946
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-1.946
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="MET".or.aminoacid_name=="NMET".or.aminoacid_name=="CMET") then
		if(obs.le.2) then
			entropy=-3.495
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-3.495
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="TRP".or.aminoacid_name=="NTRP".or.aminoacid_name=="CTRP") then
		if(obs.le.2) then
			entropy=-2.316
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.316
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	endif
	
	return
	end subroutine entropy_calculation

	
	subroutine conformation_entropy(entropy4individual, entropy)
	implicit none
	integer							:: categoryID, chainID, ic, i, j
	real							:: entropy4individual(repeated_unit,gnum), entropy4sort(repeated_unit,gnum)
	real							:: minimum, entropy

	do ic=1, gnum
		chainID=1
		minimum=entropy4individual(chainID,ic)
		do chainID=2, repeated_unit
			if(entropy4individual(chainID,ic).lt.minimum) minimum=entropy4individual(chainID,ic)
		enddo

		do chainID=1, repeated_unit
			entropy4sort(chainID,ic)=minimum
		enddo
	enddo
	
	entropy=0.0
	chainID=1
!	do chainID=1, repeated_unit
		do ic=1,gnum	
			entropy=entropy+entropy4sort(chainID,ic)
		enddo
!	enddo

	return
	end subroutine conformation_entropy	

	
	subroutine bindingenergy(group, group_para, entropy4individual, numex, inb, numex4, inb4, score, binding_energy, entropy, score4hydration, Pagg)
	implicit none
	integer							:: ligan_recep_sta(num4category), ligan_recep_end(num4category)	
	integer							:: natom, ipres, atomid(repeated_unit*atom_num)
	integer    						:: numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)
	real							:: mdcrd(repeated_unit*atom_num,3), charge(repeated_unit*atom_num), epsion(repeated_unit*atom_num)
	real							:: r(repeated_unit*atom_num), rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)
	real							:: entropy4individual(repeated_unit,gnum)
	character*4						:: lbres(repeated_unit*atom_num)

	integer							:: categoryID, i, ii, ic, ik, j, k, flag, num4peptides
	real							:: rij, vdw, ele, sgb
	real							:: comp_vdw, comp_ele, comp_sgb, comp_snp
	real							:: ligan_recep_vdw(num4category), ligan_recep_ele(num4category), ligan_recep_sgb(num4category), ligan_recep_snp(num4category)
	real							:: binding_vdw, binding_ele, binding_sgb, binding_snp
	real							:: score, binding_energy, entropy, score4hydration, Pagg
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)
	real, dimension(:), allocatable :: alpha

	call conformation_entropy(entropy4individual, entropy)	
	call hydrationscale(group, score4hydration)
    call aggregation_propensity(group, Pagg)

	natom=0
	flag=0
	num4peptides=0
	do categoryID=1, num4category
		num4peptides=num4peptides+selfassembly(categoryID)%num4peptides
		do i=1, selfassembly(categoryID)%num4peptides
			ii=selfassembly(categoryID)%peptideID(i)
			do ic=1, gnum			
				if(i==1.and.ic==1.and.flag==0) then
					ligan_recep_sta(categoryID)=natom+1
					flag=1
				endif

				ipres=natom
				do ik=1, group(ii,ic)%cnum1
					natom=natom+1
					charge(natom)=group_para(ii,ic)%charge1(ik)
					epsion(natom)=group_para(ii,ic)%epsion1(ik)
					r(natom)=group_para(ii,ic)%r1(ik)
					rborn(natom)=group_para(ii,ic)%rborn1(ik)
					fs(natom)=group_para(ii,ic)%fs1(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons1(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid1(ik)
					mdcrd(natom,1)=group(ii,ic)%coo1(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo1(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo1(ik,3)
					lbres(natom)=group(ii,ic)%gtype	
				enddo
				
				do ik=1, group(ii,ic)%cnum2
					natom=natom+1
					charge(natom)=group_para(ii,ic)%charge2(ik)
					epsion(natom)=group_para(ii,ic)%epsion2(ik)
					r(natom)=group_para(ii,ic)%r2(ik)
					rborn(natom)=group_para(ii,ic)%rborn2(ik)
					fs(natom)=group_para(ii,ic)%fs2(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons2(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid2(ik)
					mdcrd(natom,1)=group(ii,ic)%coo2(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo2(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo2(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo
	
				do ik=1, group(ii,ic)%cnum3
					natom=natom+1
					charge(natom)=group_para(ii,ic)%charge3(ik)
					epsion(natom)=group_para(ii,ic)%epsion3(ik)
					r(natom)=group_para(ii,ic)%r3(ik)
					rborn(natom)=group_para(ii,ic)%rborn3(ik)
					fs(natom)=group_para(ii,ic)%fs3(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons3(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid3(ik)
					mdcrd(natom,1)=group(ii,ic)%coo3(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo3(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo3(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo

				if(i==selfassembly(categoryID)%num4peptides.and.ic==gnum.and.flag==1) then
					ligan_recep_end(categoryID)=natom
					flag=0
				endif

			enddo
		enddo
	enddo

	allocate (alpha(natom))
	do i=1, natom
		call effgbradi(1,natom,i,rborn,fs,mdcrd,alpha(i))
	enddo

	comp_vdw=0.0; comp_ele=0.0;	comp_sgb=0.0; comp_snp=0.0
	do i=1, natom
		do j=i, natom
			vdw=0.0; ele=0.0; sgb=0.0
			rij=(mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
			if(j.ne.i) then
				do k=1,numex(atomid(i))
					if(atomid(j).eq.inb(atomid(i),k)) goto 40
				enddo
				flag=0
				do k=1,numex4(atomid(i))
					if(atomid(j).eq.inb4(atomid(i),k)) then
						flag=1
						goto 45
					endif
				enddo
45				continue
				call vdwcontri(i,j,rij,epsion,r,vdw)
				call elecontri(i,j,rij,charge,dielecons,ele)
				if(flag==1) then
					vdw=vdw/vdw14_coeff
					ele=ele/ele14_coeff
				endif
			endif
40			continue
			call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
			if(j.eq.i) sgb=sgb/2.0
			comp_vdw=comp_vdw+vdw
			comp_ele=comp_ele+ele
			comp_sgb=comp_sgb+sgb
		enddo
	enddo
	deallocate (alpha)

	do categoryID=1, num4category	
		allocate (alpha(natom))
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			call effgbradi(ligan_recep_sta(categoryID),ligan_recep_end(categoryID),i,rborn,fs,mdcrd,alpha(i))
		enddo

		ligan_recep_vdw(categoryID)=0.0; ligan_recep_ele(categoryID)=0.0; ligan_recep_sgb(categoryID)=0.0; ligan_recep_snp(categoryID)=0.0
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			do j=i, ligan_recep_end(categoryID)
				vdw=0.0; ele=0.0; sgb=0.0
				rij= (mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
				if(j.ne.i) then
					do k=1,numex(atomid(i))
						if(atomid(j).eq.inb(atomid(i),k)) goto 50
					enddo
					flag=0
					do k=1,numex4(atomid(i))
						if(atomid(j).eq.inb4(atomid(i),k)) then
							flag=1
							goto 55
						endif
					enddo
55					continue
					call vdwcontri(i,j,rij,epsion,r,vdw)
					call elecontri(i,j,rij,charge,dielecons,ele)
					if(flag==1) then
						vdw=vdw/vdw14_coeff
						ele=ele/ele14_coeff
					endif						
				endif
50				continue
				call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
				if(j.eq.i) sgb=sgb/2.0		
				ligan_recep_vdw(categoryID)=ligan_recep_vdw(categoryID)+vdw
				ligan_recep_ele(categoryID)=ligan_recep_ele(categoryID)+ele	
 				ligan_recep_sgb(categoryID)=ligan_recep_sgb(categoryID)+sgb
			enddo
		enddo
		deallocate (alpha)
		
	enddo
	
	binding_vdw=0; binding_ele=0; binding_sgb=0; binding_snp=0
	do categoryID=1, num4category
		binding_vdw=binding_vdw+ligan_recep_vdw(categoryID)
		binding_ele=binding_ele+ligan_recep_ele(categoryID)
		binding_sgb=binding_sgb+ligan_recep_sgb(categoryID)
!		binding_snp=binding_snp+ligan_recep_snp(categoryID)
	enddo

    binding_energy=((comp_vdw+comp_ele+comp_sgb+comp_snp)-(binding_vdw+binding_ele+binding_sgb+binding_snp))/real(num4peptides)
!	score=binding_energy-entropy+propensity_weighting_factor*(score4hydration+Pagg)
	score=binding_energy-entropy-propensity_weighting_factor*(score4hydration+Pagg)

	return
	end subroutine bindingenergy

	
	subroutine hydrationscale(group, score4hydration)
	implicit none
	integer							:: i, j, ic, chainID
	integer							:: group_num(20)
	real							:: gly, leu, val, ile, met, phe, tyr, trp, glu, asp
	real							:: arg, lys, ser, thr, asn, gln, hie, ala, cys, pro
	real							:: avg_pho, avg_pol, avg_oth
	real							:: score_pho, score_pol, score_oth
	real							:: Npho, Npol, Noth
	real							:: Tscore, score4hydration
	type(groupdetails)				:: group(repeated_unit,gnum)

	ala=0.42; leu=2.32; val=1.66; ile=2.46; met=1.68; phe=2.44; tyr=1.31; trp=3.07; gly=0; glu=-0.87
	asp=-1.05; arg=-1.37; lys=-1.35; ser=-0.05; thr=0.35; asn=-0.82; gln=-0.30; hie=0.18; cys=1.34; pro=0.98
	avg_pho=1.58;  avg_pol=-0.63; avg_oth=1.10	
	score4hydration=0.0

	chainID=1
	group_num=0.0
	do ic=1,gnum
		if(group(chainID,ic)%gtype=="ALA".or.group(chainID,ic)%gtype=="NALA".or.group(chainID,ic)%gtype=="CALA") then
			group_num(1)=group_num(1)+1.0
		elseif(group(chainID,ic)%gtype=="LEU".or.group(chainID,ic)%gtype=="NLEU".or.group(chainID,ic)%gtype=="CLEU") then
			group_num(2)=group_num(2)+1.0
		elseif(group(chainID,ic)%gtype=="VAL".or.group(chainID,ic)%gtype=="NVAL".or.group(chainID,ic)%gtype=="CVAL") then
			group_num(3)=group_num(3)+1.0
		elseif(group(chainID,ic)%gtype=="ILE".or.group(chainID,ic)%gtype=="NILE".or.group(chainID,ic)%gtype=="CILE") then
			group_num(4)=group_num(4)+1.0
		elseif(group(chainID,ic)%gtype=="MET".or.group(chainID,ic)%gtype=="NMET".or.group(chainID,ic)%gtype=="CMET") then
			group_num(5)=group_num(5)+1.0
		elseif(group(chainID,ic)%gtype=="PHE".or.group(chainID,ic)%gtype=="NPHE".or.group(chainID,ic)%gtype=="CPHE") then
			group_num(6)=group_num(6)+1.0
		elseif(group(chainID,ic)%gtype=="TYR".or.group(chainID,ic)%gtype=="NTYR".or.group(chainID,ic)%gtype=="CTYR") then
			group_num(7)=group_num(7)+1.0
		elseif(group(chainID,ic)%gtype=="TRP".or.group(chainID,ic)%gtype=="NTRP".or.group(chainID,ic)%gtype=="CTRP") then
			group_num(8)=group_num(8)+1.0
		elseif(group(chainID,ic)%gtype=="GLY".or.group(chainID,ic)%gtype=="NGLY".or.group(chainID,ic)%gtype=="CGLY") then
			group_num(9)=group_num(9)+1.0
		elseif(group(chainID,ic)%gtype=="GLU".or.group(chainID,ic)%gtype=="NGLU".or.group(chainID,ic)%gtype=="CGLU") then
			group_num(10)=group_num(10)+1.0
		elseif(group(chainID,ic)%gtype=="ASP".or.group(chainID,ic)%gtype=="NASP".or.group(chainID,ic)%gtype=="CASP") then
			group_num(11)=group_num(11)+1.0
		elseif(group(chainID,ic)%gtype=="ARG".or.group(chainID,ic)%gtype=="NARG".or.group(chainID,ic)%gtype=="CARG") then
			group_num(12)=group_num(12)+1.0
		elseif(group(chainID,ic)%gtype=="LYS".or.group(chainID,ic)%gtype=="NLYS".or.group(chainID,ic)%gtype=="CLYS") then
			group_num(13)=group_num(13)+1.0
		elseif(group(chainID,ic)%gtype=="ASN".or.group(chainID,ic)%gtype=="NASN".or.group(chainID,ic)%gtype=="CASN") then
			group_num(14)=group_num(14)+1.0
		elseif(group(chainID,ic)%gtype=="GLN".or.group(chainID,ic)%gtype=="NGLN".or.group(chainID,ic)%gtype=="CGLN") then
			group_num(15)=group_num(15)+1.0
		elseif(group(chainID,ic)%gtype=="SER".or.group(chainID,ic)%gtype=="NSER".or.group(chainID,ic)%gtype=="CSER") then
			group_num(16)=group_num(16)+1.0
		elseif(group(chainID,ic)%gtype=="THR".or.group(chainID,ic)%gtype=="NTHR".or.group(chainID,ic)%gtype=="CTHR") then
			group_num(17)=group_num(17)+1.0
		elseif(group(chainID,ic)%gtype=="HIE".or.group(chainID,ic)%gtype=="NHIE".or.group(chainID,ic)%gtype=="CHIE") then
			group_num(18)=group_num(18)+1.0
		elseif(group(chainID,ic)%gtype=="CYS".or.group(chainID,ic)%gtype=="NCYS".or.group(chainID,ic)%gtype=="CCYS") then
			group_num(19)=group_num(19)+1.0
		elseif(group(chainID,ic)%gtype=="PRO".or.group(chainID,ic)%gtype=="NPRO".or.group(chainID,ic)%gtype=="CPRO") then
			group_num(20)=group_num(20)+1.0
		endif
	enddo

	score_pho=real(group_num(1)*ala+group_num(2)*leu+group_num(3)*val+group_num(4)*ile+group_num(5)*met+group_num(6)*phe+group_num(7)*tyr+group_num(8)*trp+group_num(9)*gly)
	score_pol=real(group_num(10)*glu+group_num(11)*asp+group_num(12)*arg+group_num(13)*lys+group_num(14)*asn+group_num(15)*gln+group_num(16)*ser+group_num(17)*thr+group_num(18)*hie)
	score_oth=real(group_num(19)*cys+group_num(20)*pro)
	
	Npho=group_num(1)+group_num(2)+group_num(3)+group_num(4)+group_num(5)+group_num(6)+group_num(7)+group_num(8)+group_num(9)
	Npol=group_num(10)+group_num(11)+group_num(12)+group_num(13)+group_num(14)+group_num(15)+group_num(16)+group_num(17)+group_num(18)
	Noth=group_num(19)+group_num(20)

	Tscore=abs(score_pho-avg_pho*Npho)+abs(score_pol-avg_pol*Npol)+abs(score_oth-avg_oth*Noth)
!	score4hydration=score4hydration+Tscore
	score4hydration=0-Tscore

	return
	end subroutine hydrationscale

	
	subroutine aggregation_propensity(group,Pagg)
	implicit none	
	integer							:: chainID, ic
	real							:: charge, helix, sheet, pat, pat_old
	real							:: CH, Alpha, Belta, Pattern
	real, parameter					:: coeff_CH=-0.16, coeff_Alpha=-5.7, coeff_Belta=5.0, coeff_Pat=0.39
	real							:: Pagg, p
	type(groupdetails)				:: group(repeated_unit,gnum)

	Pagg=0.0
	chainID=1
	CH=0.0; Alpha=0.0; Belta=0.0; Pattern=0.0; pat_old=0.5
	do ic=1,gnum
		if(group(chainID,ic)%gtype=="ALA".or.group(chainID,ic)%gtype=="NALA".or.group(chainID,ic)%gtype=="CALA") then
			charge=0.0;   helix=0.04;    sheet=0.12;    pat=1.0
		elseif(group(chainID,ic)%gtype=="LEU".or.group(chainID,ic)%gtype=="NLEU".or.group(chainID,ic)%gtype=="CLEU") then
			charge=0.0;   helix=0.38;    sheet=-0.15;   pat=1.0
		elseif(group(chainID,ic)%gtype=="VAL".or.group(chainID,ic)%gtype=="NVAL".or.group(chainID,ic)%gtype=="CVAL") then
			charge=0.0;   helix=0.06;    sheet=0.70;    pat=1.0
		elseif(group(chainID,ic)%gtype=="ILE".or.group(chainID,ic)%gtype=="NILE".or.group(chainID,ic)%gtype=="CILE") then
			charge=0.0;   helix=0.26;    sheet=0.77;    pat=1.0
		elseif(group(chainID,ic)%gtype=="MET".or.group(chainID,ic)%gtype=="NMET".or.group(chainID,ic)%gtype=="CMET") then
			charge=0.0;   helix=0.09;    sheet=0.71;    pat=1.0
		elseif(group(chainID,ic)%gtype=="PHE".or.group(chainID,ic)%gtype=="NPHE".or.group(chainID,ic)%gtype=="CPHE") then
			charge=0.0;   helix=0.01;    sheet=0.67;    pat=1.0
		elseif(group(chainID,ic)%gtype=="TYR".or.group(chainID,ic)%gtype=="NTYR".or.group(chainID,ic)%gtype=="CTYR") then
			charge=0.0;   helix=-0.05;   sheet=0.49;    pat=1.0
		elseif(group(chainID,ic)%gtype=="TRP".or.group(chainID,ic)%gtype=="NTRP".or.group(chainID,ic)%gtype=="CTRP") then
			charge=0.0;   helix=-0.21;   sheet=0.14;    pat=1.0
		elseif(group(chainID,ic)%gtype=="GLY".or.group(chainID,ic)%gtype=="NGLY".or.group(chainID,ic)%gtype=="CGLY") then
			charge=0.0;   helix=-1.24;   sheet=-0.76;   pat=0.0
		elseif(group(chainID,ic)%gtype=="GLU".or.group(chainID,ic)%gtype=="NGLU".or.group(chainID,ic)%gtype=="CGLU") then
			charge=-1.0;  helix=0.33;    sheet=-0.91;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="ASP".or.group(chainID,ic)%gtype=="NASP".or.group(chainID,ic)%gtype=="CASP") then
			charge=-1.0;  helix=-0.27;   sheet=-1.12;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="ARG".or.group(chainID,ic)%gtype=="NARG".or.group(chainID,ic)%gtype=="CARG") then
			charge=1.0;   helix=1.30;    sheet=-1.34;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="LYS".or.group(chainID,ic)%gtype=="NLYS".or.group(chainID,ic)%gtype=="CLYS") then
			charge=1.0;   helix=0.18;    sheet=-0.29;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="ASN".or.group(chainID,ic)%gtype=="NASN".or.group(chainID,ic)%gtype=="CASN") then
			charge=0.0;   helix=-0.25;   sheet=-1.05;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="GLN".or.group(chainID,ic)%gtype=="NGLN".or.group(chainID,ic)%gtype=="CGLN") then
			charge=0.0;   helix=0.02;    sheet=-1.67;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="SER".or.group(chainID,ic)%gtype=="NSER".or.group(chainID,ic)%gtype=="CSER") then
			charge=0.0;   helix=-0.15;   sheet=-1.45;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="THR".or.group(chainID,ic)%gtype=="NTHR".or.group(chainID,ic)%gtype=="CTHR") then
			charge=0.0;   helix=-0.39;   sheet=0.70;    pat=-1.0
		elseif(group(chainID,ic)%gtype=="HIE".or.group(chainID,ic)%gtype=="NHIE".or.group(chainID,ic)%gtype=="CHIE") then
			charge=0.0;   helix=0.11;    sheet=-1.34;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="CYS".or.group(chainID,ic)%gtype=="NCYS".or.group(chainID,ic)%gtype=="CCYS") then
			charge=0.0;   helix=-0.57;   sheet=0.63;    pat=1.0
		elseif(group(chainID,ic)%gtype=="PRO".or.group(chainID,ic)%gtype=="NPRO".or.group(chainID,ic)%gtype=="CPRO") then
			charge=0.0;   helix=0.0;     sheet=0.0;     pat=1.0
		endif

		CH=CH+charge
		Alpha=Alpha+helix
		Belta=Belta+sheet
		
		if(pat/=pat_old) then
			Pattern=Pattern+1.0
			pat_old=pat
		endif
	enddo
		
	p=coeff_CH*abs(CH)+coeff_Alpha*Alpha/gnum+coeff_Belta*Belta/gnum+coeff_Pat*Pattern/gnum
	Pagg=Pagg+p
	
	return
	end subroutine aggregation_propensity

	
	subroutine vdwcontri(x,y,rxy,epsion,r,vdw)
	implicit none
	integer					:: x, y
	real					:: rxy, vdw
	real					:: epsion_xy, r_xy
	real					:: acoeff, bcoeff
	real					:: epsion(repeated_unit*atom_num), r(repeated_unit*atom_num)

	epsion_xy=sqrt(epsion(x)*epsion(y))
	r_xy=r(x)+r(y)
	acoeff=epsion_xy*(r_xy**12)
	bcoeff=epsion_xy*2*(r_xy**6)
	vdw=acoeff/(rxy**6)-bcoeff/(rxy**3)

	return
	end subroutine vdwcontri

	
	subroutine	elecontri(x,y,rxy,charge,dielecons,ele)
	implicit none
	integer					:: x, y
	real					:: rxy,ele
	real					:: qx, qy, dielecons4solute
	real					:: charge(repeated_unit*atom_num),dielecons(repeated_unit*atom_num)

	qx=charge(x)
	qy=charge(y)
	if(dielecons(x).ge.dielecons(y)) then
		dielecons4solute=dielecons(x)
	else
		dielecons4solute=dielecons(y)
	endif
	ele=(qx*qy)/(dielecons4solute*sqrt(rxy))
	
	return
	end subroutine elecontri

	
	subroutine	sgbcontri(x,y,rxy,alphax,alphay,charge,dielecons,sgb)
	implicit none
	integer					:: x, y
	real					:: rxy, sgb
	real					:: dielecons4water
	real					:: fgb, alphax, alphay
	real					:: sgb1, sgb2
	real					:: qx, qy, dielecons4solute
	real					:: charge(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)

	dielecons4water=80.0
	qx=charge(x)
	qy=charge(y)
	fgb=sqrt(rxy+alphax*alphay*exp(-rxy/(4*alphax*alphay)))
	if(dielecons(x).ge.dielecons(y)) then
		dielecons4solute=dielecons(x)
	else
		dielecons4solute=dielecons(y)
	endif
	sgb1=(1.0/dielecons4solute)-(1.0/dielecons4water)
	sgb2=qx*qy/fgb
	sgb=-sgb1*sgb2

	return
	end subroutine sgbcontri

	
	subroutine effgbradi(xstart,xend,x,rborn,fs,mdcrd,alphax)
	implicit none
	integer					:: x
	integer					:: xstart, xend
	real					:: alpha, beta, gamma
	real					:: redborn
	real					:: psi, integra, alphax
	real					:: rborn(repeated_unit*atom_num),fs(repeated_unit*atom_num),mdcrd(repeated_unit*atom_num,3)

	alpha=0.8
	beta=0.0
	gamma=2.91
	redborn=rborn(x)-0.09
	call areafract(xstart,xend,x,rborn,fs,mdcrd,integra)
	psi=integra*redborn
	alphax=1.0/(1.0/redborn-tanh(alpha*psi-beta*psi*psi+gamma*psi*psi*psi)/rborn(x))

	return
	end subroutine effgbradi

	
	subroutine areafract(xstart,xend,x,rborn,fs,mdcrd,integra)
	implicit none
	integer					:: x, y
	integer					:: xstart, xend
	real					:: integra, rxy, sum, redborn
	real					:: lxy, uxy
	real					:: rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), mdcrd(repeated_unit*atom_num,3)

	integra=0.0
	do y=xstart, xend
		if(y.ne.x) then
			rxy=(mdcrd(x,1)-mdcrd(y,1))**2+(mdcrd(x,2)-mdcrd(y,2))**2+(mdcrd(x,3)-mdcrd(y,3))**2
			rxy=sqrt(rxy)
			redborn=fs(y)*(rborn(y)-0.09)
			call evalualij(x,rxy,redborn,rborn,lxy)
			call evaluauij(x,rxy,redborn,rborn,uxy)
			sum=(1.0/lxy)-(1.0/uxy)+(1.0/(uxy*uxy)-1.0/(lxy*lxy))*rxy/4.0+log(lxy/uxy)/(2.0*rxy)+ &
				(1.0/(lxy*lxy)-1.0/(uxy*uxy))*redborn*redborn/(4*rxy)
			integra=integra+sum
		endif
	enddo
	integra=integra/2.0

	return
	end	subroutine areafract

	
	subroutine evalualij(x,rxy,redborn,rborn,lxy)
	implicit none
	integer					:: x
	real					:: rxy, redborn
	real					:: lxy
	real					:: rborn(repeated_unit*atom_num)

	if(rborn(x).le.(rxy-redborn)) then
		lxy=rxy-redborn
	elseif(rborn(x).le.(rxy+redborn)) then
		lxy=rborn(x)-0.09
	else
		lxy=1.0
	endif

	return
	end subroutine evalualij

	
	subroutine evaluauij(x,rxy,redborn,rborn,uxy)
	implicit none
	integer					:: x
	real					:: rxy, redborn, redborn1
	real					:: uxy
	real					:: rborn(repeated_unit*atom_num)

	if(rborn(x).lt.(rxy+redborn)) then
		uxy=rxy+redborn
	else
		uxy=1.0
	endif

	return
	end subroutine evaluauij

end module energy_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module advancedfunction

	use constant
	use randomgenerator
	use input
	use mathfunction
	use database
	use transplant
	use energy_calculation

	contains
	subroutine scmf_substitution(group, sub_cycle, tgroup)
	implicit none
	integer							:: sub_cycle, i, j, k, l, m, ic, ic1, ic2, flag, flag1, flag2
	integer							:: ii, chainID, ip
	integer							:: Npho, Npol, Nchg, Noth
	integer							:: positive_number, negative_number, aminoacid_number, rotanum, feedback
	integer							:: icpointer(4,gnum), pos_type(gnum), neg_type(gnum)
	integer							:: group_num(4), num(4)
	integer							:: num4peptides, peptideID(6)
	real							:: ran2
	character*4						:: aminoacid_name(10), aminoacid_name_1, group_name_1(3), group_name_2(3)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum), temp_group_1(repeated_unit,gnum), temp_group_2(repeated_unit,gnum), aa_group(repeated_unit,40)

	tgroup=group
	group_num=0
	ii=1
	do ic=1, gnum
		if(tgroup(ii,ic)%gtype=="ALA".or.tgroup(ii,ic)%gtype=="LEU".or.tgroup(ii,ic)%gtype=="VAL".or.tgroup(ii,ic)%gtype=="ILE".or.tgroup(ii,ic)%gtype=="MET".or.tgroup(ii,ic)%gtype=="PHE".or.tgroup(ii,ic)%gtype=="TYR".or.tgroup(ii,ic)%gtype=="TRP".or.tgroup(ii,ic)%gtype=="GLY"  &
			.or.tgroup(ii,ic)%gtype=="NALA".or.tgroup(ii,ic)%gtype=="NLEU".or.tgroup(ii,ic)%gtype=="NVAL".or.tgroup(ii,ic)%gtype=="NILE".or.tgroup(ii,ic)%gtype=="NMET".or.tgroup(ii,ic)%gtype=="NPHE".or.tgroup(ii,ic)%gtype=="NTYR".or.tgroup(ii,ic)%gtype=="NTRP".or.tgroup(ii,ic)%gtype=="NGLY"  &
			.or.tgroup(ii,ic)%gtype=="CALA".or.tgroup(ii,ic)%gtype=="CLEU".or.tgroup(ii,ic)%gtype=="CVAL".or.tgroup(ii,ic)%gtype=="CILE".or.tgroup(ii,ic)%gtype=="CMET".or.tgroup(ii,ic)%gtype=="CPHE".or.tgroup(ii,ic)%gtype=="CTYR".or.tgroup(ii,ic)%gtype=="CTRP".or.tgroup(ii,ic)%gtype=="CGLY") then			   			   			   
			group_num(1)=group_num(1)+1
			icpointer(1, group_num(1))=ic
		elseif(tgroup(ii,ic)%gtype=="ASN".or.tgroup(ii,ic)%gtype=="GLN".or.tgroup(ii,ic)%gtype=="SER".or.tgroup(ii,ic)%gtype=="THR".or.tgroup(ii,ic)%gtype=="HIE" &
			.or.tgroup(ii,ic)%gtype=="NASN".or.tgroup(ii,ic)%gtype=="NGLN".or.tgroup(ii,ic)%gtype=="NSER".or.tgroup(ii,ic)%gtype=="NTHR".or.tgroup(ii,ic)%gtype=="NHIE" &
			.or.tgroup(ii,ic)%gtype=="CASN".or.tgroup(ii,ic)%gtype=="CGLN".or.tgroup(ii,ic)%gtype=="CSER".or.tgroup(ii,ic)%gtype=="CTHR".or.tgroup(ii,ic)%gtype=="CHIE") then			
			group_num(2)=group_num(2)+1
			icpointer(2, group_num(2))=ic		
		elseif(tgroup(ii,ic)%gtype=="ARG".or.tgroup(ii,ic)%gtype=="LYS".or.tgroup(ii,ic)%gtype=="GLU".or.tgroup(ii,ic)%gtype=="ASP" &
			.or.tgroup(ii,ic)%gtype=="NARG".or.tgroup(ii,ic)%gtype=="NLYS".or.tgroup(ii,ic)%gtype=="NGLU".or.tgroup(ii,ic)%gtype=="NASP" &
			.or.tgroup(ii,ic)%gtype=="CARG".or.tgroup(ii,ic)%gtype=="CLYS".or.tgroup(ii,ic)%gtype=="CGLU".or.tgroup(ii,ic)%gtype=="CASP") then			
			group_num(3)=group_num(3)+1
			icpointer(3, group_num(3))=ic
		elseif(tgroup(ii,ic)%gtype=="PRO".or.tgroup(ii,ic)%gtype=="CYS"  &
		    .or.tgroup(ii,ic)%gtype=="NPRO".or.tgroup(ii,ic)%gtype=="NCYS"  &
			.or.tgroup(ii,ic)%gtype=="CPRO".or.tgroup(ii,ic)%gtype=="CCYS") then
			group_num(4)=group_num(4)+1
			icpointer(4, group_num(4))=ic
		endif
	enddo
	
	Npho=fpho; Npol=fpol; Nchg=fchg; Noth=foth
	
	if((Npho+Npol+Nchg+Noth).ne.(group_num(1)+group_num(2)+group_num(3))) then
		open(20, file="error.txt", access="append")
			write(20,*) "Input setting:       ", "Npho=", Npho, "Npol=", Npol, "Nchg=", Nchg, "Noth=", Noth
			write(20,*) "Program calculation: ", "Npho=", group_num(1), "Npol=", group_num(2), "Nchg=", group_num(3), "Noth=", group_num(4)
			write(20,*) "Please adjust the number of residue in each category of residue types in the Input file!"
		close(20)
		stop	
	endif

	num(1)=group_num(1)-Npho
	num(2)=group_num(2)-Npol
	num(3)=group_num(3)-Nchg
	num(4)=group_num(4)-Noth

	positive_number=0
	negative_number=0
	do i=1, 4
		if(num(i)>0) then
			do j=1, (group_num(i)-1)
				call ran_gen(ran2,0)
				ic1=int(ran2*group_num(i)-1.0e-3)+1
				if(ic1.gt.group_num(i)) ic1=group_num(i)
				call ran_gen(ran2,0)
				ic2=int(ran2*group_num(i)-1.0e-3)+1
				if(ic2.gt.group_num(i)) ic2=group_num(i)

				k=icpointer(i,ic1)
				icpointer(i,ic1)=icpointer(i,ic2)
				icpointer(i,ic2)=k
			enddo
			do j=1, num(i)
				pos_type(positive_number+j)=i
			enddo
			positive_number=positive_number+num(i)		
		elseif(num(i)<0) then
			do j=1, abs(num(i))
				neg_type(negative_number+j)=i
			enddo
			negative_number=negative_number+abs(num(i))
		endif
	enddo

	do i=1, (positive_number-1)
		call ran_gen(ran2,0)
		ic1=int(ran2*positive_number-1.0e-3)+1
		if(ic1.gt.positive_number) ic1=positive_number
		call ran_gen(ran2,0)
		ic2=int(ran2*positive_number-1.0e-3)+1
		if(ic2.gt.positive_number) ic2=positive_number

		k=pos_type(ic1)
		pos_type(ic1)=pos_type(ic2)
		pos_type(ic2)=k
	enddo
	do i=1, (negative_number-1)
		call ran_gen(ran2,0)
		ic1=int(ran2*negative_number-1.0e-3)+1
		if(ic1.gt.negative_number) ic1=negative_number
		call ran_gen(ran2,0)
		ic2=int(ran2*negative_number-1.0e-3)+1
		if(ic2.gt.negative_number) ic2=negative_number

		k=neg_type(ic1)
		neg_type(ic1)=neg_type(ic2)
		neg_type(ic2)=k
	enddo

	temp_group_1=tgroup
	sub_cycle=min(positive_number, negative_number)
	do i=1, sub_cycle
		call scmf_choose_aminoacid(neg_type(i), aminoacid_number, aminoacid_name)
		l=1
		do while(l<=group_num(pos_type(i)))
			do j=1, aminoacid_number
				call groupinfo(temp_group_1(ii,icpointer(pos_type(i),l))%gtype, group_name_1, flag1)
				call groupinfo(aminoacid_name(j), group_name_2, flag2)

				aminoacid_name_1=group_name_2(flag1)

				call findrotamer(icpointer(pos_type(i),l), temp_group_1, aminoacid_name_1, rotanum, aa_group, ip)				

				flag=0
				chainID=1
				do m=1, rotanum
					call residue_replace(chainID, icpointer(pos_type(i),l), temp_group_1, m, aa_group, temp_group_2)

					call check_transplant(chainID, icpointer(pos_type(i),l), temp_group_2, feedback)
			
					if(feedback==1) then
						temp_group_1=temp_group_2
						flag=flag+1
						goto 10
					endif
				enddo
10				continue
				if(flag==1) goto 20				
			enddo
			l=l+1		
		enddo
20		continue

		if(flag==1) then
			do k=1, (group_num(pos_type(i))+l)
				if(k<=l) then
					icpointer(pos_type(i),(group_num(pos_type(i))+k))=icpointer(pos_type(i),k)
				else
					icpointer(pos_type(i),(k-l))=icpointer(pos_type(i),k)
				endif
			enddo
			group_num(pos_type(i))=group_num(pos_type(i))-1
		else
			open(20, file="error.txt", access="append")
				write(20,*) "scmf_substitution is wrong!"
				write(20,*) "Please check whether the code in the module MAINTENACE is right or not!"
			close(20)
			stop
		endif
	enddo

	tgroup=temp_group_1

	return 
	end subroutine scmf_substitution


	subroutine sidechain_optimization(stage, chainID, ic, group, group_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)
	implicit none
	integer								:: grade, grade_num(6), monitor(6)
	integer								:: chainID, i, j, k, ic, account_num, flag, stage, trial_count
	integer								:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60, 60)
	integer								:: dihedral_num	
	real								:: delta_chi, cos_angle, sin_angle, error, t, Tenergy_min, Tenergy
	real								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	real								:: h2_denominator, h3_denominator
	real								:: Dihedral4entropy(4)
	type(groupdetails)					:: group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	type(energyparameters)				:: group_para(repeated_unit,gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6), class_old(6), class_new(6), class_min(6)
	type(dihedralparameters)			:: dihedral	
	
	real, dimension(:), allocatable 	:: energy_forward, energy_backward
	real, dimension(:,:), allocatable 	:: gradient_old, gradient_new, Hessian_old, Hessian_new
	real, dimension(:,:), allocatable   :: H2, H3, H31
	real, dimension(:,:), allocatable	:: d, y, s, Tchi

	call sidechain_category(chainID, ic, group, Iclass, grade, grade_num, index, monitor)
	call dihedralangle_reading(group(chainID,ic)%gtype, dihedral_num, dihedral)

	allocate(energy_forward(grade)); allocate(energy_backward(grade))
	allocate(gradient_old(grade,1)); allocate(gradient_new(grade,1))
	allocate(Hessian_old(grade,grade)); allocate(Hessian_new(grade,grade))
	allocate(H2(grade,grade)); allocate(H3(grade,grade)); allocate(H31(grade,1))
	allocate(d(grade,1)); allocate(y(grade,1)); allocate(s(grade,1))

	Tgroup=group
	do i=1, Tgroup(chainID,ic)%cnum1
		if(Tgroup(chainID,ic)%atype1(i)=="CA") then
			CA(1)=Tgroup(chainID,ic)%coo1(i,1); CA(2)=Tgroup(chainID,ic)%coo1(i,2); CA(3)=Tgroup(chainID,ic)%coo1(i,3)
		endif
	enddo
	s=0.0
	class_new=Iclass

30	continue
	if(stage==0) then	
		delta_chi=5
	elseif(stage==1) then
		delta_chi=1
	endif
	cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		endif

		call axisrotation(rotaxis, cos_angle, sin_angle, m)
		
		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo
			
			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember
			
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000				
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo
		
		do j=1, Tgroup(chainID,ic)%cnum2
			Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_forward(i))
	enddo
	
	cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		endif

		call axisrotation(rotaxis, cos_angle, sin_angle, m)
		
		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo
			
			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember
			
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000				
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo
		
		do j=1, Tgroup(chainID,ic)%cnum2
			Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_backward(i))
	enddo
	
	do i=1, grade
		gradient_new(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)	
	enddo

	do i=1, grade
		do j=1, grade
			if(i==j) then
				Hessian_new(i,j)=1
			else
				Hessian_new(i,j)=0
			endif
		enddo
	enddo

	account_num=0
	t=0.0
	do while(.true.)
		Hessian_old=Hessian_new
		gradient_old=gradient_new
		d=-matmul(Hessian_old, gradient_old)
		
		allocate(Tchi(grade,1))
		trial_count=0
		do while(.true.)
			class_old=Iclass
			Tchi=t*d
			flag=0
			if(stage==0) then
				do i=1, grade
					if(Tchi(i,1).gt.5.0) then
						Tchi(i,1)=5.0
						flag=1
					elseif(Tchi(i,1).lt.(-5.0)) then
						Tchi(i,1)=-5.0
						flag=1
					endif
				enddo
			elseif(stage==1) then
				do i=1, grade
					if(Tchi(i,1).gt.1.0) then
						Tchi(i,1)=1.0
						flag=1
					elseif(Tchi(i,1).lt.(-1.0)) then
						Tchi(i,1)=-1.0
						flag=1
					endif
				enddo
			endif
			if(flag==1) account_num=account_num+1
			
			s=s+Tchi
			do i=1, grade
				cos_angle=cosd(s(i,1)); sin_angle=sind(s(i,1))
				if(i==1) then
					rotaxis_x=class_old(i)%member(monitor(i),1)-CA(1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-CA(2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-CA(3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				else
					rotaxis_x=class_old(i)%member(monitor(i),1)-class_old(i-1)%member(monitor(i-1),1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-class_old(i-1)%member(monitor(i-1),2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-class_old(i-1)%member(monitor(i-1),3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				endif

				call axisrotation(rotaxis, cos_angle, sin_angle, m)
					
				do j=(i+1), (grade+1)
					do k=1, grade_num(j)
						class_old(j)%member(k,1)=class_old(j)%member(k,1)-class_old(i)%member(monitor(i),1)
						class_old(j)%member(k,2)=class_old(j)%member(k,2)-class_old(i)%member(monitor(i),2)
						class_old(j)%member(k,3)=class_old(j)%member(k,3)-class_old(i)%member(monitor(i),3)
					enddo
						
					Tmember=matmul(class_old(j)%member, m)
					class_old(j)%member=Tmember
						
					do k=1, grade_num(j)
						class_old(j)%member(k,1)=anint((class_old(j)%member(k,1)+class_old(i)%member(monitor(i),1))*1000)/1000
						class_old(j)%member(k,2)=anint((class_old(j)%member(k,2)+class_old(i)%member(monitor(i),2))*1000)/1000				
						class_old(j)%member(k,3)=anint((class_old(j)%member(k,3)+class_old(i)%member(monitor(i),3))*1000)/1000
					enddo
				enddo
			enddo
				
			do j=1, Tgroup(chainID,ic)%cnum2
				Tgroup(chainID,ic)%coo2(j,1)=class_old(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(chainID,ic)%coo2(j,2)=class_old(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(chainID,ic)%coo2(j,3)=class_old(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, Tenergy)
			
			if(t==0.0) then
				Tenergy_min=Tenergy
				class_min=class_old
			else
				if(Tenergy.lt.Tenergy_min) then
					Tenergy_min=Tenergy
					class_min=class_old
					trial_count=trial_count+1
				else
					s=s-Tchi
					goto 10
				endif
			endif				
			t=delta_chi
		enddo
10		continue
		deallocate(Tchi)
		if(stage==0) then
			if(account_num.gt.20) goto 20
		elseif(stage==1) then
			if(account_num.gt.40) goto 20
		endif

		if(trial_count==0) goto 20
		error=0.0
		do i=1, grade
			error=error+abs(d(i,1))
		enddo
		if(error.lt.0.01) goto 20

		cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
		do i=1, grade
			Tclass=class_min
			if(i==1) then
				rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			else
				rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			endif

			call axisrotation(rotaxis, cos_angle, sin_angle, m)
			
			do j=(i+1), (grade+1)
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
					Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
					Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
				enddo
				
				Tmember=matmul(Tclass(j)%member, m)
				Tclass(j)%member=Tmember
				
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
					Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000				
					Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
				enddo
			enddo
			
			do j=1, Tgroup(chainID,ic)%cnum2
				Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_forward(i))
		enddo
		
		cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
		do i=1, grade
			Tclass=class_min
			if(i==1) then
				rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			else
				rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			endif

			call axisrotation(rotaxis, cos_angle, sin_angle, m)
			
			do j=(i+1), (grade+1)
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
					Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
					Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
				enddo
				
				Tmember=matmul(Tclass(j)%member, m)
				Tclass(j)%member=Tmember
				
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
					Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000				
					Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
				enddo
			enddo
			
			do j=1, Tgroup(chainID,ic)%cnum2
				Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_backward(i))
		enddo
		
		do i=1, grade
			gradient_new(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)			
		enddo
				
		y=gradient_new-gradient_old

		H2=matmul(s, transpose(s))
		h2_denominator=0.0
		do i=1, grade
			h2_denominator=h2_denominator+s(i,1)*y(i,1)
		enddo

		H3=matmul(Hessian_old, matmul(y, matmul(transpose(y), Hessian_old)))
		H31=matmul(Hessian_old, y)
		h3_denominator=0.0			
		do i=1, grade
			h3_denominator=h3_denominator+y(i,1)*H31(i,1)
		enddo			
		
		Hessian_new=Hessian_old+H2/h2_denominator-H3/h3_denominator
 
	enddo
20	continue
	class_new=class_min
	
	if(stage==0.and.Tenergy_min.lt.100.0) then
		stage=1
		goto 30
	endif

	Iclass=class_min

	deallocate(energy_forward); deallocate(energy_backward)
	deallocate(gradient_old); deallocate(gradient_new)
	deallocate(Hessian_old); deallocate(Hessian_new)
	deallocate(H2); deallocate(H3); deallocate(H31)	
	deallocate(d); deallocate(y); deallocate(s)
	
	do i=1, group(chainID,ic)%cnum2
		group(chainID,ic)%coo2(i,1)=Iclass(index(i)%class_No)%member(index(i)%member_No,1)
		group(chainID,ic)%coo2(i,2)=Iclass(index(i)%class_No)%member(index(i)%member_No,2)
		group(chainID,ic)%coo2(i,3)=Iclass(index(i)%class_No)%member(index(i)%member_No,3)
	enddo

	if(stage==1) then
		call torsionangle4sidechain(group, chainID, ic, grade, dihedral, Dihedral4entropy)
	endif	

	return
	end	subroutine sidechain_optimization

end module advancedfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module optimization_techniques

	use constant
	use pdbfile
	use mathfunction	
	use database
	use transplant
	use energy_calculation
	use advancedfunction	

	contains
	subroutine MC_technique(score_old, energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual)
	implicit none
	real							:: ran2, score_old, score_new, score_change
	real							:: energy_old, energy_new
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, Pagg_old
	real							:: score4hydration_new, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	score_change=score_new-score_old
	if(score_change.le.0) then
		group=tgroup
		score_old=score_new
		energy_old=energy_new
		entropy_old=entropy_new
		score4hydration_old=score4hydration_new
		Pagg_old=Pagg_new
		entropy4individual=Tentropy4individual
	else
		call ran_gen(ran2,0)
		if(ran2.le.exp(-score_change/ekt)) then
			group=tgroup
			score_old=score_new
			energy_old=energy_new
			entropy_old=entropy_new
			score4hydration_old=score4hydration_new
			Pagg_old=Pagg_new
			entropy4individual=Tentropy4individual
		endif
	endif
	
	return
	end subroutine MC_technique


	subroutine sequence_mutation_nonthermal(ic, group, aminoacid_name, tgroup, flag)
	implicit none
	integer							:: flag, i, j, m, m_best, rotanum, feedback
	integer							:: chainID, ic, ic_1, ic_2, stage, ip
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	real							:: ran2
	real							:: vdw_energy, vdw_energy_min, score, score_min
	real							:: binding_energy, entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para


	allocate(aa_group(repeated_unit,40))
	call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO") then	

		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		do chainID=1, repeated_unit
				
			m_best=0
			vdw_energy_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)

				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
				endif

				call check_transplant(chainID, ic, temp_group, feedback)

				if(feedback==1) then
					call vdwenergy(chainID, ic, temp_group, tgroup_para, vdw_energy)
					if(vdw_energy.lt.vdw_energy_min) then
						vdw_energy_min=vdw_energy
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo

			if(m_best==0) goto 10
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
10		continue
		deallocate(tgroup_para)
		deallocate(temp_group)

	else
		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		do chainID=1, repeated_unit

			m_best=0
			score_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)
				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
					call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
				endif
					
				stage=0
				call sidechain_optimization(stage, chainID, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

				if(stage==1) then
					call bindingenergy_noentropy(temp_group, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, binding_energy, score4hydration, Pagg)
					if(score.lt.score_min) then
						score_min=score
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo
				
			if(m_best==0) goto 20
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
20		continue
		deallocate(tgroup_para)
		deallocate(temp_group)
	endif
	deallocate(aa_group)

	return
	end subroutine sequence_mutation_nonthermal


	subroutine sequence_mutation(ic, group, entropy4individual, aminoacid_name, tgroup, Tentropy4individual, flag)
	implicit none
	integer							:: flag, i, j, m, m_best, feedback
	integer							:: categoryID, chainID, ic, ic_1, ic_2, ip, stage
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	integer							:: rotanum, obs, grade
	real							:: ran2
	real							:: vdw_energy, vdw_energy_min, score, score_min
	real							:: binding_energy, entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4), matrix(34,4)
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para

	allocate(aa_group(repeated_unit,40))
	call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
	grade=aa_lib(ip)%grade
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO") then	

		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		Tentropy4individual=entropy4individual
		do chainID=1, repeated_unit
				
			m_best=0
			vdw_energy_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)

				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
				endif

				call check_transplant(chainID, ic, temp_group, feedback)

				if(feedback==1) then
					call vdwenergy(chainID, ic, temp_group, tgroup_para, vdw_energy)
					if(vdw_energy.lt.vdw_energy_min) then
						vdw_energy_min=vdw_energy
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo

			if(m_best==0) goto 10						
			matrix=0.0; obs=4
			call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)			
			Tentropy4individual(chainID,ic)=entropy
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
10		continue
		deallocate(tgroup_para)
		deallocate(temp_group)

	else
		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		Tentropy4individual=entropy4individual
		do chainID=1, repeated_unit

			m_best=0
			score_min=500.0
			obs=0
			matrix=0.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)
				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
					call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
				endif

				stage=0
				call sidechain_optimization(stage, chainID, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

				if(stage==1) then
					obs=obs+1
					do j=1, grade
						matrix(obs,j)=Dihedral4entropy(j)
					enddo

					call bindingenergy_noentropy(temp_group, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, binding_energy, score4hydration, Pagg)
					if(score.lt.score_min) then
						score_min=score
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo

			if(m_best==0) goto 20			
			call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
			Tentropy4individual(chainID,ic)=entropy
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
20		continue
		deallocate(tgroup_para)
		deallocate(temp_group)
	endif
	deallocate(aa_group)

	return
	end subroutine sequence_mutation


!	subroutine scmf_loop(group)
!	implicit none
!	integer							:: ic, m, flag
!	integer							:: rotanum, chainID, ip, feedback
!	character*4						:: aminoacid_name
!	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum), temp_group(repeated_unit,gnum), aa_group(repeated_unit,40)

!	do ic=1, gnum
!		call mc_choose_aminoacid(ic, group, aminoacid_name)

!		call sequence_mutation_nonthermal(ic, group, aminoacid_name, tgroup, flag)
!		call findrotamer(ic, tgroup, aminoacid_name, rotanum, aa_group, ip)
		
!		flag=0
!		do chainID=1, repeated_unit
!			if(aminoacid_name==tgroup(chainID,ic)%gtype) then
!				flag=flag+1
!				goto 10
!			else
!				do m=1, rotanum
!					call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)
	
!					call check_transplant(chainID, ic, temp_group, feedback)
				
!					if(feedback==1) then
!						tgroup=temp_group
!						flag=flag+1
!						goto 10
!					endif
!				enddo
!			endif
!10			continue
!		enddo
		
!		if(flag==1) then
!			group=tgroup
!			goto 20
!		else
!			open(5, file="error.txt", access="append")
!				write(5,*) "scmf_loop runs wrong!"
!				write(5,*) "Errors take place in repacking the same residue on all the backbone scaffolds."
!			close(5)
!			stop			
!		endif
!20		continue	
!	enddo
!	group=tgroup
	
!	return
!	end subroutine scmf_loop


	subroutine initial_entropy_individual(group, entropy4individual)
	implicit none
	integer							:: i, j, ii, ip, m, stage
	integer							:: ic, rotanum, obs, grade
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	real							:: Dihedral4entropy(4), matrix(34,4)
	real							:: entropy,entropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name

	type(groupdetails)				:: group(repeated_unit,gnum)
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para

	entropy4individual=0.0
	do ii=1, repeated_unit
		do ic=1,gnum
			aminoacid_name=group(ii,ic)%gtype
			
			if(aminoacid_name=="GLY".or.aminoacid_name=="NGLY".or.aminoacid_name=="CGLY".or.aminoacid_name=="PRO".or.aminoacid_name=="NPRO".or.aminoacid_name=="CPRO".or.  &
			   aminoacid_name=="CYS".or.aminoacid_name=="NCYS".or.aminoacid_name=="CCYS".or.aminoacid_name=="ALA".or.aminoacid_name=="NALA".or.aminoacid_name=="CALA".or.  &
			   aminoacid_name=="VAL".or.aminoacid_name=="NVAL".or.aminoacid_name=="CVAL".or.aminoacid_name=="SER".or.aminoacid_name=="NSER".or.aminoacid_name=="CSER".or.  &
			   aminoacid_name=="THR".or.aminoacid_name=="NTHR".or.aminoacid_name=="CTHR") then

				rotanum=4; matrix=0.0; obs=4; grade=4
				call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
				entropy4individual(ii,ic)=entropy
			else
				allocate(aa_group(repeated_unit,40))
				call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
				grade=aa_lib(ip)%grade

				allocate(temp_group(repeated_unit,gnum))
				allocate(tgroup_para(repeated_unit,gnum))

				obs=0
				matrix=0.0
				do m=1, rotanum
					call residue_replace(ii, ic, group, m, aa_group, temp_group)
					if(m==1) then
						call energy_parameter(temp_group, tgroup_para)
						call atom_links4sidechain(ii, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
						call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
					endif

					stage=0
					call sidechain_optimization(stage, ii, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

					if(stage==1) then
						obs=obs+1
						do j=1, grade
							matrix(obs,j)=Dihedral4entropy(j)								
						enddo
					endif
				enddo

				call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
				entropy4individual(ii,ic)=entropy
		
				deallocate(tgroup_para)
				deallocate(temp_group)
				deallocate(aa_group)
			endif		
		enddo
	enddo

	return
	end subroutine initial_entropy_individual

	
	subroutine sequence_optimization_nonthermal(group)
	implicit none
	integer							:: attempt, ic_1, ic_2
	integer							:: i, ii, feedback_1, feedback_2, flag1, flag2
	real							:: ran2
	character*4						:: aminoacid_name_1, aminoacid_name_2
	character*4						:: group_name_1(3), group_name_2(3)

	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	do attempt=1, gnum
		call ran_gen(ran2,0)
		if(ran2.le.scmfswitch1) then
			call ran_gen(ran2,0)
			ic_1=int(ran2*gnum-1.0e-3)+1
			if(ic_1.gt.gnum) ic_1=gnum
			
			call mc_choose_aminoacid(ic_1, group, aminoacid_name_1)			

			call sequence_mutation_nonthermal(ic_1, group, aminoacid_name_1, tgroup, feedback_1)

			if(feedback_1==1) then
				group=tgroup
			endif

		else
			call ran_gen(ran2,0)
			ic_1=int(ran2*gnum-1.0e-3)+1
			if(ic_1.gt.gnum) ic_1=gnum

			do while(.true.)
				call ran_gen(ran2,0)
				ic_2=int(ran2*gnum-1.0e-3)+1
				if(ic_2.gt.gnum) ic_2=gnum
				if(ic_1.ne.ic_2) goto 20
			enddo
20			continue

			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			call groupinfo(group(1,ic_2)%gtype, group_name_2, flag2)

			aminoacid_name_1=group_name_2(flag1)
			aminoacid_name_2=group_name_1(flag2)
			
			call sequence_mutation_nonthermal(ic_1, group, aminoacid_name_1, temp_group, feedback_1)
	
			if(feedback_1==1) then
				call sequence_mutation_nonthermal(ic_2, temp_group, aminoacid_name_2, tgroup, feedback_2)
				if(feedback_2==1) then
					group=tgroup
				endif
			endif			

		endif
	enddo

	return
	end subroutine sequence_optimization_nonthermal	


	subroutine sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
	implicit none
	integer							:: step, attempt, Num4attempt, ic_1, ic_2, categoryID_1, categoryID_2, chainID_1, chainID_2
	integer							:: i, ii, j, feedback_1, feedback_2, flag1, flag2
	real							:: ran2
	real							:: score_old, score_new
	real							:: binding_energy_old, binding_energy_new
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, score4hydration_new
	real							:: Pagg_old, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum), temp_entropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name_1, aminoacid_name_2
	character*4						:: group_name_1(3), group_name_2(3)

	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(energyparameters), dimension(:,:), allocatable &
									:: group_para
	integer, dimension(:), allocatable &
									:: W_numex, W_numex4
	integer, dimension(:,:), allocatable &
									:: W_inb, W_inb4

	do attempt=1, 7
		call ran_gen(ran2,0)
		if(ran2.le.scmfswitch1) then
			call ran_gen(ran2,0)
			ic_1=int(ran2*gnum-1.0e-3)+1
			if(ic_1.gt.gnum) ic_1=gnum

			call mc_choose_aminoacid(ic_1, group, aminoacid_name_1)

			call sequence_mutation(ic_1, group, entropy4individual, aminoacid_name_1, tgroup, Tentropy4individual, feedback_1)

			if(feedback_1==1) then
				allocate(group_para(repeated_unit,gnum))
				allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
				allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))
				
				call atom_links(tgroup, W_numex, W_inb, W_numex4, W_inb4)
				call energy_parameter(tgroup, group_para)

				call bindingenergy(tgroup, group_para, Tentropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new)
				call MC_technique(score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual)

				deallocate(W_inb); deallocate(W_inb4)
				deallocate(W_numex); deallocate(W_numex4)
				deallocate(group_para)
			endif

			open(3, file="energydetails.txt",  access="append")
				write(3,4) step, attempt, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, "Not SCMF"
				write(3,"(<gnum>(a4))") (group(1,j)%gtype, j=1, gnum)
				write(3,*) "*******************************"
			close(3)

		else
			call ran_gen(ran2,0)
			ic_1=int(ran2*gnum-1.0e-3)+1
			if(ic_1.gt.gnum) ic_1=gnum

			do while(.true.)
				call ran_gen(ran2,0)
				ic_2=int(ran2*gnum-1.0e-3)+1
				if(ic_2.gt.gnum) ic_2=gnum
				if(ic_1.ne.ic_2) goto 20
			enddo
20			continue

			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			call groupinfo(group(1,ic_2)%gtype, group_name_2, flag2)

			aminoacid_name_1=group_name_2(flag1)
			aminoacid_name_2=group_name_1(flag2)
			
			call sequence_mutation(ic_1, group, entropy4individual, aminoacid_name_1, temp_group, temp_entropy4individual, feedback_1)
			
			if(feedback_1==1) then

				call sequence_mutation(ic_2, temp_group, temp_entropy4individual, aminoacid_name_2, tgroup, Tentropy4individual, feedback_2)

				if(feedback_2==1) then
					allocate(group_para(repeated_unit,gnum))
					allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
					allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))
					
					call atom_links(tgroup, W_numex, W_inb, W_numex4, W_inb4)
					call energy_parameter(tgroup, group_para)
		
					call bindingenergy(tgroup, group_para, Tentropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new)
					call MC_technique(score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual)

					deallocate(W_inb); deallocate(W_inb4)
					deallocate(W_numex); deallocate(W_numex4)
					deallocate(group_para)
				endif
			endif

			open(3, file="energydetails.txt",  access="append")
				write(3,4) step, attempt, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, "intra-SCMF"
				write(3,"(<gnum>(a4))") (group(1,j)%gtype, j=1, gnum)
				write(3,*) "*******************************"
			close(3)
		endif
		
		if(score_old.lt.energy_min(1)) then
			do i=(num4pdbsave-1),1,-1
				energy_min(i+1)=energy_min(i)
			enddo
			energy_min(1)=score_old
			open(2, file="minimum_energy.txt", access="append")
				write(2,*) "step=", step, "attempt=", attempt, "score=", score_old
			close(2)
			call generatepdb(step, attempt, group)
		elseif(score_old.lt.energy_min(num4pdbsave)) then
			do i=1,num4pdbsave
				if(score_old.eq.energy_min(i)) then
					goto 50
				elseif(score_old.lt.energy_min(i)) then
					do j=(num4pdbsave-1),i,-1
						energy_min(j+1)=energy_min(j)
					enddo
					energy_min(i)=score_old
					goto 60
				endif
			enddo
60			continue
			call generatepdb(step, attempt, group)
50			continue
		endif	
	enddo
4	format(i7,i7,5f20.13,a15)

	return
	end subroutine sequence_optimization

end module optimization_techniques							

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program ProteinDesign

	use constant
	use randomgenerator
	use input
	use pdbfile
	use mathfunction
	use database
	use energy_calculation
	use advancedfunction
	use optimization_techniques

	implicit none
	integer							:: step, sub_circle, i, j
	real							:: ran2
	real							:: score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old
	real							:: entropy4individual(repeated_unit,gnum)
	type(groupdetails)				:: group(repeated_unit,gnum)

	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group
	type(energyparameters), dimension(:,:), allocatable &
									:: group_para
	integer, dimension(:), allocatable &
									:: W_numex, W_numex4
	integer, dimension(:,:), allocatable &
									:: W_inb, W_inb4
	
	call inputfile
	call readpdb(group)
	call rotamerlib
	
	if(recalcu_switch==0) then
		allocate(temp_group(repeated_unit,gnum))
		call scmf_substitution(group, sub_circle, temp_group)
		if(sub_circle.ne.0) then
			group=temp_group
		endif
		deallocate(temp_group)

!		call scmf_loop(group)
		
		do step=1, 2
			call sequence_optimization_nonthermal(group)
		enddo

		call generatepdb(0, 0, group)
		energy_min=0.0

		allocate(group_para(repeated_unit,gnum))
		allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
		allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))

		call energy_parameter(group, group_para)
		call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
	
		call initial_entropy_individual(group, entropy4individual)
		call bindingenergy(group, group_para, entropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)

	
		deallocate(W_inb); deallocate(W_inb4)
		deallocate(W_numex); deallocate(W_numex4)
		deallocate(group_para)
	
		open(5, file="energyprofile.txt", access="append")
			write(5,6) 0, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old
		close(5)
		if(score_old.lt.0) energy_min(1)=score_old		
	
	else
		open(5, file="backup4backbone.txt", status="old")
			do i=1, num4pdbsave
				read(5,*) energy_min(i)
			enddo
			read(5,*) score_old
			read(5,*) binding_energy_old
			read(5,*) entropy_old
			read(5,*) score4hydration_old
			read(5,*) Pagg_old
			do i=1, repeated_unit
				do j=1, gnum
					read(5,*) entropy4individual(i,j)
				enddo
			enddo
		close(5)
	endif	
	
	do step=nstep_start, nstep_terminal
		call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)

		open(5, file="energyprofile.txt", access="append")
			write(5,6) step, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old
		close(5)
		
		if(mod(step,50).eq.0) then
			call generatepdb((step+1), 0, group) 
			call ran_gen(ran2,1)	
			open(5, file="backup4backbone.txt", status="replace")
				do i=1, num4pdbsave
					write(5,*) energy_min(i)
				enddo
				write(5,*) score_old
				write(5,*) binding_energy_old
				write(5,*) entropy_old
				write(5,*) score4hydration_old
				write(5,*) pagg_old
				do i=1, repeated_unit
					do j=1, gnum
						write(5,*) entropy4individual(i,j)
					enddo
				enddo
				write(5,*) "final step=", step
			close(5)
		endif
		
		
	enddo
6	format(i5,5f10.4)

end
	
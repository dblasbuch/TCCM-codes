 program MCLJNL
 Implicit none
 Real*8 :: T,sigma,L,dens,max,difr,skin,x,y,rn,dist,Tr,Dr,Pot,iPot,a,w,diff,e,Fin,Fin2,S1,R1,S2,R2,S3,R3                                         
 Integer :: Natom,rc,opt,prod,p,rmax,naccum,r100,M,i,j,k,at,zeros,Ndif,s,diff2,f,nmax,accept,c,val,n                                
 Real*8, Allocatable :: InitCoor(:,:),FinalCoor(:,:),Iter(:,:),NL(:,:),rdf(:,:),SV(:,:)
 Real*8, Parameter :: pi=acos(-1.0) !define the value of pi
 Integer, Allocatable :: Ni(:),ncount(:)
 Real*8, Allocatable :: g(:)
 Real*8 :: V(3)                                                     

 !------------------- File generation:

 open(10,file="result.out")               !results's file
 open(15,file="points.dat")               !initial points' file for gnuplot
 open( 20,file="setdata.plt")              !plot's 
 

 !------------------- Variables:

 T=1.2            !Temeprature in reduced units
 rc=3             !Cut off dsitance
 sigma=2.1        !Sigma for reduced units
 L=10.14          !Length of the cell in reduced units
 Natom=729            !Number of initial molecules
 opt=500*Natom !steps for optimization
 prod=125*Natom !steps for production
 p=opt+prod       !Number of steps of the montecarlo simulation
 dens=0.7         !density in reduced units
 max=0.2          !maximum movement in a change
 rmax=4           !max distnace for the rdf
 difr=0.01        !dr for the rdf
 skin=0.4         !rskin for neighbour-list procedure
 rn=rc+skin
 nmax=INT((dens*(4./3.)*pi*((rn)**3))*2)
 
 !------------------- Radial distrib. function variables:

 naccum=0               
 r100=INT(rmax/difr)
 allocate(Ni(r100))
 allocate(g(r100))
 allocate(rdf(r100,3))
 Ni=0
 Fin2=0
 
 !------------------- First structure generation (cubic):

 M=int(Natom**(1.d0/3.d0)) + 1                            !nearest integer of the sqrt3 of Natom
 
 M=M+1                                                !nearest integer +1
 Allocate(InitCoor(3,M**3))

 !loop that generates xyz coordinates from 0-1
 InitCoor=0
 at=1

 !Do i is for X and it changes everytime i changes
 !Do j is for Y and so same
 !Do k is for Z
 !First create 001 then 002 and so on up to 00M then j changes
 !So 010 011 and so on up to 0MM that then change i so 100 ... 10M ...1M0.... MMM
 !at is a counter of the number of the atom that is actually be placing
  
 do i = 1,M
  V(1)=(i-1)*L/M
  do j=1,M
   V(2)=(j-1)*L/M
   do k=1,M
    V(3)=(k-1)*L/M
    InitCoor(:,at) = V(:)
    at=at+1
   end do
  end do
 end do

 write(10,*) "InitCoor:"   !Write the generated coordinates in the file
 do i=1,M 
     write(10,*) InitCoor(i,:)
 end do
 
 !------------------- Delete extra random points:

 zeros=1
 do while (zeros < (M**3-Natom))
  call random_number(x)    
  At=int(x*M**3)+1
  if (sum(InitCoor(:,At)) /= 0) then   !avoid repetition
   InitCoor(:,At)= 0
   zeros = zeros+1
  end if
 end do

 Allocate(FinalCoor(3,Natom))

 !------------------- Obtain the final structure without "extra atoms":
 
 i=1
 do while(i<Natom)
  do j=1,M**3
   if (sum(InitCoor(:,j)) /= 0.d0 ) then
    FinalCoor(:,i)=InitCoor(:,j)
    i=i+1
   end if
  end do
 end do

 write (10,*) "FinalCoor:"   !write the final structure to initialize the MC procedure
 do i=1,Natom
  
     write(10,*) FinalCoor(i,:)
 end do

 !------------------- Monte Carlo (MC) program:
 accept=1
 s=0
 Allocate(Iter(1,3))
 Allocate(NL(nmax,Natom))
 Allocate(ncount(Natom))
 Allocate(SV(3,Natom))
 
 do while (s <= p)         !Loop for optimization and production steps

  write(6,*) "Calculating,step",s

  call random_number(x)
  Ndif = NINT((x*Natom)+1)            !+1 in order not to obtain 0 as Natom(random atom selection)
  do i=1,3
    Iter(1,i)=FinalCoor(i,Ndif)     !Change in structure to work with
  end do

  do i=1,3                         !Loop along the three coordinates
    call Random_number(x)          !Define the sign of the movement
    call Random_number(y)          !Value for random movement
    
    !Random modification:

    if (x>=0.5) then
        Iter(1,i)= Iter(1,i)+y*max    !situation for the ">L"
        if (Iter(1,i)>=L) then           !return it inside the cell
            Iter(1,i)=Iter(1,i)-L
        end if
    else
        Iter(1,i)= Iter(1,i)-y*max    !situation for the "<L"
        if (Iter(1,i)<=0) then           !return it inside the cell
            Iter(1,i)=Iter(1,i)+L        
        end if
    end if
  end do
  
  !Neighbour-list correction for efficiency:

  if (accept == 1) then
    if (s /= 0) then
      Fin2=0
      do i=1,3
        Fin=(Iter(1,i)-SV(i,Ndif))
        if (ABS(Fin) > L/2) then            !Periodic Boundary Conditions
          Fin=L-(ABS(Fin))
        end if
        Fin2=Fin2+(Fin**2)
      end do
    end if

      if (Fin2 < (skin/2)**2) then
        NL=0
        ncount=0
        SV=FinalCoor
        do i=1,Natom
          n=0
          do j=1,Natom
            dist=0
            do f=1,3
              e = FinalCoor(f,i)-FinalCoor(f,j)
              if (ABS(e) > L/2) then            !Periodic Boundary Conditions
                e = L-(ABS(e))
              end if
              dist = dist + ((e)**2)
            end do
            if (dist < (rn**2) .and. i/=j) then
              n=n+1 
              NL(n,i)=j         !atoms of the neighbour-list
            end if
          end do
          ncount(i)=n
        end do
      end if
  end if

    !Acumulation, used in order not to enter everytime in the production procedure (have better results in rdf)

    if ((s-opt)==25*naccum) then         !it doesn't actualize every step
        naccum=naccum+1
        do i=1,Natom
            do j=i+1,Natom
               dist=0
                do k=1,3
                    e = FinalCoor(k,i)-FinalCoor(k,j)
                    if (ABS(e) > L/2) then             !minimum image convetion
                        e = L- ABS(e)
                    end if
                    dist=dist+(e)**2
                end do
                 diff=sqrt(dist)
                if (diff<=rmax) then                  !condition to be close enough
                    diff2=INT(diff/difr)+1
                    Ni(diff2)=Ni(diff2)+1     !store the amount of times we have a distance
                end if
            end do
        end do
     end if

  !------------------- Lennard-Jones potential:
  
  iPot=0
  val=ncount(Ndif)        !number of neighbours of the moved atom
  
  do c=1,val
    do i=1,3                           !Periodic Boundary Conditions
      S1=(Iter(1,i)-FinalCoor(i,c))
      R1=(FinalCoor(i,Ndif)-FinalCoor(i,c))
      if (ABS(R1) > L/2) then
        R2=L-(ABS(R1))
      end if
      if (ABS(S1) > L/2) then            
        S2=L-(ABS(S1))
      end if
      S3=Tr+(S2**2)
      R3=Dr+(R2**2)
    end do
    Pot=((((1/(S3))**6)-(1/(S3)**3))-((1/(R3)**6)-(1/(R3)**3)))
    iPot=iPot+Pot
  end do
    
  iPot=4*iPot
  if (iPot > 0) then
    a=EXP(-iPot/T)                     !Accepted using alpha (a)
    call random_number(w)
    if (a >= w) then
        FinalCoor(:,Natom)=Iter(:,Natom)
        accept=1
    else
      accept=2
    end if
  else
    FinalCoor(:,Natom)=Iter(:,Natom)           !Automatically accepted
    accept=1
  end if
  
  s=s+1             !Counter for the number of steps

 end do

 Deallocate(InitCoor)
 Deallocate(FinalCoor)
 Deallocate(Iter)
 Deallocate(NL)
 
 !------------------- Write the final conformation obtained at the end of the program:
 write(6,*) Ni
 write(6,*) "Montecarlo finished"
 
 !------------------- Radial Distribution Function (RDF):
 
 do i=1,r100
  if (Ni(i) /= 0) then
      g(i)=3.d0*(Real(Ni(i)))/(Real(naccum)*(Real(Natom)/2)*4*pi*((((Real(i)*difr)+difr)**3)-((Real(i)*difr)**3))*dens)
  end if
  rdf(i,1)=Real(i)*difr
  rdf(i,2)=g(i)
  rdf(i,3)=Ni(i)
 end do

 do i=1,r100
   write(15,*) rdf(i,:)
 end do

 write(6,*) "The program has finished!"

 close(10)
 close(15)

 write(6,*) "Creation of the file with the information needed to execute gnuplot:"
 write(6,*) "-The data in order to print the obtained plot is in the points.dat file"
 write(6,*) "You just have to initialize gnuplot in your terminal and introduce the following commands:"
 write(6,*) "> plot 'points.dat' using 1:2 with line"

 end program MCLJNL
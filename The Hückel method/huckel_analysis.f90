!Author: Dani Blas Buch with some help

      subroutine btop_mat(filenm,n_carb,mat_top)   
      implicit double precision(a-h,o-z)
      dimension mat_top(200,200),dst_mat(200,200),rd_file(200,200), m_c_posit(200), rd_dst_mat(200,200)
      character*40 filenm, filetitle, atom_nm
      n_carb = 0      
      open(unit=9,file=filenm)
      read(9,*) natoms
      read(9,*) filetitle
      write(6,*) " "
      write(6,*) "             The molecule analyzed: ", filetitle

      do i = 1, 200
          read(9,*,end=25) atom_nm, (rd_file(i,j), j=1,3)
          if (atom_nm.eq."C") then           
              n_carb = n_carb+1        ! Counting number of carbon atoms 
              m_c_posit(n_carb) = i    ! Tracing positions of carbon atoms in the list of coordinates. 
          endif              
          length = i 
      enddo    
   25 continue 
      close(9)
      
      
      ! Building the distance matrix
      do i=1,length      ! Run across the reference atom 
          do j=1,length  ! Run across other atoms (including reference atom)
              dst_sq = 0
                  do k = 1,3      ! Run across x,y,z coordinates
                      crd_1 = rd_file(i,k)
                      crd_2 = rd_file(j,k)
                      dst_sq = dst_sq + (crd_1 - crd_2)**2
                  enddo
                  dst_mat(i,j) = dsqrt(dst_sq)    
          enddo
      enddo 
      
      ! Building reduced matrix from distance matrix (Using positions of carbon atoms in lists and matrices). 
      ind_1 = 1
      do i = 1,length
          if (any(m_c_posit==i)) then
              ind_2 = 1
              do j = 1, length
                  if (any(m_c_posit==j)) then
                      rd_dst_mat(ind_1,ind_2) = dst_mat(i,j)
                      ind_2 = ind_2 + 1
                  endif
              enddo
              ind_1 = ind_1 + 1 
          endif
      enddo    
      
      
      ! Building topological matrix from reduced distance matrix
      do i = 1, n_carb
          do j = 1, n_carb
              dst = rd_dst_mat(i,j)
              if ((dst.gt.1.20).and.(dst.lt.1.54)) then ! Using limiting cases as bounds (single and triple bonds) to decide if bonded. 
                  mat_top(i,j) = 1            ! Bonded atoms.  
              else
                  mat_top(i,j) = 0            ! Non-bonded atoms.
              endif
          enddo
      enddo            
                                                   
      end  
       
      
      
      subroutine bld_hamilt(n_carb, mat_top, alpha, beta, hamilt)

      implicit double precision(a-h,o-z)
      dimension mat_top(200,200), hamilt(200,200)
      do i = 1, n_carb
          do j = 1, n_carb
              if (i.eq.j) then
                  hamilt(i,j) = alpha
              else
                  if (mat_top(i,j).eq.1) then
                      hamilt(i,j) = beta
                  else
                      hamilt(i,j) = 0
                  endif
              endif
          enddo    
      enddo                
      end                        
                     
      
      subroutine diagonalize(n_carb,hamilt,d_hamilt,p_eig_nm)

      implicit double precision(a-h,o-z)
      dimension hamilt(200,200),d_hamilt(200,200),o_mat(200,200),o_mat_t(200,200),pr_int(200,200),p_eig(200,200),p_eig_nm(200,200)    
      
      val_max = 1
      conv = 1D-8
      iter = 0
      
      d_hamilt = hamilt        ! Copy hamiltonian into "diagonal" (not-yet) hamiltonian to be iterated
      do while (dabs(val_max).gt.conv) 
          val_max = 0
          iter = iter + 1
      
          do i = 1, n_carb        ! Find the largest non-diagonal element
              do j = 1, n_carb                        
                  if (i.ne.j) then
                      if (dabs(d_hamilt(i,j)).gt.dabs(val_max)) then
                          val_max = d_hamilt(i,j)
                          np = i
                          nq = j
                      endif
                  endif          
             enddo
          enddo      
          
         ! Start of creating transformation matrices
         ! write(6,*) "Maximum value is:", val_max
         ! write(6,*) " "
          
          theta = 0.5*datan(2*d_hamilt(np,nq)/(d_hamilt(np,np)-d_hamilt(nq,nq)))  ! Calculate angle for matrix setup          
          
          do i = 1, n_carb       ! Set up identity matrix first (some values to be substituted)
              do j = 1, n_carb
                  if (i.eq.j) then
                      o_mat(i,j) = 1
                  else 
                      o_mat(i,j) = 0
                  endif
              enddo
          enddo                  
                      
          o_mat(np,np) = dcos(theta)        ! Substitute appropriate values in created identity matrix. 
          o_mat(nq,nq) = dcos(theta)
          o_mat(np,nq) = -1*dsin(theta)
          o_mat(nq,np) = dsin(theta)
          
          o_mat_t = transpose(o_mat)
          
          ! End of creating transformation matrices
          
          pr_int = matmul(o_mat_t, d_hamilt)      ! Intermediate product of matrices 
          
          d_hamilt = matmul(pr_int,o_mat)        ! Similarity transformation O^T_pq x d-hamilt x O_pq completed. 
           
          ! Iterative computation of matrix containing eigenvectors in columns. 
                  
          if (iter.eq.1) then 
             p_eig = o_mat      ! Start of the iterative computation 
          else
             p_eig = matmul(p_eig,o_mat)    ! For second iter, older o_mat multiplies by new o_mat, and for n_iteration, product
          endif                             ! of all previous iteration is multiplied by new o_mat. 
       

          
      enddo
      ! Normalization of all eigenvectors
      
      do j = 1, n_carb
          v_norm_sq = 0                 ! Normalization constant squared 
          do i = 1, n_carb
              v_norm_sq = v_norm_sq + p_eig(i,j)**2   ! Accumulation of normalization constant squared accors rows
          enddo
          v_norm = dsqrt(v_norm_sq)     ! Calculation of normalization constant
          do i = 1, n_carb
              p_eig_nm(i,j) = p_eig(i,j)/v_norm      ! Normalization of eigenvector (in column j). 
          enddo
      enddo            
                  
      end            
                          
          
      subroutine order_eigs(n_carb,d_hamilt,p_eig_nm,eigvals,eigvecs)

      implicit double precision(a-h,o-z)
      dimension d_hamilt(200,200), eigvals(200), eigvecs(200,200), p_eig_nm(200,200), column(200)
      
      do i = 1,n_carb                  ! Copy diagonal elements into the 1D array. 
          eigvals(i) = d_hamilt(i,i)
      enddo
      
      eigvecs=p_eig_nm    

      
      do i = 1,n_carb                  ! Order by energies ascending 
          do k = i+1, n_carb
              if (eigvals(i).gt.eigvals(k)) then
                  r = eigvals(i)
                  eigvals(i) = eigvals(k)
                  eigvals(k) = r
                  column = eigvecs(:,i)
                  eigvecs(:,i) = eigvecs(:,k)
                  eigvecs(:,k) = column
                   
              endif
          enddo
      enddo 

      end
      

      subroutine occupancy(n_carb,n_pi_e,eigvals,occup)   ! Creates nx2 occupancy matrix (energy : n_elec)

      implicit double precision(a-h,o-z)
      dimension eigvals(200), occup(200,200)
      ! Figure out if number of electrons is odd or even,
      n_pi_e_av = n_pi_e
      level = 1
      do while(n_pi_e_av.gt.0)            ! Fill up the energy levels while there are electrons available
          n_deg = 0                       ! Initial definition of degeneracy to be counted
          eig = eigvals(level)            ! Get the energy of the level
          eig_l = eig - 0.001           ! Lower bound for interval in which energy for this level falls into. 
          eig_u = eig + 0.001           ! Upper bound for interval in which energy for this level falls into. 
          
          
          do i = 1, n_carb                   ! Count degeneracy of the energy level. 
              if ((eigvals(i).gt.eig_l).and.(eigvals(i).lt.eig_u)) then 
                  n_deg = n_deg + 1
              endif
          enddo
          
          n_cap = n_deg * 2               ! Capacity of the energy level for the electrons
          
          if (n_cap.lt.n_pi_e_av) then    
              n_el_lvl = n_cap            ! If there are more/equal electrons than level can accomodate, we fill the level completely.
          else
              n_el_lvl = n_pi_e_av        ! If there are less electrons available than level can accomodate, we use all the available electrons     
          endif 
          
          do i=level,level+n_deg-1          ! Filling all the degenerate orbitals within the level. 
              el_number = 1D0*n_el_lvl/n_deg  ! (This is to guarantee fractional occupancy in case of open-shell/SOMO systems). 
              occup(i,1) = eig             ! Eigenvalue (energy level) 
              occup(i,2) = el_number       ! (Potentially) fractional number of electrons on that orbital
          enddo
          
          n_pi_e_av = n_pi_e_av - n_el_lvl
          level = level + n_deg 
                                                                   
      enddo
      
      do i = level,n_carb
          occup(i,1) = eigvals(i)    
          occup(i,2) = 0
      enddo    
      end        
              
      ! Note that for rigorous Mulliken analysis, treating the cases of partial occupancies of level, in principle we might have
      ! fractional occupancies of levels. 
      
      

      subroutine an_mulliken(n_carb,eigvals,eigvecs,occup,popul)

      implicit double precision(a-h,o-z)
      dimension eigvals(200), eigvecs(200,200), occup(200,200), popul(200)
      do i = 1, n_carb
          pop_i = 0
          do j = 1, n_carb
              pop_i = pop_i + occup(j,2) * eigvecs(i,j)**2
          enddo
          popul(i) = pop_i
      enddo
      end    
      

      subroutine bond_order(n_carb,eigvecs,occup,bond_ord_mat, mat_top)

      implicit double precision(a-h,o-z)
      dimension eigvecs(200,200), occup(200,200), bond_ord_mat(200,200), mat_top(200,200)
      do j = 1, n_carb
          do k = 1, n_carb
              if (mat_top(j,k).eq.0) then
                   b_ord_jk = 0
              else 
                  b_ord_jk = 0
                  do i = 1, n_carb
                      b_ord_jk = b_ord_jk + occup(i,2)*eigvecs(j,i)*eigvecs(k,i)
                  enddo
              endif          
              bond_ord_mat(j,k) = b_ord_jk
          enddo
      enddo
      end
    
                  
      subroutine calc_energy(n_carb,occup,energy)
      implicit double precision(a-h,o-z)
      dimension occup(200,200)
      energy = 0
      do i = 1, n_carb
          energy = energy + occup(i,1) * occup(i,2)
      enddo
      end            
           
      
      
      implicit double precision(a-h,o-z)
      dimension mat_top(200,200),hamilt(200,200),d_hamilt(200,200),p_eig_nm(200,200),eigvals(200),eigvecs(200,200), &
      occup(200,200), popul(200), bond_ord_mat(200,200)
      character*40 filenm
      write(6,*) " "
      write(6,"(A)") "-----------------------------------------------------------------------------------------------------------"
      write(6,"(A)", advance="no") " Write name of the file of xyz file containing structure: "
      read (5,*) filenm
      write(6,"(A)", advance="no") " Write charge of the molecule: "
      read(5,*) n_charge
      
      alpha = -11.4      !(Units of eV) ------- Hückel parameter 
      beta = -0.8        !(Units of eV) ------- Hückel parameter 
      
      write(6,"(A)") "-----------------------------------------------------------------------------------------------------------" 
      
      call btop_mat(filenm,n_carb,mat_top)    
      call bld_hamilt(n_carb,mat_top,alpha,beta,hamilt)
      write(6,*) ""
      write(6,*) "             Number of carbon atoms:", n_carb
      n_pi_e = n_carb - n_charge
      write(6,*) "             Number of π electrons: ", n_pi_e

      call diagonalize(n_carb,hamilt,d_hamilt,p_eig_nm)
      
      call order_eigs(n_carb,d_hamilt,p_eig_nm,eigvals,eigvecs)  
      

      write(6,*) ""
      write(6,*) "" 
      write(6,"(A)", advance="no")"||---------------Ordered Eigenvalues (units of eV) and corresponding occupancies (N electrons)"
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"     

      write(6,"(A)", advance="no")"||--------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------" 	
      

      write(6,"(A)", advance="no")"||--------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------" 	               
      write(6,"(A, 80F12.4)") "|| Eigenvalue||", (eigvals(i),i=1,n_carb)           ! Printing ordered eigenvalues. 
      
      write(6,"(A)", advance="no")"||--------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------" 

      
      call occupancy(n_carb,n_pi_e,eigvals,occup)                         

      write(6,"(A, 80F12.4)") "|| Occupancy ||", (occup(i,2),i=1,n_carb)

      write(6,"(A)", advance="no")"||--------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------" 
      write(6,"(A)") ""  							
        
      write(6,"(A)") "" 
      
      write(6,"(A)", advance="no")"||--------------Ordered eigenvectors corresponding to order in eigenvalues (atom: coefficient)"
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"                    
      write(6,"(A)", advance="no")"||-----------||-------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"      
      write(6,"(A)", advance="no")"||   Atom #  ||----------------------------Coefficients---------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"    
      write(6,"(A)", advance="no")"||-----------||-------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"       
         
      do i=1,n_carb
          write(6,"(A,I3,A,80F12.4)") "||   ",i,"     ||", (eigvecs(i,j),j=1,n_carb)            
      enddo  
      write(6,"(A)", advance="no")"||-----------||-------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"     
      write(6,"(A)", advance="no")"||-----------||-------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"          
      write(6,*) "" 
      write(6,*) ""
      
     
      
      call an_mulliken(n_carb,eigvals,eigvecs,occup,popul)

      write(6,"(A)", advance="no")"----------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"      
      
      write(6,"(A)", advance="no")"---------------------------π-electron population on the atoms---------------------------------"
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"
             
      write(6,"(A)", advance="no")"----------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"      
            
      write(6,"(A, 80I12)") "||       Atom #      ||", (i, i=1,n_carb) 
      write(6,"(A)", advance="no")"||-------------------||-----------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"       
      write(6,"(A, 80F12.4)") "|| π-elec population ||", (popul(i),i=1,n_carb)
      write(6,"(A)", advance="no")"||-------------------||-----------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"       
       
 
                                            
      write(6,*) "" 
      write(6,*) "" 
      
      
      call bond_order(n_carb,eigvecs,occup,bond_ord_mat, mat_top)
      
             
      write(6,"(A)", advance="no")"----------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"  
      write(6,"(A)", advance="no")"---------------------------π-bond order for pairs of atoms -----------------------------------"
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"   
      write(6,"(A)", advance="no")"----------------------------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"        
      write(6,"(A, 80I12)") "||    Atom #   ||", (i, i=1,n_carb) 
      write(6,"(A)", advance="no")"||-------------||-----------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"  
      
      do i=1,n_carb
          write(6,"(A,I3,A,80F12.4)") "||    ",i,"      ||", (bond_ord_mat(i,j),j=1,n_carb)            
      enddo        
      write(6,"(A)", advance="no")"||-------------||-----------------------------------------------------------------------------" 
      write(6,"(A)")"------------------------------------------------------------------------------------------------------------"       
           
                                                 
      write(6,*) "" 

                                      
      write(6,*) "" 
      call calc_energy(n_carb,occup,energy)
      
      write(6,"(A46, F11.4, A3)") "Total π-electron energy of the system =", energy, "eV"
      write(6,*) " "      

  200 format(80f12.4)         
             
      end    

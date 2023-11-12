program main
  !$ use omp_lib
   implicit none
   integer,parameter:: md=200, nd = 200, od = 200   ! md, nd > grid size (m,n)
   real:: dx, dy, dz, dt
   real:: xnue, density, width, height, depth, time, inlet_velocity, outlet_pressure, AoA, thickness
   real,dimension(0:md,0:nd,0:od):: u, v, w, p, u_old, v_old, w_old 
   real,dimension(0:md,0:nd,0:od):: porosity
   real,dimension(0:md):: xp
   real,dimension(0:nd):: yp
   real,dimension(0:od):: zp
   real,dimension(0:md,0:nd,0:od)::ap, ae, aw, an, as, at, ab, bb
   integer:: m, n, o, istep, istep_max, iset, istep_out
   integer:: i, j, k
  !  real:: t1, t2, t3, t4
  !  real:: tp_start, tp_end, tu_start, tu_end, tv_start, tv_end, tw_start, tw_end
  !  real:: tbc_start, tbc_end, ttmpout_start, ttmpout_end
  !  real:: tp, tu, tv, tw, tbc, ttmpout
   character(len=50) :: output_folder
   character(len=50) :: csv_file
  ! ----------------
  ! read input data by using namelist 
  ! by Nobuto Nakamichi 4/7/2023
  namelist /file_control/istep_out
  namelist /grid_control/istep_max
  namelist /directory_control/output_folder, csv_file
  open(11,file="controlDict.txt",status="old",action="read")
  read(11,nml=file_control)
  read(11,nml=grid_control)
  read(11,nml=directory_control)
  close(11)
  ! ----------------
  ! call cpu_time(t1)
  ! write(*,*)'porosity setting:0 or calculation start:1 ?'
  ! read(*,*) iset
  ! make output directory
  
  call system('mkdir -p '//trim(output_folder))
  ! -----------------
  
  iset=1
  
  !-----------------
  ! porosity setting
  
  if (iset==0)then
    m=0         ! setup switch for grid conditions
    density=0.  ! setup switch for physical conditions
  
    call  physical_conditions (xnue, density, width, height, depth, time &
                          , inlet_velocity, outlet_pressure, AoA, m, n, o)
    call  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time &
                          , inlet_velocity, AoA, porosity, m, n, o, istep_max, iset)
    ! call  output_grid_list (xp, yp, m, n, angle_of_attack)
    stop
  end if
  
  ! ----------------
  ! calculation start  (if iest=!0)
  ! ----------------
  ! set up condistions
  m=0         ! setup switch for grid conditions
  density=0.  ! setup switch for physical conditions
  
  call  physical_conditions (xnue, density, width, height, depth, time &
                         , inlet_velocity, outlet_pressure, AoA, m, n, o)
  call  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time &
                        , inlet_velocity, AoA, porosity, m, n, o, istep_max, iset)
  call  output_grid (xp, yp, zp, m, n, o)
  
  ! write(*,*) "check", (porosity(i,10), i=1,m)
  
  istep = 0
  time = istep * dt
  
  ! ----------------
  
  write(*,*) 'istep_max= ', istep_max,'   istep_out= ', istep_out
  
  call  initial_conditions (p, u, v, w, xp, yp, zp, width, height, depth &
                         , inlet_velocity, outlet_pressure, AoA, m, n, o)
  call  boundary (p, u, v, w, xp, yp, zp, width, height, depth &
                         , inlet_velocity, outlet_pressure, AoA, porosity, m, n, o)
  
  ! print initial conditions
  ! call  output_solution (p, u, v, m, n)
  
  ! call cpu_time(t2)
  ! write(*,*) 'Initialization time= ', t2 - t1, '[s]'
  
  ! ----------------
  ! MAC algorithm start
  ! tp = 0.0
  ! tu = 0.0
  ! tv = 0.0
  ! tbc = 0.0
  ! ttmpout = 0.0
  do istep = 1, istep_max
  
    time=istep* dt
    write(*,*)'--- time_steps= ',istep, ' --  time = ',time
  
    do i = 0, m+1
      do j = 0, n+1
        do k = 0, o+1
          u_old(i,j,k) = u(i,j,k)
          v_old(i,j,k) = v(i,j,k)
          w_old(i,j,k) = w(i,j,k)
        end do
      end do
    end do
  
    ! call cpu_time(tp_start)
    call solve_p (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, height, thickness, yp, dx, dy, dz, dt, m, n, o)
    ! call cpu_time(tp_end)
    ! tp = tp + (tp_end - tp_start)
  
    ! call cpu_time(tu_start)
    call solve_u (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, dx, dy, dz, dt, m, n, o)
    ! call cpu_time(tu_end)
    ! tu = tu + (tu_end - tu_start)
  
    ! call cpu_time(tv_start)
    call solve_v (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, dx, dy, dz, dt, m, n, o)
    ! call cpu_time(tv_end)
    ! tv = tv + (tv_end - tv_start)

    ! call cpu_time(tw_start)
    call solve_w (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, dx, dy, dz, dt, m, n, o)
    ! call cpu_time(tw_end)
    ! tw = tw + (tw_end - tw_start)
  
    ! call cpu_time(tbc_start)
    call boundary(p, u, v, w, xp, yp, zp, width, height, depth  &
                      , inlet_velocity, outlet_pressure, AoA, porosity, m, n, o)
    ! call cpu_time(tbc_end)
    ! tbc = tbc + (tbc_end - tbc_start)
  
    ! call cpu_time(ttmpout_start)
    !call output_solution (p, u, v, m, n)
    if(mod(istep,istep_out)==0) call  output_paraview_temp (p, u, v, w, porosity, xp, yp, zp, m, n, o, istep)
    ! call cpu_time(ttmpout_end)
    ! ttmpout = ttmpout + (ttmpout_end - ttmpout_start)
  
  end do
  ! MAC algorithm end
  ! ----------------
  
  ! call cpu_time(t3)
  ! write(*,*) 'Pressure solver time= ', tp, '[s]'
  ! write(*,*) 'Velocity(1) solver time= ', tu, '[s]'
  ! write(*,*) 'Velocity(2) solver time= ', tv, '[s]'
  ! write(*,*) 'Velocity(3) solver time= ', tw, '[s]'
  ! write(*,*) 'Boundary condition time= ', tbc, '[s]'
  ! write(*,*) 'Total Mac algorithm time= ', t3 - t2, '[s]'
  ! print conditions (recall)
  ! call  physical_conditions (xnue, density, width, height, depth, time, inlet_velocity, outlet_pressure, m, n)
  ! call  grid_conditions (xp, yp, dx, dy, dt, xnue, density, width, height, depth, time, inlet_velocity, porosity, m, n, istep_max)
  
  ! print solutions 
  call  output_solution_post (p, u, v, w, xp, yp, zp, porosity, m, n, o)
  call  output_divergent (p, u, v, w, porosity, dx, dy, dz, m, n, o)
  call  output_paraview (p, u, v, w, porosity, xp, yp, zp, m, n, o)
  
  ! write(*,*) 'Total output time= ', ttmpout + (t4 -t3), '[s]'
  write(*,*) 'program finished'
  
  end program main
  !******************
  
  !  solve variables  
  
  !******************
  subroutine  solve_p (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, height, thickness, yp, dx, dy, dz, dt, m, n, o)
   implicit none
   integer,parameter:: md = 200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in):: dx, dy, dz, dt
   real,intent(in):: xnue, density, height, thickness
   real,intent(inout),dimension(0:md,0:nd,0:od):: u, v, w, p, u_old, v_old, w_old 
   real,intent(in),dimension(0:md,0:nd,0:od):: porosity
   real,intent(in),dimension(0:nd):: yp
   integer,intent(in):: m, n, o
  
  !-----------------
  ! local variables 
   real, parameter:: small = 1.e-6, big = 1.e6, zero = 0.
   real:: u_stg, v_stg
   real,dimension(0:md,0:nd,0:od):: ap, ae, aw, an, as, at, ab, bb, div
   integer:: i, j, k
   real:: fc, poro_grad
  
  !-----------------
  !  divergence term  div(u)
  !-----------------
  ! ----------------
  ! read input data by using namelist 
  ! by Nobuto Nakamichi 27/7/2023
  logical::nonslip
  namelist /calculation_method/nonslip
  open(11,file="controlDict.txt",status="old",action="read")
  read(11,nml=calculation_method)
  close(11)
  ! ----------------
  do i = 1, m
    do j = 1, n
      do k = 1, o
        div(i,j,k)= (u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5 &
                  + (v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5 &
                  + (w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5 
      end do
    end do
  end do
  
  do j = 1, n
    do k = 1, o
      div(0,j,k)  = 0.  ! inlet zx
      div(m+1,j,k)= 0.  ! outlet zx
    end do
  end do
  
  do i = 1, m
    do k = 1, o
      div(i,0,k)  = div(i,n,k)  ! periodic condition xz
      div(i,n+1,k)= div(i,1,k)
    end do
  end do

  do i = 1, m
    do j = 1, n
      div(i,j,0)  = div(i,j,o)  ! periodic condition xy
      div(i,j,o+1)= div(i,j,1)
    end do
  end do
  
  ! ----------------
  fc=0.
  
  do i = 1, m
  do j = 1, n
  do k = 1, o
  !poro_grad= sqrt( ((porosity(i+1,j)-porosity(i-1,j))/dx/2.)**2 &
  !                +((porosity(i,j+1)-porosity(i,j-1))/dy/2.)**2 )
  ! ----------------
  !   velocity u
  ! ----------------
  ! convection_x  (1st upwind scheme)
  !u(i,j)=u_old(i,j)						&
  !      -dt*max(u_old(i,j),0.)*(u_old(i,j)-u_old(i-1,j))/dx	&  ! u>0 1st upwind scheme 
  !      -dt*min(u_old(i,j),0.)*(u_old(i+1,j)-u_old(i,j))/dx	   ! u<0 1st upwind scheme
  !      -dt*u_old(i,j)*(u_old(i+1,j)-u_old(i-1,j))/dx*.5         ! 2nd central scheme (canceled)
  u(i,j,k)=u_old(i,j,k)-dt*( &
          fc*( max(u_old(i,j,k),0.)*(u_old(i,j,k)-u_old(i-1,j,k))/dx   &      ! u>0 1st upwind scheme 
              +min(u_old(i,j,k),0.)*(u_old(i+1,j,k)-u_old(i,j,k))/dx ) &      ! u<0 1st upwind scheme
  +(1.-fc)* u_old(i,j,k)*(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5)    ! 2nd central scheme
  
  ! convection_y
  !u(i,j)=u(i,j)							&
  !      -dt*max(v_old(i,j),0.)*(u_old(i,j)-u_old(i,j-1))/dy	&  ! v>0 1st upwind scheme 
  !      -dt*min(v_old(i,j),0.)*(u_old(i,j+1)-u_old(i,j))/dy	   ! v<0 1st upwind scheme
  !      -dt*v_old(i,j)*(u_old(i,j+1)-u_old(i,j-1))/dx*.5         ! 2nd central scheme (canceled)
  u(i,j,k)=u(i,j,k)-dt*( &
         fc*( max(v_old(i,j,k),0.)*(u_old(i,j,k)-u_old(i,j-1,k))/dy   &   ! v>0 1st upwind scheme 
              +min(v_old(i,j,k),0.)*(u_old(i,j+1,k)-u_old(i,j,k))/dy) &   ! v<0 1st upwind scheme
  +(1.-fc)* v_old(i,j,k)*(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5) ! 2nd central scheme

  ! convection_w
  u(i,j,k)=u(i,j,k)-dt*( &
          fc*( max(w_old(i,j,k),0.)*(u_old(i,j,k)-u_old(i,j,k-1))/dz   &      ! w>0 1st upwind scheme 
              +min(w_old(i,j,k),0.)*(u_old(i,j,k-1)-u_old(i,j,k))/dz ) &      ! w<0 1st upwind scheme
  +(1.-fc)* w_old(i,j,k)*(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5)    ! 2nd central scheme
  
  ! diffusion_x
  u(i,j,k)=u(i,j,k) +dt*xnue*(u_old(i+1,j,k)-2.*u_old(i,j,k)+u_old(i-1,j,k))/dx/dx 
  !      +dt*xnue/(small+porosity(i,j))*(u_old(i+1,j)-u_old(i-1,j))*(porosity(i+1,j)-porosity(i-1,j))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  u(i,j,k)=u(i,j,k) +dt*xnue*(u_old(i,j+1,k)-2.*u_old(i,j,k)+u_old(i,j-1,k))/dy/dy
  !      +dt*xnue/(small+porosity(i,j))*(u_old(i,j+1)-u_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! diffusion_z
  u(i,j,k)=u(i,j,k) +dt*xnue*(u_old(i,j,k+1)-2.*u_old(i,j,k)+u_old(i,j,k-1))/dz/dz
  
  ! divergence term
  u(i,j,k)=u(i,j,k) +dt*xnue*(3./3.)*(div(i+1,j,k)-div(i-1,j,k))/dx*.5

  ! additional terms by porosity profile
  u(i,j,k)=u(i,j,k)                 &
        +dt*( ( (u_old(i+1,j,k)-u_old(i-1,j,k))/dx*.5+(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*.5) &
                *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*.5                            &
             +( (u_old(i,j+1,k)-u_old(i,j-1,k))/dy*.5+(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*.5) &
                *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*.5                            &
             +( (u_old(i,j,k+1)-u_old(i,j,k-1))/dz*.5+(w_old(i+1,j,k)-w_old(i-1,j,k))/dx*.5) &
                *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*.5                            &
             + div(i,j,k)*xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5*(-0./3.)         &
         )/porosity(i,j,k)
  ! force on wall
  if (nonslip) then
    u(i,j,k)=u(i,j,k)- dt*xnue*u_old(i,j,k)/(thickness*dx)**2 *32.*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
  end if
  ! ----------------
  !   velocity v
  ! ----------------
  ! convection_x  (1st upwind scheme)
  !v(i,j)=v_old(i,j)						&
  !      -dt*max(u_old(i,j),0.)*(v_old(i,j)-v_old(i-1,j))/dx	&  ! u>0 1st upwind scheme
  !      -dt*min(u_old(i,j),0.)*(v_old(i+1,j)-v_old(i,j))/dx	   ! u<0 1st upwind scheme
  !      -dt*u_old(i,j)*(v_old(i+1,j)-v_old(i-1,j))/dx*.5         ! 2nd central scheme (canceled)
  v(i,j,k)=v_old(i,j,k)-dt*(          &
        fc *(max(u_old(i,j,k),0.)*(v_old(i,j,k)-v_old(i-1,j,k))/dx        &  ! u>0 1st upwind scheme
            +min(u_old(i,j,k),0.)*(v_old(i+1,j,k)-v_old(i,j,k))/dx)       &  ! u<0 1st upwind scheme
   +(1.-fc)* u_old(i,j,k)*(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5)      ! 2nd central scheme
  
  ! convection_y
  !v(i,j)=v(i,j)							&
  !      -dt*max(v_old(i,j),0.)*(v_old(i,j)-v_old(i,j-1))/dy	&  ! v>0 
  !      -dt*min(v_old(i,j),0.)*(v_old(i,j+1)-v_old(i,j))/dy	   ! v<0
  !      -dt*v_old(i,j)*(v_old(i,j+1)-v_old(i-1,j-1))/dx*.5         ! 2nd central scheme (canceled)
  v(i,j,k)=v(i,j,k)-dt*(          &
        fc *(max(v_old(i,j,k),0.)*(v_old(i,j,k)-v_old(i,j-1,k))/dy   &  ! v>0 
            +min(v_old(i,j,k),0.)*(v_old(i,j+1,k)-v_old(i,j,k))/dy)  &  ! v<0
   +(1.-fc)* v_old(i,j,k)*(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5)       ! 2nd central scheme

  ! convection_z
  v(i,j,k)=v(i,j,k)-dt*(          &
        fc *(max(w_old(i,j,k),0.)*(v_old(i,j,k)-v_old(i,j,k-1))/dz   &  ! v>0 
            +min(w_old(i,j,k),0.)*(v_old(i,j,k+1)-v_old(i,j,k))/dz)  &  ! v<0
   +(1.-fc)* w_old(i,j,k)*(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*0.5)       ! 2nd central scheme
  
  ! diffusion_x
  v(i,j,k)=v(i,j,k) +dt*xnue*(v_old(i+1,j,k)-2.*v_old(i,j,k)+v_old(i-1,j,k))/dx/dx
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i+1,j)-v_old(i-1,j))*(porosity(i+1,j)-porosity(i-1,j))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  v(i,j,k)=v(i,j,k) +dt*xnue*(v_old(i,j+1,k)-2.*v_old(i,j,k)+v_old(i,j-1,k))/dy/dy
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i,j+1)-v_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! diffusion_z
  v(i,j,k)=v(i,j,k) +dt*xnue*(v_old(i,j,k+1)-2.*v_old(i,j,k)+v_old(i,j,k-1))/dz/dz
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i,j+1)-v_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! divergence term   ! L+(2/3)N = (1/3)N;(2/3) or 0(1/3)
  v(i,j,k)=v(i,j,k) +dt*xnue*(3./3.)*(div(i,j+1,k)-div(i,j-1,k))/dy*.5
  ! additional terms by porosity profile   ! canceled for non-slip condition    ! L+(2/3)N = (1/3)N;(-1/3) or 0:(-2/3)N
  v(i,j,k)=v(i,j,k)               &
        +dt*( ( (v_old(i+1,j,k)-v_old(i-1,j,k))/dx*.5+(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*.5) &
                *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*.5                            &
             +( (v_old(i,j+1,k)-v_old(i,j-1,k))/dy*.5+(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*.5) &
                *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*.5                            &
             +( (v_old(i,j,k+1)-v_old(i,j,k-1))/dz*.5+(w_old(i,j+1,k)-w_old(i,j-1,k))/dy*.5) &
                *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*.5                            &
           + div(i,j,k)*xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5*(-0./3.)           &
         )/porosity(i,j,k)
  ! force on wall
  if (nonslip) then
    v(i,j,k)=v(i,j,k)- dt*xnue*v_old(i,j,k)/(thickness*dy)**2 *32.*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
  end if
  ! ----------------
  !   velocity w
  ! ----------------
  ! convection_x  (1st upwind scheme)
  w(i,j,k)=w_old(i,j,k)-dt*(          &
        fc *(max(u_old(i,j,k),0.)*(w_old(i,j,k)-w_old(i-1,j,k))/dx        &  ! w>0 1st upwind scheme
            +min(u_old(i,j,k),0.)*(w_old(i+1,j,k)-w_old(i,j,k))/dx)       &  ! w<0 1st upwind scheme
   +(1.-fc)* u_old(i,j,k)*(w_old(i+1,j,k)-w_old(i-1,j,k))/dx/2.)             ! 2nd central scheme
  
  ! convection_y
  w(i,j,k)=w(i,j,k)-dt*(          &
        fc *(max(v_old(i,j,k),0.)*(w_old(i,j,k)-w_old(i,j-1,k))/dy   &  ! w>0 
            +min(v_old(i,j,k),0.)*(w_old(i,j+1,k)-w_old(i,j,k))/dy)  &  ! w<0
   +(1.-fc)* v_old(i,j,k)*(w_old(i,j+1,k)-w_old(i,j-1,k))/dy/2.)        ! 2nd central scheme

  ! convection_z
  w(i,j,k)=w(i,j,k)-dt*(          &
        fc *(max(w_old(i,j,k),0.)*(w_old(i,j,k)-w_old(i,j,k-1))/dz   &  ! w>0 
            +min(w_old(i,j,k),0.)*(w_old(i,j,k+1)-w_old(i,j,k))/dz)  &  ! w<0
   +(1.-fc)* w_old(i,j,k)*(w_old(i,j,k+1)-w_old(i,j,k-1))/dz/2.   )     ! 2nd central scheme
  
  ! diffusion_x
  w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i+1,j,k)-2.*w_old(i,j,k)+w_old(i-1,j,k))/dx/dx
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i+1,j)-v_old(i-1,j))*(porosity(i+1,j)-porosity(i-1,j))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i,j+1,k)-2.*w_old(i,j,k)+w_old(i,j-1,k))/dy/dy
  !      +dt*xnue/(small+porosity(i,j))*(w_old(i,j+1)-w_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! diffusion_z
  w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i,j,k+1)-2.*w_old(i,j,k)+w_old(i,j,k-1))/dz/dz
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i,j+1)-v_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! divergence term   ! L+(2/3)N = (1/3)N;(2/3) or 0(1/3)
  w(i,j,k)=w(i,j,k) +dt*xnue*(3./3.)*(div(i,j,k+1)-div(i,j,k-1))/dz*.5
  ! additional terms by porosity profile   ! canceled for non-slip condition    ! L+(2/3)N = (1/3)N;(-1/3) or 0:(-2/3)N
  w(i,j,k)=w(i,j,k)               &
        +dt*( ( (w_old(i+1,j,k)-w_old(i-1,j,k))/dx*.5+(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*.5) &
                *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*.5                            &
             +( (w_old(i,j+1,k)-w_old(i,j-1,k))/dy*.5+(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*.5) &
                *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*.5                            &
             +( (w_old(i,j,k+1)-w_old(i,j,k-1))/dz*.5+(w_old(i,j,k+1)-w_old(i,j,k-1))/dz*.5) &
                *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*.5                            &
           + div(i,j,k)*xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5*(-0./3.)           &
         )/porosity(i,j,k)
  ! force on wall
  if (nonslip) then
    w(i,j,k)=w(i,j,k)- dt*xnue*w_old(i,j,k)/(thickness*dz)**2 *32.*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
  end if
  end do
  end do
  end do
  
  ! ----------------
  ! matrix solution  !  formulation of porous media
  
  do i = 1, m
  do j = 1, n
  do k = 1, o
  ae(i,j,k)= dt*max(small,(porosity(i+1,j,k)+porosity(i,j,k))*0.5)/dx/dx
  aw(i,j,k)= dt*max(small,(porosity(i,j,k)+porosity(i-1,j,k))*0.5)/dx/dx
  an(i,j,k)= dt*max(small,(porosity(i,j+1,k)+porosity(i,j,k))*0.5)/dy/dy
  as(i,j,k)= dt*max(small,(porosity(i,j,k)+porosity(i,j-1,k))*0.5)/dy/dy
  at(i,j,k)= dt*max(small,(porosity(i,j,k+1)+porosity(i,j,k))*0.5)/dz/dz
  ab(i,j,k)= dt*max(small,(porosity(i,j,k)+porosity(i,j,k-1))*0.5)/dz/dz
  ap(i,j,k)= -ae(i,j,k)-aw(i,j,k)-an(i,j,k)-as(i,j,k)-at(i,j,k)-ab(i,j,k)
  
  bb(i,j,k)= ((porosity(i+1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i+1,j,k))*0.5             &
             -(porosity(i-1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i-1,j,k))*0.5)*density/dx &
            +((porosity(i,j+1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j+1,k))*0.5             &
             -(porosity(i,j-1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j-1,k))*0.5)*density/dy &
            +((porosity(i,j,k+1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k+1))*0.5             &
             -(porosity(i,j,k-1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k-1))*0.5)*density/dz 
  
  !if (porosity(i,j) <small) then   !in solid (dummy solution)
  ! ap(i,j)=-1.
  ! bb(i,j)= 0.
  ! ae(i,j)= 0.25
  ! aw(i,j)= 0.25
  ! an(i,j)= 0.25
  ! as(i,j)= 0.25
  !end if
  
  end do
  end do
  end do
  
  call boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, o, height, yp)

  ! call solve_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, o)
  call solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, m, n, o)
  ! call solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, m, n, o)
  ! ----------------
  ! ----------------
  return
  end subroutine solve_p
  !******************
  
  !******************
  ! OpenMP Parallelized
  ! Written only for CPU machine
  ! No efficiency ensured on GPU machine 
  subroutine  solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, m, n, o)
   implicit none
   integer,parameter:: md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd,0:od):: p
   real,intent(in),dimension(0:md,0:nd,0:od):: ap, ae, aw, an, as, at, ab, bb
   integer,intent(in):: m, n, o
  
  ! local variables
  real:: relux_factor, error
  real,dimension(0:md,0:nd,0:od):: p_old
  integer::i, j, k, iter, iter_max, l

  !$omp parallel private(iter, i, j, k, l) &
  !$omp & shared(iter_max, relux_factor, m, n, o) &
  !$omp & shared(error, p_old, p, ap, ae, aw, an, as, at, ab, bb) &
  !$omp & default(none)
  
  ! ----------------
  !   SOR algorithm
  ! ----------------
  !$omp single
  iter_max = 100 ! SOR max interation steps
  relux_factor=1.7 ! SOR reluxation factor
  !$omp end single

  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'

  !$omp single 
  error=0.
  !$omp end single

  ! default periodic condition in yz-direction
  !$omp do
  do i = 1, m
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$omp end do

  ! default periodic condition in xy-direction
  !$omp do
  do i = 1, m
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp do
  do i = 0, m+1
    do j = 0, n+1
      do k = 0, o+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$omp end do
  
  !-- EVEN SPACE process
  !$omp do reduction(max:error)
  do l = 2, m*n*o, 2 ! evenspace
  k = (l - 1) / (m * n)  + 1
  j = ((l - 1) / m + 1) - (k - 1) * n
  i = (l - (j - 1) * m) - (k-1)*m*n

  !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
  if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i - 1
  p(i,j,k) = (bb(i,j,k) &
            - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
            - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
            - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
            / ap(i,j,k)*relux_factor &
            + p_old(i,j,k)*(1.-relux_factor)
  
  error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  !$omp end do

  ! default periodic condition in yz-direction
  !$omp do
  do i = 1, m
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$omp end do

  ! default periodic condition in xy-direction
  !$omp do
  do i = 1, m
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp do
  do i = 0, m+1
    do j = 0, n+1
      do k = 0, o+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$omp end do
  

  !-- ODD SPACE process
  !$omp do reduction(max:error)
  do l = 1, m*n*o, 2 ! odd space
    k = (l - 1) / (m * n)  + 1
    j = ((l - 1) / m + 1) - (k - 1) * n
    i = (l - (j - 1) * m) - (k-1)*m*n
  
    !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
    if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i + 1
    p(i,j,k) = (bb(i,j,k) &
              - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
              - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
              - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
              / ap(i,j,k)*relux_factor &
              + p_old(i,j,k)*(1.-relux_factor)
    
    error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  !$omp end do

  ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  
  end do

  ! default periodic condition in yz-direction
  !$omp do
  do i = 1, m
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$omp end do

  ! default periodic condition in xy-direction
  !$omp do
  do i = 1, m
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp master
  write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  !$omp end master
  !$omp end parallel

  ! if (error > 1e5) then
  !   write(*,*)'Error value diverges. Terminate the process.'
  !   call exit(0)
  ! end if
  ! ----------------
  
  return
  end subroutine solve_matrix_vec_omp
  !******************
 
  
  !******************
  ! OpenACC Parallelized
  ! Written only for GPU machine
  ! No efficiency ensured on CPU machine 
  subroutine  solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, m, n, o)
   implicit none
   integer,parameter:: md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd,0:od):: p
   real,intent(in),dimension(0:md,0:nd,0:od):: ap, ae, aw, an, as, at, ab, bb
   integer,intent(in):: m, n, o
  
  ! local variables
  real:: relux_factor, error
  real,dimension(0:md,0:nd,0:od):: p_old
  integer::i, j, k, iter, iter_max, l

  !$acc data copy(p_old, p, error) &
  !$acc & copyin(ap, ae, aw, an, as, at, ab, bb, relux_factor) 
  
  ! ----------------
  !   SOR algorithm
  ! ----------------

  iter_max = 100 ! SOR max interation steps
  relux_factor=1.7 ! SOR reluxation factor

  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'

  error=0.

  ! default periodic condition in yz-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$acc end kernels

  ! default periodic condition in xy-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  !$acc end kernels

  !$acc kernels
  !$acc loop independent
  do i = 0, m+1
  !$acc loop independent
    do j = 0, n+1
  !$acc loop independent
      do k = 0, o+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$acc end kernels
  
  !-- EVEN SPACE process
  !$acc kernels
  !$acc loop reduction(max:error)
  do l = 2, m*n*o, 2 ! evenspace
  k = (l - 1) / (m * n)  + 1
  j = ((l - 1) / m + 1) - (k - 1) * n
  i = (l - (j - 1) * m) - (k-1)*m*n

  !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
  if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i - 1
  p(i,j,k) = (bb(i,j,k) &
            - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
            - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
            - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
            / ap(i,j,k)*relux_factor &
            + p_old(i,j,k)*(1.-relux_factor)
  
  error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  !$acc end kernels

  ! default periodic condition in yz-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$acc end kernels

  ! default periodic condition in xy-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  !$acc end kernels

  !$acc kernels
  !$acc loop independent
  do i = 0, m+1
  !$acc loop independent
    do j = 0, n+1
  !$acc loop independent
      do k = 0, o+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$acc end kernels
  

  !-- ODD SPACE process
  !$acc kernels
  !$acc loop reduction(max:error)
  do l = 1, m*n*o, 2 ! odd space
    k = (l - 1) / (m * n)  + 1
    j = ((l - 1) / m + 1) - (k - 1) * n
    i = (l - (j - 1) * m) - (k-1)*m*n
  
    !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
    if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i + 1
    p(i,j,k) = (bb(i,j,k) &
              - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
              - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
              - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
              / ap(i,j,k)*relux_factor &
              + p_old(i,j,k)*(1.-relux_factor)
    
    error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  !$acc end kernels

  ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  
  end do

  ! default periodic condition in yz-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$acc end kernels

  ! default periodic condition in xy-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  !$acc end kernels

  !$acc end data 

  write(*,*)'SOR iteration no.', iter-1,'  -- error=', error

  ! if (error > 1e5) then
  !   write(*,*)'Error value diverges. Terminate the process.'
  !   call exit(0)
  ! end if
  ! ----------------
  
  return
  end subroutine solve_matrix_vec_oacc
  !******************
 
  !******************
  subroutine  solve_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, o)
   implicit none
   integer,parameter:: md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd,0:od):: p
   real,intent(in),dimension(0:md,0:nd,0:od):: ap, ae, aw, an, as, at, ab, bb
   integer,intent(in):: m, n, o
  
  ! local variables
  real:: relux_factor, error
  real,dimension(0:md,0:nd,0:od):: p_old
  integer::i, j, k, iter, iter_max
  
  ! ----------------
  !   SOR algorithm
  ! ----------------
  iter_max = 100 ! SOR max interation steps
  relux_factor=1.7 ! SOR reluxation factor
  
  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'
  error=0.
  
  ! default periodic condition in yz-direction
  do i = 1, m
    do k = 1, o
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do

  ! default periodic condition in xy-direction
  do i = 1, m
    do j = 1, n
      p(i,j,0)  =p(i,j,o)
      p(i,j,o+1)=p(i,j,1)
    end do
  end do
  
  do i = 0, m+1
    do j = 0, n+1
      do k = 0, o+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  
  do i = 1, m
  do j = 1, n
  do k = 1, o
  p(i,j,k) = (bb(i,j,k) &
            - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
            - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
            - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
            / ap(i,j,k)*relux_factor &
            + p_old(i,j,k)*(1.-relux_factor)
  
  error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  end do
  end do
  
  
  ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  
  end do
  
  write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  if (error > 1e5) then
    write(*,*)'Error value diverges. Terminate the process.'
    call exit(0)
  end if
  
  ! ----------------
  
  return
  end subroutine solve_matrix
  !******************
  
  !******************
  subroutine  boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, o, height, yp)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in)::height
   real,intent(in),dimension(0:md,0:nd,0:od)::p
   real,intent(inout),dimension(0:md,0:nd,0:od)::ap, ae, aw, an, as, at, ab, bb
   real,intent(in),dimension(0:nd)::yp
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  ! inlet (dp/x=0 at i=1)
  do j= 1, n
    do k = 1, o
      ae(1,j,k) =ae(1,j,k)+aw(1,j,k)
      aw(1,j,k) =0.
    end do
  end do
  
  ! outlet (p=outlet_pressure at i=m)
  do j= 1, n
    do k = 1, o
      bb(m,j,k) = bb(m,j,k)+ae(m,j,k)*p(m+1,j,k)
      ae(m,j,k) = 0.
      aw(m,j,k) = 0.
      an(m,j,k) = 0.
      as(m,j,k) = 0.
      at(m,j,k) = 0.
      ab(m,j,k) = 0.
    end do
  end do
  
  ! default : periodic condition in matrix solver
  
  ! symmetry or wall (dp/dy=0. at j=1)   xp>0
  !do i= 1,m
  ! an(i,1) =an(i,1)+as(i,1)
  ! as(i,1) = 0.
  !end do
  
  ! symmetry or wall  (dp/dy=0. at j=n)  xp>0
  !do i= 1,m
  ! as(i,n) =as(i,n)+an(i,n)
  ! an(i,n) = 0.
  !end do
  ! ----------------
  
  return
  end subroutine  boundrary_matrix
  !******************
  
  !******************
  subroutine  solve_u (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, dx, dy, dz, dt, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in)::dx, dy, dz, dt
   real,intent(in)::xnue, density
   real,intent(inout),dimension(0:md,0:nd,0:od)::u, v, w, p, u_old, v_old, w_old
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  do k = 1, o
  ! convection_x  (1st upwind scheme)
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)

  ! convection_z
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! diffusion_z
  ! (already calculated in solve_p)
  
  ! pressure
  u(i,j,k)=u(i,j,k) -dt/density*(p(i+1,j,k)-p(i-1,j,k))/dx*0.5
  end do
  end do
  end do
  
  ! ----------------
  return
  end subroutine solve_u
  !******************
  
  !******************
  subroutine  solve_v (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, dx, dy, dz, dt, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in)::dx, dy, dz, dt
   real,intent(in)::xnue, density
   real,intent(inout),dimension(0:md,0:nd,0:od)::u, v, w, p, u_old, v_old, w_old
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  do k = 1, o
  ! convection_x  (1st upwind scheme)
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)

  ! convection_z
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! diffusion_z
  ! (already calculated in solve_p)
  
  ! pressure
  v(i,j,k)=v(i,j,k) -dt/density*(p(i,j+1,k)-p(i,j-1,k))/dy*.5
  end do
  end do
  end do
  ! ----------------
  return
  end subroutine solve_v
  !******************

  !******************
  subroutine  solve_w (p, u, v, w, u_old, v_old, w_old, porosity, xnue, density, dx, dy, dz, dt, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in)::dx, dy, dz, dt
   real,intent(in)::xnue, density
   real,intent(inout),dimension(0:md,0:nd,0:od)::u, v, w, p, u_old, v_old, w_old
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  do k = 1, o
  ! convection_x  (1st upwind scheme)
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)

  ! convection_z
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! diffusion_z
  ! (already calculated in solve_p)

  ! pressure
  w(i,j,k)=w(i,j,k) -dt/density*(p(i,j,k+1)-p(i,j,k-1))/dz*.5
  end do
  end do
  end do
  ! ----------------
  return
  end subroutine solve_w
  !******************
  
  !  conditions  
  
  !******************
  subroutine  boundary(p, u, v, w, xp, yp, zp, width, height, depth    &
                       , inlet_velocity, outlet_pressure, AoA, porosity, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in)::width, height, depth, inlet_velocity, outlet_pressure, AoA
   real,intent(inout),dimension(0:md,0:nd,0:od)::u, v, w, p
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   integer,intent(in)::m, n, o
  
  ! local variables
   real, parameter::small=1.e-6, big=1.e6, zero=0., pai=atan(1.)*4.
   integer::i, j, k
  
  ! ----------------
  ! inlet (u=inlet_velocity, v=0., dp/dx=0 at i=1)
  do j = 1, n
    do k = 1, o
      u(1,j,k) =inlet_velocity*cos(AoA/180.*pai)
      v(1,j,k) =inlet_velocity*sin(AoA/180.*pai)
      w(1,j,k) =0.
      u(0,j,k) =u(1,j,k)    ! dummy
      v(0,j,k) =v(1,j,k)    ! dummy
      w(0,j,k) =w(1,j,k)    ! dummy
      p(0,j,k) =p(2,j,k)
    end do
  end do
  
  ! outlet (du/dx=0., dv/dx=0., p=outlet_pressure at i=m)
  do j = 1, n
    do k = 1, o
      u(m+1,j,k) =u(m-1,j,k)
      v(m+1,j,k) =v(m-1,j,k)
      w(m+1,j,k) =w(m-1,j,k)
      ! p(m,j) =outlet_pressure
      p(m+1,j,k)=outlet_pressure   ! dummy
    end do
  end do

  ! default: periodic condition (xz-direction at j=1 & n)
  do i= 0, m+1
    do k = 0, o+1
      u(i,0,k)   = u(i,n,k)
      v(i,0,k)   = v(i,n,k)
      w(i,0,k)   = w(i,n,k)
      p(i,0,k)   = p(i,n,k)
      u(i,n+1,k) = u(i,1,k)
      v(i,n+1,k) = v(i,1,k)
      w(i,n+1,k) = w(i,1,k)
      p(i,n+1,k) = p(i,1,k)
    end do
  end do

  ! default: periodic condition (xz-direction at j=1 & n)
  do i= 0, m+1
    do j = 0, n+1
      u(i,j,0)   = u(i,j,o)
      v(i,j,0)   = v(i,j,o)
      w(i,j,0)   = w(i,j,o)
      p(i,j,0)   = p(i,j,o)
      u(i,j,o+1) = u(i,j,1)
      v(i,j,o+1) = v(i,j,1)
      w(i,j,o+1) = w(i,j,1)
      p(i,j,o+1) = p(i,j,1)
    end do
  end do
  
  ! option: lower wall (u=0., v=0., dp/dy=0. at j=1)
  !do i= 0, m+1
  ! u(i,1) =0.
  ! v(i,1) =0.
  ! u(i,0) =0.					! dummy
  ! v(i,0) = -v(i,2)		  	! dummy
  ! p(i,0) =p(i,2)
  !end do
  
  ! option: symmetry (du/dy=0., v=0., dp/dy=0. at j=1)  xp>0
  !do i= 1, m
  ! u(i,0) = u(i,2)
  ! v(i,1) =0.
  ! v(i,0) = -v(i,2)		  	! dummy
  ! p(i,0) =p(i,2)
  !end do
  
  ! option: symmetry  (du/dy=0., v=0., dp/dy=0. at j=n)   xp>0
  !do i= 1, m
  ! u(i,n+1) = u(i,n-1)  
  ! v(i,n) =0.
  ! v(i,n+1) = -v(i,n-1)		! dummy 
  ! p(i,n+1) =p(i,n-1)
  !end do
  ! ----------------
  !do i= 0, m+1
  !do j= 0, n+1
  !if (porosity(i,j) <small) then   !in solid (dummy solution)
  ! u(i,j) = 0.
  ! v(i,j) = 0.
  !end if
  !end do
  !end do
  
  return
  end subroutine boundary
  !*****************************
  
  !*****************************
  subroutine physical_conditions(xnue, density, width, height, depth, time &
                                , inlet_velocity, outlet_pressure, AoA, m, n, o)
   implicit none
   integer,parameter:: md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(inout):: xnue, density, width, height, depth, time  &
                         ,inlet_velocity, outlet_pressure, AoA
   integer,intent(in):: m, n, o
  ! local variables
   real:: reynolds_no, wing_width
  !  integer:: i, j, k
  
  ! ----------------
  
  ! ----------------
  ! read input file 
  ! by Nobuto Nakamichi 4/7/2023
  namelist /physical/xnue, density, width, height, depth, time  &
                   ,inlet_velocity, outlet_pressure, AoA
  
  if (density == 0.)then
  
    open(11,file="controlDict.txt",status="old",action="read")
    read(11,nml=physical)
    close(11)
  
  end if
  
  ! ----------------
  
  !wing_width =1.    ! (m)
  !reynolds_no=wing_width*inlet_velocity/xnue
  
  write(*,*) 
  write(*,*) 'xnue =', xnue
  write(*,*) 'density =', density
  write(*,*) 'width =', width
  write(*,*) 'height =', height
  write(*,*) 'depth =', depth
  write(*,*) 'time =', time
  write(*,*) 'inlet_velocity =', inlet_velocity
  write(*,*) 'outlet_pressure =', outlet_pressure
  write(*,*) 'Angle of inlet_velocity (AoA) =', AoA
  !write(*,*) 'reynolds_no='	, reynolds_no
  
  ! ----------------
  
  return
  end subroutine physical_conditions
  !******************
  
  !******************
  subroutine  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time &
                              , inlet_velocity, AoA, porosity, m, n, o, istep_max, iset)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(inout)::dx, dy, dz, dt, AoA, thickness
   real,intent(in)::xnue, density, width, height, depth, time, inlet_velocity
   real,intent(inout),dimension(0:md,0:nd,0:od)::porosity
   real,intent(inout),dimension(0:md):: xp
   real,intent(inout),dimension(0:nd):: yp
   real,intent(inout),dimension(0:od):: zp
   integer,intent(inout):: m, n, o, istep_max, iset
   character(len = 50) :: csv_file
   character(len = 50) :: output_folder
  
  ! local variables
  !real,dimension(0:md,0:nd)::	distance
  real:: cfl_no, pecret_no, diffusion_factor, reynolds_no
  real:: pai, distance, center_x, center_y, radius
  integer:: i, j, k
  integer:: x, y, z
  real:: poro_val
  real, parameter:: small=1.e-6, big=1.e6, zero=0.
  ! --- 
  
  ! ----------------
  ! namelist 
  ! by Nobuto Nakamichi 4/7/2023
  namelist /grid_control/istep_max
  namelist /porosity_control/thickness
  namelist /directory_control/csv_file, output_folder
  open(11,file="controlDict.txt",status="old",action="read")
  read(11,nml=grid_control)
  read(11,nml=porosity_control)
  read(11,nml=directory_control)
  close(11)
  !-----------------
  
  ! read pixel data
  open(52,file=csv_file, form='formatted')
  
  read(52,*) m,n,o

  do k = 1, o
    do j = 1, n
      do i = 1, m
        read(52, *) x, y, z, poro_val
        porosity(x+1, y+1, z+1) = poro_val
      end do
    end do
  end do
  
  close(52) 
  
  ! thickness = 2.5
  dx = width / real(m-1)
  dy = height / real(n-1)
  dz = depth / real(o-1)
  dt = time / real(istep_max)
  
  radius = height*0.25
  
  cfl_no           = inlet_velocity * dt / dx
  pecret_no        = inlet_velocity * dx / xnue
  diffusion_factor = xnue * dt / dy / dy
  reynolds_no      = radius * inlet_velocity / xnue
  
  !----- check print out
  write(*,*)
  write(*,*) 'm, n, o =', m, n, o
  write(*,*) 'istep_max =', istep_max
  write(*,*) 'dx, dy, dz =', dx, dy, dz
  write(*,*) 'dt =', dt
  write(*,*) 'cfl_no =', cfl_no
  write(*,*) 'pecret_no =', pecret_no
  write(*,*) 'diffusion_factor =', diffusion_factor
  write(*,*) 'reynolds_no=', reynolds_no
  write(*,*) 'thickness =', thickness
  
  do i = 0, m+1
    xp(i) = dx * real(i-1) - width*0.5
  end do
  
  do j = 0, n+1
    yp(j) = dy * real(j-1) - height*0.5
  end do

  do k = 0, o+1
    zp(k) = dz * real(k-1) - depth*0.5
  end do
  
  ! set porosity (fluid area difinition) by distance function data
  !
  
  ! set porosity (fluid area difinition)
  ! pai=atan(1.0)*4.
  ! center_x=0.*width
  ! center_y=0.*height
  
  do i = 1, m
  do j = 1, n
  do k = 1, o
  porosity(i,j,k) = max(small, porosity(i,j,k))
  end do
  end do
  end do
  
  ! default: far field condtion in x-direction
  !do j = 1, n
  ! porosity(0,j)   = 1.
  ! porosity(m+1,j) = 1.
  !end do
  
  ! default: outlet condtion in x-direction
  do j = 1, n+1
  do k = 1, o+1
  porosity(0,j,k) = porosity(1,j,k)
  porosity(m+1,j,k) = porosity(m,j,k)
  end do
  end do
  
  ! default: periodic condtion in y-direction
  do i = 0, m+1
  do k = 0, o+1
  porosity(i,0,k)   = porosity(i,n,k)
  porosity(i,n+1,k) = porosity(i,1,k)
  end do
  end do
  
  ! ----------------
  return
  end subroutine  grid_conditions
  !******************
  
  !******************
  subroutine  initial_conditions (p, u, v, w, xp, yp, zp, width, height, depth  &
                                 , inlet_velocity, outlet_pressure, AoA, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in)::width, height, depth, inlet_velocity, outlet_pressure, AoA
   real,intent(out),dimension(0:md,0:nd,0:od)::u, v, w, p 
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  real, parameter :: pai=atan(1.)*4.   
  
  ! ----------------
  do j = 1, n
  do i = 1, m
  do k = 1, o
   u(i,j,k)=inlet_velocity*cos(AoA/180*pai)
   v(i,j,k)=inlet_velocity*sin(AoA/180*pai)
   w(i,j,k)=0.
   p(i,j,k)=outlet_pressure
  end do
  end do
  end do
  ! ----------------
  
  return
  end subroutine initial_conditions
  !******************
  
  ! output
  
  !******************
  subroutine  output_solution (p, u, v, w, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd,0:od)::u, v, w, p 
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  write(*,*)
  
  write(*,*)'velocity u '
  do k = 0, o+1
  do j = 0, n+1
  write(*,*) (u(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  
  write(*,*)'velocity v '
  do k = 0, o+1
  do j = 0, n+1
  write(*,*) (v(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  
  write(*,*)'velocity w '
  do k = 0, o+1
  do j = 0, n+1
  write(*,*) (w(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  
  write(*,*)'pressure'
  do k = 0, o+1
  do j = 0, n+1
  write(*,*) (p(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  ! ----------------
  
  return
  end subroutine output_solution
  !******************
  
  !******************
  subroutine  output_grid (xp, yp, zp, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  
  open (60, file='grid.dat', status='replace')
  ! ----------------
  write(60,*)'m, n, o =', m, n, o
  write(60,*)'grid points ='
  write(60,*) (xp(i), i=1,m)
  write(60,*) (yp(j), j=1,n)
  write(60,*) (zp(k), k=1,o)
  ! ----------------
  close (60)
  return
  end subroutine output_grid
  !******************
  
  !******************
  subroutine  output_grid_list (xp, yp, zp, m, n, o, angle_of_attack)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   integer,intent(in)::m, n, o
   real,intent(in):: angle_of_attack
  
  ! local variables
  integer::i, j
  real::pai=atan(1.)*4.
  real::x, y, z, th
  
  open (60, file='cellcenter.dat', status='replace')
  ! ----------------
  th = angle_of_attack/180.*pai
  do i=1,m
  do j=1,n
  x=xp(i)*cos(th)-yp(j)*sin(th)
  y=xp(i)*sin(th)+yp(j)*cos(th)
  z=zp(i)
  write(60,*) x,y,z
  end do
  end do
  ! ----------------
  close (60)
  return
  end subroutine output_grid_list
  !******************
  
  !******************
  subroutine  output_solution_post (p, u, v, w, xp, yp, zp, porosity, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd,0:od)::u, v, w, p
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   integer,intent(in)::m, n, o
  
  ! local variables
   real, parameter::small=1.e-6, big=1.e6, zero=0.
   real, parameter::pmin=0.25, pmax=0.75
   integer::i, j, k
   real,dimension(0:md, 0:nd, 0:od)::u_cnt, v_cnt, w_cnt, p_cnt
  
  open (61, file='solution_uvp.dat', status='replace')
  
  ! ----------------
  ! interpolation at p-center grid
  
  do i = 1, m
  do j = 1, n
  do k = 1, o
   u_cnt(i,j,k)=u(i,j,k)*porosity(i,j,k)
   v_cnt(i,j,k)=v(i,j,k)*porosity(i,j,k)
   w_cnt(i,j,k)=w(i,j,k)*porosity(i,j,k)
   if (porosity(i,j,k) > small) then
    p_cnt(i,j,k)=p(i,j,k)
   else
    p_cnt(i,j,k)=zero
   end if 
  end do
  end do
  end do
  
  do j = 1, n
  do k = 1, o
   u_cnt(0,j,k)=u_cnt(1,j,k)
   v_cnt(0,j,k)=v_cnt(1,j,k)
   w_cnt(0,j,k)=w_cnt(1,j,k)
   p_cnt(0,j,k)=p_cnt(1,j,k)
   u_cnt(m+1,j,k)=u_cnt(m,j,k)
   v_cnt(m+1,j,k)=v_cnt(m,j,k)
   w_cnt(m+1,j,k)=w_cnt(m,j,k)
   p_cnt(m+1,j,k)=p_cnt(m,j,k)
  end do
  end do
  
  do i = 0, m+1
   u_cnt(i,0,k)=u_cnt(i,1,k)
   v_cnt(i,0,k)=v_cnt(i,1,k)
   w_cnt(i,0,k)=w_cnt(i,1,k)
   p_cnt(i,0,k)=p_cnt(i,1,k)
   u_cnt(i,n+1,k)=u_cnt(i,n,k)
   v_cnt(i,n+1,k)=v_cnt(i,n,k)
   w_cnt(i,n+1,k)=w_cnt(i,n,k)
   p_cnt(i,n+1,k)=p_cnt(i,n,k)
  end do
  
  !-----------------
  write(61,*)'m, n, o =', m, n, o
  
  write(61,*)'velocity u_bulk '
  do k = 1, o
  write(61,*) ((u_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity v_bulk '
  do k = 1, o
  write(61,*) ((v_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity w_bulk '
  do k = 1, o
  write(61,*) ((w_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity u_inst '
  do k = 1, o
  write(61,*) ((u(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity v_inst '
  do k = 1, o
  write(61,*) ((v(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity w_inst '
  do k = 1, o
  write(61,*) ((w(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'pressure p_fluid'
  do k = 1, o
  write(61,*) ((p_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'pressure P_all'
  do k = 1, o
  write(61,*) ((p(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'porosity'
  do k = 1, o
  write(61,*) ((porosity(i,j,k), i=1,m),j=1,n)
  end do
  
  !write(*,*) "check", (porosity(i,10), i=1,m)
  
  close (61)
  ! ----------------
  
  ! ----------------
  ! surface profile
  open (62, file='surface_profile.dat', status='replace')
  
  do k=1,o
  do j=1,n
  do i=1,m
  
  if( porosity(i,j,k) < pmax .and. porosity(i,j,k) > pmin )then
   write(62,*) xp(i), yp(j), zp(i), p_cnt(i,j,k), porosity(i,j,k)
  end if
  
  end do
  end do
  end do
  close (62)
  ! ----------------
  
  
  return
  end subroutine output_solution_post
  !******************
  
  !******************
  subroutine  output_paraview (p, u, v, w, porosity, xp, yp, zp, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   real,intent(in),dimension(0:md, 0:nd, 0:od)::u, v, w, p
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   integer,intent(in)::m, n, o
   integer::i, j, k
  
  ! local variables
   real,dimension(0:md,0:nd,0:od):: div
  
   character(len=50)::csv_file
   character(len=50)::output_folder
  
   namelist /directory_control/csv_file, output_folder
   open(11,file="controlDict.txt",status="old",action="read")
   read(11,nml=directory_control)
   close(11)
  
  open(50,file=trim(output_folder)//'/output_paraview.vtk',status="unknown",form="formatted",position="rewind")
  !open(*,file='solution.vtk',status="replace")
  ! ----------------
  
      write(50,"('# vtk DataFile Version 3.0')")
      write(50,"('3D flow')")
      write(50,"('ASCII ')")
      
      write(50,"('DATASET STRUCTURED_GRID')")
      write(50,"('DIMENSIONS ',3(1x,i4))") m, n, o
      
      write(50,"('POINTS ',i9,' float')") m*n*o
      do k=1,o
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") xp(i), yp(j), zp(k)
      enddo
      enddo
      enddo
      
      write(50,"('POINT_DATA ',i9)") m*n*o
      
  !! velocity vector
      write(50,"('VECTORS velocity float')")
      do k=1,o
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j,k), v(i,j,k), w(i,j,k)
      enddo
      enddo
      enddo
  
  !! velocity vector
      write(50,"('VECTORS velocityInFluid float')")
      do k=1,o
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j,k)*porosity(i,j,k), v(i,j,k)*porosity(i,j,k), w(i,j,k)*porosity(i,j,k)
      enddo
      enddo
      enddo
        
  !! pressure
      write(50,"('SCALARS pressure float')")
      write(50,"('LOOKUP_TABLE default')")
      do k=1,o
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") p(i,j,k)
      enddo
      enddo
      enddo
  
  do i=1,m
  do j=1,n
  do k=1,o
   div(i,j,k)= (u(i+1,j,k)-u(i-1,j,k))/(xp(i+1)-xp(i-1)) &
              +(v(i,j+1,k)-v(i,j-1,k))/(yp(j+1)-yp(j-1)) &
              +(v(i,j,k+1)-v(i,j,k-1))/(zp(k+1)-zp(k-1))
  end do
  end do
  end do
  
  !! divergent velocity
      write(50,"('SCALARS VelocityDivergent float')")
      write(50,"('LOOKUP_TABLE default')")
      do k=1,o
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") div(i,j,k)
      enddo
      enddo
      enddo
    
  !! porosity
      write(50,"('SCALARS porosity float')")
      write(50,"('LOOKUP_TABLE default')")
      do k=1,o
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") porosity(i,j,k)
      enddo
      enddo
      enddo
  
  ! ----------------
  close(50)
  
  return
  end subroutine  output_paraview
  !******************
  
  !******************
  subroutine  output_divergent (p, u, v, w, porosity, dx, dy, dz, m, n, o)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd,0:od)::u, v, w, p 
   real,intent(in),dimension(0:md,0:nd,0:od)::porosity
   real,intent(in)::dx, dy, dz
   integer,intent(in)::m, n, o
  
  ! local variables
  integer::i, j, k
  real,dimension(0:md,0:nd,0:od)::div
  
  open (62, file='divergent.dat', status='replace')
  ! ----------------
  
  do i = 1, m
  do j = 1, n
  do k = 1, o
  div(i,j,k)= ((porosity(i+1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i+1,j,k))/2      &
              -(porosity(i-1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i-1,j,k))/2 )/dx &
             +((porosity(i,j+1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j+1,k))/2      &
              -(porosity(i,j-1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j-1,k))/2 )/dy &
             +((porosity(i,j,k+1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k+1))/2      &
              -(porosity(i,j,k-1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k-1))/2 )/dz
  end do
  end do
  end do
  
  write(62,*)
  write(62,*)'porosity'
  do k = 1, o
  do j = 1, n
  write(62,*) (porosity(i,j,k), i=1,m)
  end do
  end do
  
  write(62,*)
  write(62,*)'divergent velocity'
  do k = 1, o
  do j = 1, n
  write(62,*) (div(i,j,k), i=1,m)
  end do
  end do
  write(62,*)
  
  ! ----------------
  close (62)
  
  end subroutine  output_divergent
  !******************
  
  !******************
  subroutine  output_paraview_temp (p, u, v, w, porosity, xp, yp, zp, m, n, o, istep)
   implicit none
   integer,parameter::md=200, nd = 200, od = 200     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:od)::zp
   real,intent(in),dimension(0:md, 0:nd, 0:od)::u, v, w, p
   real,intent(in),dimension(0:md, 0:nd, 0:od)::porosity
   integer,intent(in)::m, n, o, istep
  
  ! -- local variable
   real,dimension(0:md,0:nd,0:od):: div
   integer::i, j, k
   character(5)::number
   character(len=50)::csv_file
   character(len=50)::output_folder
  ! -- open file
  
   namelist /directory_control/csv_file, output_folder
   open(11,file="controlDict.txt",status="old",action="read")
   read(11,nml=directory_control)
   close(11)
  
  write(number,"(I5.5)")istep
  
  open(65,file=trim(output_folder)//"/output_"//number//".vtk",status="unknown",form="formatted",position="rewind")
  !open(*,file='solution.vtk',status="replace")
  ! ----------------
  
  write(65,"('# vtk DataFile Version 3.0')")
  write(65,"('3D flow')")
  write(65,"('ASCII ')")
  
  write(65,"('DATASET STRUCTURED_GRID')")
  write(65,"('DIMENSIONS ',3(1x,i4))") m, n, o
  
  write(65,"('POINTS ',i9,' float')") m*n*o
  do k=1,o
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") xp(i), yp(j), zp(k)
  enddo
  enddo
  enddo
  
  write(65,"('POINT_DATA ',i9)") m*n*o
  
!! velocity vector
  write(65,"('VECTORS velocity float')")
  do k=1,o
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") u(i,j,k), v(i,j,k), w(i,j,k)
  enddo
  enddo
  enddo

!! velocity vector
  write(65,"('VECTORS velocityInFluid float')")
  do k=1,o
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") u(i,j,k)*porosity(i,j,k), v(i,j,k)*porosity(i,j,k), w(i,j,k)*porosity(i,j,k)
  enddo
  enddo
  enddo
    
!! pressure
  write(65,"('SCALARS pressure float')")
  write(65,"('LOOKUP_TABLE default')")
  do k=1,o
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") p(i,j,k)
  enddo
  enddo
  enddo

  do i=1,m
  do j=1,n
  do k=1,o
  div(i,j,k)= (u(i+1,j,k)-u(i-1,j,k))/(xp(i+1)-xp(i-1)) &
            +(v(i,j+1,k)-v(i,j-1,k))/(yp(j+1)-yp(j-1)) &
            +(v(i,j,k+1)-v(i,j,k-1))/(zp(k+1)-zp(k-1))
  end do
  end do
  end do

!! divergent velocity
  write(65,"('SCALARS VelocityDivergent float')")
  write(65,"('LOOKUP_TABLE default')")
  do k=1,o
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") div(i,j,k)
  enddo
  enddo
  enddo

!! porosity
  write(65,"('SCALARS porosity float')")
  write(65,"('LOOKUP_TABLE default')")
  do k=1,o
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") porosity(i,j,k)
  enddo
  enddo
  enddo

  ! ----------------
  close(65)
  
  return
  end subroutine  output_paraview_temp
  !******************

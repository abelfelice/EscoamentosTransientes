program acoplamentoPV_Abel
  implicit none  

  ! 1. Declaração dos Parametros de entrada e do problema

  integer :: nx, ny, i, j, cont, k, count_start, count_end, count_rate, snap_freq, snap_id
  double precision :: H, L, dx, dy, dt, rho, nu, U_tampa, U_paredeEs, U_paredeDir, U_paredebx, &
  V_tampa, V_paredeEs, V_paredeDir, V_paredebx, CFL, ae, ad, ac, ab, ap, tol, erro, Re_dx, Pe_dx, &
  dt_advc, dt_dfu, t_total, t_simulacao, Reynolds, tempo_execucao
  double precision, dimension (:,:), allocatable :: u,v,p, u_estrela, v_estrela, P_linha, B, P_linha_erro
  character(len=30) :: snap_filename


  ! 2. Implemementação dos valores nas variáveis declaradas e condições de contorno e valor inicial
  call system_clock(count_start, count_rate)
  H = 1.0d0
  L = 1.0d0
  U_tampa = 1.0d0
  U_paredeEs = 0.0d0 
  U_paredeDir = 0.0d0
  U_paredebx = 0.0d0
  V_tampa = 0.0d0
  V_paredeEs = 0.0d0
  V_paredeDir = 0.0d0
  V_paredebx = 0.0d0
  rho = 1.0d0
  nu = 1.0d-2/10.0d0 ! Controle o reynolds pelo denominador
  nx = 258 ! Usar um numero par
  ny = 258 ! Usar um numero par
  CFL = 0.1d0
  dx = L/dble(nx)
  dy = H/dble(ny)
  dt_dfu = (1.0d0 / (2.0d0 * nu)) * ((dx*dy)**2) / (dx**2 + dy**2)
  dt_advc = cfl*(min(dx,dy)/max(U_tampa, U_paredeEs, U_paredeDir, U_paredebx, &
  V_tampa, V_paredeEs, V_paredeDir, V_paredebx))
  dt = min(dt_advc,dt_dfu)

  Re_dx = 0.0d0
  Pe_dx = 0.0d0
  t_total = 120.0d0
  t_simulacao = 0.0d0
  Reynolds = (max(U_tampa, U_paredeEs, U_paredeDir, U_paredebx, &
  V_tampa, V_paredeEs, V_paredeDir, V_paredebx)*max(L,h))/nu


  Re_dx = (max(U_tampa, U_paredeEs, U_paredeDir, U_paredebx, &
  V_tampa, V_paredeEs, V_paredeDir, V_paredebx)*dx)/nu


  snap_freq = 200       ! salvar a cada 200 passos do loop temporal
  snap_id = 0           ! contador de snapshots (começa em zero)

 
  call execute_command_line('if not exist snapshots mkdir snapshots', wait=.true.)

  write(*,*) "O Numero de Reynolds do escoamento eh", Reynolds
  write(*,*) "O Reynolds de celula eh", Re_dx
  write(*,*) "O passo de tempo eh ", dt,"s"


  ae = 0.0d0
  ad = 0.0d0
  ab = 0.0d0
  ac = 0.0d0
  ap = 0.0d0
  tol = 1.0d-6
  erro = 1.0d0
  cont = 0 
  k = 0 

  allocate (u(1:ny, 1:nx+1))
  allocate (v(1:ny+1, 1:nx))
  allocate (p(1:ny, 1:nx))
  allocate (u_estrela(1:ny, 1:nx+1))
  allocate (v_estrela(1:ny+1, 1:nx))
  allocate (P_linha(1:ny, 1:nx))
  allocate (P_linha_erro(1:ny, 1:nx))
  allocate (B(1:ny,1:nx))


  u = 0.0d0
  v = 0.0d0
  p = 0.0d0

  P_linha = 0.0d0
  P_linha_erro = 0.0d0

  u(:,1) = U_paredeEs 
  u(:,nx+1) = U_paredeDir
  v(1,:) = V_paredebx
  v(ny+1,:) = v_tampa

  u_estrela = u
  v_estrela = v


  do while (t_simulacao<=t_total .and. k<=1.0d8) 

    k = k + 1
    t_simulacao = t_simulacao + dt
    write (*,*) "tempo de simulacao", t_simulacao

   ! Equação de u estrela para os meios 
    do i = 2,nx
      do j = 2,ny-1

          u_estrela(j,i) = u(j,i) - dt*((+u(j,i)*((u(j,i+1)-u(j,i-1))/(2.0d0*dx)))  + &
          ((v(j,i)+v(j,i-1)+v(j+1,i)+v(j+1,i-1))/4.0d0)*((u(j+1,i)-u(j-1,i))/(2.0d0*dy)) + &
          (+1.0d0/rho)*((P(j,i)-P(j,i-1))/dx) + & 
          -nu*((u(j,i-1) - 2.0d0*u(j,i) + u(j,i+1))/(dx**2.0d0) + &
          (u(j+1,i) - 2.0d0*u(j,i) + u(j-1,i))/(dy**2.0d0)))

      end do
    end do


  ! Equação de u estrela para a parede inferior
    j = 1
    do i = 2,nx
      
          u_estrela(j,i) = u(j,i) - dt*((+u(j,i)*((u(j,i+1)-u(j,i-1))/(2.0d0*dx)))  + &
          ((2.0d0*V_paredebx + v(j+1,i) + v(j+1,i-1))/4.0d0)*((u(j+1,i)-2.0d0*U_paredebx+u(j,i))/(2.0d0*dy)) + &
          (+1.0d0/rho)*((P(j,i)-P(j,i-1))/dx) + & 
          -nu*((u(j,i-1) - 2.0d0*u(j,i) + u(j,i+1))/(dx**2.0d0) + &
          (u(j+1,i) - 3.0d0*u(j,i) + 2.0d0*U_paredebx)/(dy**2.0d0)))
    end do 

  ! Equação de u estrela para a parede superior (tampa)
    j = ny
    do i = 2,nx
      
          u_estrela(j,i) = u(j,i) - dt*((+u(j,i)*((u(j,i+1)-u(j,i-1))/(2.0d0*dx)) + &
          ((2.0d0*V_tampa + v(j,i-1) + v(j,i))/4.0d0)*((2.0d0*U_tampa-u(ny,i)-u(ny-1,i))/(2.0d0*dy))) + &
          (+1.0d0/rho)*((P(j,i)-P(j,i-1))/dx) + & 
          -nu*((u(j,i-1) - 2.0d0*u(j,i) + u(j,i+1))/(dx**2.0d0) + &
          (2.0d0*U_tampa-3.0d0*(u(ny,i))+u(ny-1,i))/(dy**2.0d0)))

    end do 

  ! Equação de v estrela para os meios 
    do i = 2,nx-1
      do j = 2,ny 

          v_estrela(j,i) = v(j,i) -dt*(((u(j,i)+u(j,i+1)+u(j-1,i)+u(j-1,i+1))/4.0d0)*&
          ((v(j,i+1)-v(j,i-1))/(2.0d0*dx)) + &
          v(j,i)*((v(j+1,i)-v(j-1,i))/(2.0d0*dy)) + &
          (1.0d0/rho)*((P(j,i)-P(j-1,i))/dy) + & 
          -nu*((v(j,i-1) - 2.0d0*v(j,i) + v(j,i+1))/(dx**2.0d0) + &
          (v(j-1,i) - 2.0d0*v(j,i) + v(j+1,i))/(dy**2.0d0)))

      end do
    end do


! V_estrela par LATERAL ESQUERDA 
    
    i = 1
    do j = 2, ny
          
          v_estrela(j,i) = v(j,i) - dt * ( &
          ((2.0d0*U_paredeEs + u(j,i+1) + u(j-1,i+1))/4.0d0) * &
          ((v(j,i+1) - (2.0d0*V_paredeEs - v(j,i)))/(2.0d0*dx)) + &
          v(j,i) * ((v(j+1,i)-v(j-1,i))/(2.0d0*dy)) + &
          (1.0d0/rho)*((P(j,i)-P(j-1,i))/dy) + & 
          -nu * ( &
             ((2.0d0*V_paredeEs - 3.0d0*v(j,i) + v(j,i+1))/(dx**2.0d0)) + &
             ((v(j-1,i) - 2.0d0*v(j,i) + v(j+1,i))/(dy**2.0d0)) &
          ))
    end do

! V_estrela: LATERAL DIREITA 

    i = nx
    do j = 2, ny

          v_estrela(j,i) = v(j,i) - dt * ( &
          ((2.0d0*U_paredeDir + u(j,i) + u(j-1,i))/4.0d0) * &
          (((2.0d0*V_paredeDir - v(j,i)) - v(j,i-1))/(2.0d0*dx)) + &
          v(j,i) * ((v(j+1,i)-v(j-1,i))/(2.0d0*dy)) + &
          (1.0d0/rho)*((P(j,i)-P(j-1,i))/dy) + & 
          -nu * ( &
             ((v(j,i-1) - 3.0d0*v(j,i) + 2.0d0*V_paredeDir)/(dx**2.0d0)) + &
             ((v(j-1,i) - 2.0d0*v(j,i) + v(j+1,i))/(dy**2.0d0)) &
          ))
    end do


! Equação de Poison
    
    ! Termo B

    do i = 1,nx
      do j = 1,ny

        B(j,i) = (rho/dt) * ( (u_estrela(j,i+1) - u_estrela(j,i))/dx + &
        (v_estrela(j+1,i) - v_estrela(j,i))/dy )

      end do
    end do
    
    erro = 1.0d0
    cont = 0
    P_linha = 0.0d0
    ! Loop do método iterativo Gauss-Seidel
    do while (erro >= tol)

      cont = cont + 1

      P_linha_erro = P_linha

    ! Quina inferior Esquerda

      i = 1
      j = 1
      ae = 0.0d0
      ad = 1.0d0/dx**2.0d0
      ab = 0.0d0
      ac = 1.0d0/dy**2.0d0
      ap = (-1.0d0/dx**2.0d0)+(-1.0d0/dy**2.0d0)

      P_linha(j,i) = (B(j,i) - (ad * P_linha(j, i+1) + ac * P_linha(j+1, i))) / ap

    ! Quina inferior Direita

      i = nx
      j = 1
      ae = 1.0d0/dx**2.0d0
      ad = 0.0d0
      ab = 0.0d0
      ac = 1.0d0/dy**2.0d0
      ap = (-1.0d0/dx**2.0d0)+(-1.0d0/dy**2.0d0)

      P_linha(j,i) = (B(j,i) - (ae * P_linha(j, i-1) + ac * P_linha(j+1, i))) / ap


    ! Quina Superior Esquerda

      i = 1
      j = ny
      ae = 0.0d0
      ad = 1.0d0/dx**2.0d0
      ac = 0.0d0
      ab = 1.0d0/dy**2.0d0
      ap = (-1.0d0/dx**2.0d0)+(-1.0d0/dy**2.0d0)

      P_linha(j,i) = (B(j,i) - (+ad*P_linha(j,i+1)+ab*P_linha(j-1,i))) / ap


    ! Quina Superior Direita

      i = nx
      j = ny
      ae = 1.0d0/dx**2.0d0
      ad = 0.0d0
      ac = 0.0d0
      ab = 1.0d0/dy**2.0d0
      ap = (-1.0d0/dx**2.0d0)+(-1.0d0/dy**2.0d0)

      P_linha(j,i) = (B(j,i) - (ae * P_linha(j, i-1) + ab * P_linha(j-1, i))) / ap


    ! Parede Inferior 

      ae = 1.0d0/dx**2.0d0
      ad = 1.0d0/dx**2.0d0
      ac = 1.0d0/dy**2.0d0
      ab = 0.0d0
      ap = (-2.0d0/dx**2.0d0)+(-1.0d0/dy**2.0d0)
      j = 1
      do i = 2,(nx-1)

        P_linha(j,i) = (B(j,i) - (ae*P_linha(j,i-1) + ad*P_linha(j,i+1) + ac*P_linha(j+1,i)))/ap

      end do


    ! Parede Superior

      ae = 1.0d0/dx**2.0d0
      ad = 1.0d0/dx**2.0d0
      ac = 0.0d0
      ab = 1.0d0/dy**2.0d0
      ap = (-2.0d0/dx**2.0d0)+(-1.0d0/dy**2.0d0)
      j = ny
      do i = 2,(nx-1)

        P_linha(j,i) = (B(j,i) - (ae*P_linha(j,i-1) + ad*P_linha(j,i+1) + ab*P_linha(j-1,i)))/ap

      end do

    ! Parede Lateral Esquerda

      ae = 0.0d0
      ad = 1.0d0/dx**2.0d0
      ac = 1.0d0/dy**2.0d0
      ab = 1.0d0/dy**2.0d0
      ap = (-1.0d0/dx**2.0d0)+(-2.0d0/dy**2.0d0)
      i = 1 
      do j = 2,(ny-1)

        P_linha(j,i) = (B(j,i) - (ad*P_linha(j,i+1) + ac*P_linha(j+1,i) + ab*P_linha(j-1,i)))/ap

      end do


    ! Parede Lateral Direita

      ae = 1.0d0/dx**2.0d0
      ad = 0.0d0
      ac = 1.0d0/dy**2.0d0
      ab = 1.0d0/dy**2.0d0
      ap = (-1.0d0/dx**2.0d0)+(-2.0d0/dy**2.0d0)
      i = nx 
      do j = 2,(ny-1)

        P_linha(j,i) = (B(j,i) - (ae*P_linha(j,i-1) + ac*P_linha(j+1,i) + ab*P_linha(j-1,i)))/ap

      end do


    ! Poison Central 

      ae = 1.0d0/dx**2.0d0
      ad = 1.0d0/dx**2.0d0
      ac = 1.0d0/dy**2.0d0
      ab = 1.0d0/dy**2.0d0
      ap = (-2.0d0/dx**2.0d0)+(-2.0d0/dy**2.0d0)
      do i = 2,nx-1
        do j = 2,ny-1

          P_linha(j,i) = (B(j,i) - (ae*P_linha(j,i-1) + ad*P_linha(j,i+1) + &
                          ac*P_linha(j+1,i) + ab*P_linha(j-1,i)))/ap

        end do 
      end do 

      P_linha(2, 2) = 0.0d0
  
      erro = maxval(dabs(P_linha - P_linha_erro))

      
    end do

    

    P = P_linha + P


! Calculo de u verdadeiro

    do i = 2,nx
      do j = 1,ny

        u(j,i) = u_estrela(j,i) - (dt/rho) * ((P_linha(j,i) - P_linha(j,i-1)) / dx)
        

      end do 
    end do

! Calculo de v verdadeiro

    do i = 1,nx
      do j = 2,ny

        v(j,i) = v_estrela(j,i) - (dt/rho) * ((P_linha(j,i) - P_linha(j-1,i)) / dy)
        
      end do 
    end do

    write(*,*) "it sys:", cont, "it time:", k
    write(*,*) "p_linha:", maxval(dabs(P_linha))
    write(*,*) "max_v:", maxval(dabs(v)), maxval(abs((v-v_estrela)))
    write(*,*) "max_u:", maxval(dabs(u)), maxval(abs((u-u_estrela)))


      !Salvar snapshot a cada snap_freq passos
    if (mod(k, snap_freq) == 0) then
      snap_id = snap_id + 1
      write(snap_filename, '("snapshots/snap_", i5.5, ".dat")') snap_id
      open(unit=20, file=snap_filename)
      write(20, *) "# tempo = ", t_simulacao
      do i = 1, nx
        do j = 1, ny
         write(20, *) (i-0.5d0)*dx, (j-0.5d0)*dy,    &
                       (u(j,i) + u(j,i+1)) / 2.0d0,    &
                       (v(j,i) + v(j+1,i)) / 2.0d0,    &
                       P(j,i)
        end do
      end do
      close(20)
    end if


  end do 

write(*,*) "Simulacao concluida"


!Impressão dos reslutados 

open(unit = 2 , file = "resultadosacpv.dat")
write (2,*) "Tempo total de simulacao" , t_total
write (2,*) "Reynolds do escoamento" , Reynolds
write (2,*) "        x                              y                                u                             &
  v                    P " 


do i = 1,nx
  do j = 1,ny

    write(2,*) ((i-0.5d0)*dx) , ((j-0.5d0)*dy), (u(j,i+1)+u(j,i))/2.0d0 , (v(j+1,i)+v(j,i))/2.0d0, P(j,i)

  end do 
end do

close(unit = 2)

open(unit = 3 , file = "resultadosmeiosU.dat")

do j = 1, ny

  write(3,*) u(j,(Nx/2))

end do

close(unit = 3)

open(unit = 4 , file = "resultadosmeiosV.dat")

do i = 1, nx

  write(4,*) v(Ny/2,i)

end do

close(unit = 4)

write(*,*) "Dados exportados para resultadosacpv.dat"

call system_clock(count_end)
tempo_execucao = real(count_end - count_start, 8) / real(count_rate, 8)
write(*,*) "tempo de execucao em segundos", tempo_execucao 

end program acoplamentoPV_Abel

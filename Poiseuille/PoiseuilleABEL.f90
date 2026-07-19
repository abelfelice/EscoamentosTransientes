program Poiseuille2D
    implicit none


    ! 1. Parâmetros de entrada e do problema

    ! Parâmetros da malha e domínio
    integer, parameter :: N = 100          ! Número de nós verticais
    integer, parameter :: M = N           ! Número de nós horizontais
    real(kind=8), parameter :: h = 0.01   ! Distância vertical [m]
    real(kind=8), parameter :: L = h      ! Distância horizontal [m]

    ! Parâmetros da simulação
    real(kind=8), parameter :: tf = 60   ! Tempo total de simulação [s]
    real(kind=8), parameter :: CFL = 1.0 ! Não precisa de CFL para esse método

    ! Propriedades do fluido e escoamento
    real(kind=8), parameter :: nu = 1.0e-6   ! Viscosidade cintempomática [m^2/s]
    real(kind=8), parameter :: rho = 1000    ! Massa específica [kg/m^3]
    real(kind=8), parameter :: gradp = -10  ! Gradiente de Pressão [Pa/m]
    
    ! Parâmetros do Gauss-Seidel
    integer, parameter :: max_iter = 1000
    real(kind=8), parameter :: tol = 1.0e-6



    ! 2. Declaração de Variáveis

    ! Variáveis calculadas
    real(kind=8) :: delta, dt
    integer      :: Nt
    real(kind=8) :: y(N), z(M)
    
    ! Coeficientes da equação discretizada
    real(kind=8) :: c, a0, ae, ad, ab, ac
    
    ! Matrizes de velocidade
    real(kind=8), allocatable :: u(:,:,:)
    real(kind=8) :: u_old(N, M) ! Array para checar convergência
    
    ! Variáveis de controle e temporárias
    integer :: ntempo, i, j, iter
    real(kind=8) :: erro, b

    ! Variáveis do temporizador
    integer :: count_start, count_end, count_rate
    real(kind=8) :: tempo_execucao

    ! 3. Cálculo de Parâmetros Derivados e Coeficientes

    delta = h / real(N - 1, 8)
    dt = CFL * (delta**2) / nu
    Nt = ceiling(tf / dt)

    ! Vetores de coordenadas
    do i = 1, N
        y(i) = real(i-1, 8) * delta
        z(i) = real(i-1, 8) * delta
    end do
    
    ! Coeficientes
    c = (1.0_8 / rho) * gradp
    a0 = 1.0_8 + 4.0_8 * (nu * dt / delta**2)
    ae = -nu * dt / delta**2
    ad = -nu * dt / delta**2
    ab = -nu * dt / delta**2
    ac = -nu * dt / delta**2
    
    ! Impressão dos parâmetros para verificação
    print *, "Simulacao 2D com Fortran 90 (Gauss-Seidel)"
    print '("Dimensao da malha: ", I0, " x ", I0)', N, M
    print '("Passo de tempo (dt): ", E12.5, " s")', dt
    print '("Numero de passos de tempo (Nt): ", I0)', Nt


    ! 4. Alocação e Inicialização

    print *, "Alocando memoria para o campo de velocidades..."
    allocate(u(N, M, Nt))
    
    ! Inicializa toda a matriz com zero
    u = 0.0
    print *, "Memoria alocada e condicoes inicializadas."
    

    
    ! 5. Loop Principal no Tempo
    print *, "Iniciando loop temporal..."
    call system_clock(count_start, count_rate)
    
    do ntempo = 2, Nt
        u(:,:,ntempo) = u(:,:,ntempo-1)
        
        do iter = 1, max_iter
            u_old = u(:,:,ntempo) 
            
            do j = 2, M - 1 
                do i = 2, N - 1
                    b = -c * dt + u(i, j, ntempo-1)
                    u(i, j, ntempo) = (b - (ac * u(i-1, j, ntempo) + ab * u(i+1, j, ntempo) + &
                                      ad * u(i, j-1, ntempo) + ae * u(i, j+1, ntempo))) / a0
                end do
            end do
            erro = maxval(abs(u(:,:,ntempo) - u_old)) !Isso pega o maior valor em módulo, que é no meio (crista)
            if (erro < tol) then
                exit
            end if
        end do
        !write(*,*) iter
    end do
    
    call system_clock(count_end)
    tempo_execucao = real(count_end - count_start, 8) / real(count_rate, 8)

    print *, "Loop temporal finalizado."
    print '("Tempo de execucao do loop: ", F10.3, " segundos")', tempo_execucao
    print *, "--------------------------------------------"
    
    ! 6. Pós-processamento: Salvar Resultados Finais
    print *, "Salvando o campo de velocidades final em 'poiseuille_fortran_final.dat'..."
    open(unit=10, file='poiseuille_fortran_final.dat', status='replace', action='write')
    
    write(10, '(A, A, A, A, A)') '# Coordenada Y (m)', 'Coordenada Z (m)', 'Velocidade U (m/s)'
    
    do j = 1, M
        do i = 1, N
            write(10, '(E15.7, 2x, E15.7, 2x, E15.7)') y(i), z(j), u(i, j, Nt)
        end do
        write(10, *) 
    end do
    
    close(10)
    print *, "Resultados salvos."

    ! 7. Finalização
    deallocate(u)
    print *, "Simulacao concluida com sucesso."
    
end program Poiseuille2D

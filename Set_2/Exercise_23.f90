program Askisi_23

    implicit none

    integer, parameter :: N = 10., Ne = N - 1, L = 5., omega = 2.
    real, parameter :: dt = 0.1, dx = float(L) / Ne

    real :: tm, Ml(2, 2), Ll(2, 2), A(N, N), b(N), Z(N, 2), T(N), Ta(N), ipivot(N), ierr, k, x(N)
    integer :: ie, i, j, il, jl

    k = sqrt(omega / 2.)

    ! Fill connective matrix
    do ie = 1, Ne
        Z(ie, 1) = ie
        Z(ie, 2) = ie + 1
    end do

    do i = 1, N
        x(i) = (i - 1) * dx
    end do

    ! Initial condition
    T = 0.

    ! Fill local mass matrix
    Ll(1, 1) = 1.
    Ll(1, 2) = -1.
    Ll(2, 1) = -1.
    Ll(2, 2) = 1.

    Ll = Ll / dx

    ! Fill local Laplace matrix
    Ml(1, 1) = 2.
    Ml(1, 2) = 1.
    Ml(2, 1) = 1.
    Ml(2, 2) = 2.

    Ml = Ml * dx / 6.

    ! Time step
    tm = 0.
    10 tm = tm + dt

    do i = 1, N
        Ta(i) = exp(-k * x(i)) * sin(omega * tm - k * x(i))
    end do

    A = 0.
    b = 0.

    ! Fill matrix A
    do ie = 1, Ne
        do i = 1, 2
            il = Z(ie, i)
            do j = 1, 2
                jl = Z(ie, j)
                A(il, jl) = A(il, jl) + 2 * Ml(i, j) + Ll(i, j) * dt
                b(il) = b(il) + (2 * Ml(i, j) - Ll(i, j) * dt) * T(jl)
            end do
        end do
    end do

    ! Dirichlet boundary conditions
    A(1, 1:N) = 0
    A(1, 1) = 1
    b(1) = sin(omega * tm)

    A(N, 1:N) = 0
    A(N, N) = 1
    b(N) = 0

    call sgesv(N, 1, A, N, ipivot, b, N, ierr)
    T = b

    if (ierr .ne. 0) then
        write(*, *) "Error while solving the system"
        stop
    end if

    write(*, *) "t = ", tm
    write(*, *) "x: ", x
    write(*, *) "T: ", T
    write(*, *) "Ta: ", Ta
    write(*, *) " "

    if (tm < 15) then
        goto 10
    end if

end program

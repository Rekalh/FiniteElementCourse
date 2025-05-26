program Askisi_21

    implicit none

    integer, parameter :: N = 20, Ne = N - 1

    integer :: i, j, pivot(N), ierr, il, jl, ie
    real :: D(N, N), Tfe(N), Tfd(N), b(N), beta(N), dx, lambda, x, Z(Ne, 2), Ml(2, 2), Ll(2, 2)
    data Ll/1, -1, -1, 1/, Ml/2, 1, 1, 2/

    lambda = 1
    dx = 1. / (N - 1)

    Ll = Ll / dx
    Ml = Ml * dx / 6

    D = 0.
    Tfe = 0.
    Tfd = 0.
    b = 0.
    beta = 0.

    ! Finite elements
    do ie = 1, Ne
        Z(ie, 1) = ie
        Z(ie, 2) = ie + 1
    end do

    do i = 2, (N - 1)
        beta(i) = exp(-lambda * (i - 1) * dx)
    end do

    do ie = 1, Ne
        do i = 1, 2
            il = Z(ie, i)

            do j = 1, 2
                jl = Z(ie, j)

                D(il, jl) = D(il, jl) + Ll(i, j)
                b(il) = b(il) + Ml(i, j) * beta(jl)
            end do
        end do
    end do

    ! Boundary conditions

    ! Dirichlet (left)
    D(1, 1:N) = 0.
    D(1, 1) = 1.
    b(1) = 0.

    ! Neumann (right)
    b(N) = b(N) + 0

    call sgesv(N, 1, D, N, pivot, b, N, ierr)
    Tfe = b

    if (ierr .ne. 0) then
        write(*, *) "Error while solving the system"
        stop
    end if

    ! Finite differences
    pivot = 0.
    ierr = 0

    D = 0.
    b = 0.

    do i = 2, (N - 1)
        D(i, i) = -2.
        D(i, i - 1) = 1.
        D(i, i + 1) = 1.

        b(i) = -dx ** 2 * exp(-lambda * (i - 1) * dx)
    end do

    ! Boundary conditions

    ! Dirichlet (left)
    D(1, 1) = 1.
    b(1) = 0.

    ! Neumann (right)
    D(N, N) = 1.
    D(N, N - 1) = -1.
    b(N) = 0.

    call sgesv(N, 1, D, N, pivot, b, N, ierr)
    Tfd = b

    if (ierr .ne. 0) then
        write(*, *) "Error while solving the system"
        stop
    end if

    do i = 1, N
        x = (i - 1) * dx
        write(*, *) x, Tfe(i), Tfd(i), -(lambda * exp(-lambda) * x + exp(-lambda * x) - 1) / (lambda ** 2)
    end do

end program

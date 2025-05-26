! The file contains 3 separate programs

program Askisi_22_ImplicitEuler

    implicit none

    integer, parameter :: N = 100
    real, parameter :: dt = 0.1, lambda = 2

    real :: A(2, 2), f(2, N)
    integer :: i

    A(1, 1) = 1
    A(1, 2) = dt
    A(2, 1) = -(lambda ** 2) * dt
    A(2, 2) = 1

    f(1, 1) = 0
    f(2, 1) = lambda

    do i = 2, N
        f(:, i) = matmul(A, f(:, i - 1))

        write(*, *) f(1, i-1)
    end do

end program

program Askisi_22_ExplicitEuler

    implicit none

    integer, parameter :: N = 100
    real, parameter :: dt = 0.1, lambda = 2

    real :: A(2, 2), f(2, N)
    integer :: i

    A(1, 1) = 1
    A(1, 2) = dt
    A(2, 1) = -(lambda ** 2) * dt
    A(2, 2) = 1

    A = A / (1 + (lambda ** 2) * (dt ** 2))

    f(1, 1) = 0
    f(2, 1) = lambda

    do i = 2, N
        f(:, i) = matmul(A, f(:, i - 1))

        write(*, *) f(1, i-1)
    end do

end program

program Askisi_22_CrankNicholson

    implicit none

    integer, parameter :: N = 100
    real, parameter :: dt = 0.1, lambda = 2

    real :: A(2, 2), f(2, N)
    integer :: i

    A(1, 1) = 2
    A(1, 2) = dt
    A(2, 1) = -(lambda ** 2) * dt
    A(2, 2) = 2

    A = matmul(A, A)

    A = A / (4 + (lambda ** 2) * (dt ** 2))

    f(1, 1) = 0
    f(2, 1) = lambda

    do i = 2, N
        f(:, i) = matmul(A, f(:, i - 1))

        write(*, *) f(1, i-1)
    end do

end program


program hello_openmp
    ! compile and run with: gfortran -fopenmp hello_openmp.f90 -o hello_openmp.exe && ./hello_openmp.exe
    use omp_lib

    implicit none

    integer, parameter :: num_iterations = 10
    integer :: i

    !$omp parallel do
    do i = 1, num_iterations
        print *, "Hello on iteration ", i, " from thread ", omp_get_thread_num()
    end do
end program
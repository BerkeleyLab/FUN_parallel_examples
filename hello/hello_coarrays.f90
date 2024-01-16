program hello_coarrays
    ! compile and run with:
    !   $ ftn hello_coarrays.f90 -o hello_coarrays.exe
    !   $ srun -N 1 -n 4 -t 00:00:10 --account=your_project --qos=debug --constraint=cpu ./hello_coarrays.exe
    implicit none

    print *, "Hello from image ", this_image(), " of ", num_images()
end program
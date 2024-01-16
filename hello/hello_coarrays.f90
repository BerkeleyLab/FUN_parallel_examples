program hello_coarrays
    ! compile and run with: caf hello_coarrays.f90 -o hello_coarrays.exe && cafrun -n 4 ./hello_coarrays.exe
    implicit none

    print *, "Hello from image ", this_image(), " of ", num_images()
end program
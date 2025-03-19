!PUT - puts command line arguments to stdout, for use in batch files

program put
    integer :: i
    character(len=256) :: arg

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        write(*, '(a)') trim(adjustl(arg))
    end do
end program

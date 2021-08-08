! Created by mus on 06/08/2021.

subroutine read_sim_parameters(A)
    implicit none
    character(len=40), intent(out)       :: A
    character(len=40)                    :: filename
    integer                              :: stat, i

    filename = "parameters.in"
    open(50,file=filename,iostat=stat)
    print*,stat
        read(50,*) A
        print*,A

    close(50)

end subroutine read_sim_parameters
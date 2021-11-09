! Created by mus on 30/07/2021.

subroutine readmodelfiles1D(v1D, rho1D, ne, modnameprefix)
    implicit none
    integer, intent(in)                           :: ne
    real (kind=8), dimension(ne), intent(out)     :: v1D, rho1D
    character (len=40)                            :: modname_vp, modname_rho, modnameprefix

    modname_vp     = TRIM(modnameprefix)//'_vs'
    modname_rho    = TRIM(modnameprefix)//'_rho'

    open(10,file=modname_vp,access='direct',recl=ne*8)
    read(10,rec=1) v1D
    close(10)

    open(11,file=modname_rho,access='direct',recl=ne*8)
    read(11,rec=1) rho1D
    close(11)



end subroutine readmodelfiles1D
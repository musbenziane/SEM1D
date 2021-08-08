


do i=1,ngll
    sum = 0.
    do j=1,ngll
        sum = sum + Kg(i,j) * u(j)
    end do
    KU (i)= sum
end do
F_KU = F - KU

do i=1,ngll
    sum = 0.
    do j=1,ngll
        sum = sum + Minv(i,j) * F_KU(j)
    end do
    Minv_x_F_KU(i) = sum
    !print*,i, sum
end do
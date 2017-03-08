Module proc
Use var
Implicit None

Contains
  ! 计算u子程序======================================================================
  Subroutine cal_u
  Implicit None
    ! 定义cal_u子程序用到的局部变量
    Real(8) :: CuP(0:nx+1, 0:ny+1), CuN(0:nx+1, 0:ny+1)
    Real(8) :: CvP(0:nx+1, 0:ny+1), CvN(0:nx+1, 0:ny+1)
    Real(8) :: du(0:nx+1, 0:ny+1)
    Real(8) :: Rx(0:nx+1, 0:ny+1)
    Real(8) :: hu, uu, vu
    Real(8) :: pgrdx, tx, corx
    Real(8) :: advx(0:nx+1, 0:ny+1)
    Real(8) :: div, div1, div2

    ! 先计算对流项advx
    Do j = 0, ny+1
      Do i = 0, nx
        CuP(i, j) = 0.25d0*(u_cur(i, j)+u_cur(i+1, j)+abs(u_cur(i, j))+abs(u_cur(i+1, j)))*dt/dx
        CuN(i, j) = 0.25d0*(u_cur(i, j)+u_cur(i+1, j)-abs(u_cur(i, j))-abs(u_cur(i+1, j)))*dt/dx
        CvP(i, j) = 0.25d0*(v_cur(i, j)+v_cur(i+1, j)+abs(v_cur(i, j))+abs(v_cur(i+1, j)))*dt/dy
        CvN(i, j) = 0.25d0*(v_cur(i, j)+v_cur(i+1, j)-abs(v_cur(i, j))-abs(v_cur(i+1, j)))*dt/dy
      End Do
    End Do

    B_cur(:, :) = u_cur(:, :)

    Call advect(CuP, CuN, CvP, CvN)

    Do j = 1, ny
      Do i = 1, nx
        div1 = 0.5d0*(u_cur(i+1, j)-u_cur(i-1, j))/dx
        div2 = 0.5d0*(v_cur(i, j)+v_cur(i+1, j)-v_cur(i, j-1)-v_cur(i+1, j-1))/dy
        div = dt*B_cur(i, j)*(div1+div2)  
        advx(i, j)= B_nex(i, j)+div
      End Do
    End Do

    ! 计算tx及pgrdx,与advx相加得到du
    Do j = 1, ny
      Do i = 1, nx
        ! 计算u的中间变量    
        pgrdx = -dt*g*(eta_cur(i+1, j)-eta_cur(i, j))/dx
        hu    = 0.5d0*(h(i, j)+h(i+1, j))
        uu    = u_cur(i, j)
        vu    = 0.25d0*(v_cur(i, j)+v_cur(i, j-1)+v_cur(i+1, j)+v_cur(i+1, j-1))
        If(hu>0.0) Then
          tx = dt*taux/(rho*hu)
          Rx(i, j) = dt*r*sqrt(uu**2+vu**2)/hu 
        End If
        du(i, j) = tx+pgrdx+advx(i, j)

        ! 计算u
        If (mask(i, j)==1) Then
          If ((mask(i+1, j)==1) .or. (du(i, j)>0.0)) u_nex(i, j) = (u_cur(i, j)+du(i, j))/(1.0d0+Rx(i, j))
        Else
          If ((mask(i+1, j)==1) .and. (du(i, j)<0.0)) u_nex(i, j) = (u_cur(i, j)+du(i, j))/(1.0d0+Rx(i, j))
        End If
      End Do
    End Do
  End Subroutine cal_u

  ! 计算v子程序======================================================================
  Subroutine cal_v
  Implicit None
    ! 定义cal_v子程序用到的局部变量
    Real(8) :: CuP(0:nx+1, 0:ny+1), CuN(0:nx+1, 0:ny+1)
    Real(8) :: CvP(0:nx+1, 0:ny+1), CvN(0:nx+1, 0:ny+1)
    Real(8) :: dv(0:nx+1, 0:ny+1)
    Real(8) :: Ry(0:nx+1, 0:ny+1)
    Real(8) :: hv, vv, uv
    Real(8) :: pgrdy, ty, cory
    Real(8) :: advy(0:nx+1, 0:ny+1)
    Real(8) :: div, div1, div2

    ! 先计算对流项advy
    Do j = 0, ny
      Do i = 0, nx+1
        CuP(i, j) = 0.25d0*(u_cur(i, j)+u_cur(i, j+1)+abs(u_cur(i, j))+abs(u_cur(i, j+1)))*dt/dx
        CuN(i, j) = 0.25d0*(u_cur(i, j)+u_cur(i, j+1)-abs(u_cur(i, j))-abs(u_cur(i, j+1)))*dt/dx
        CvP(i, j) = 0.25d0*(v_cur(i, j)+v_cur(i, j+1)+abs(v_cur(i, j))+abs(v_cur(i, j+1)))*dt/dy
        CvN(i, j) = 0.25d0*(v_cur(i, j)+v_cur(i, j+1)-abs(v_cur(i, j))-abs(v_cur(i, j+1)))*dt/dy
      End Do
    End Do

    B_cur(:, :) = v_cur(:, :)

    Call advect(CuP, CuN, CvP, CvN)

    Do j = 1, ny
      Do i = 1, nx
        div1 = 0.5d0*(u_cur(i, j)+u_cur(i, j+1)-u_cur(i-1, j)-u_cur(i-1, j+1))/dx
        div2 = 0.5d0*(v_cur(i, j+1)-v_cur(i, j-1))/dy
        div = dt*B_cur(i, j)*(div1+div2)  
        advy(i, j)= B_nex(i, j)+div
      End Do
    End Do

    ! 计算ty及pgrdy,与advy相加得到dv
    Do j = 1, ny
      Do i = 1, nx
        ! 计算v的中间变量 
        pgrdy = -dt*g*(eta_cur(i, j+1)-eta_cur(i, j))/dy
        hv    = 0.5d0*(h(i, j)+h(i, j+1))
        vv    = v_cur(i, j)
        uv    = 0.25d0*(u_cur(i, j)+u_cur(i, j+1)+u_cur(i-1, j)+u_cur(i-1, j+1))
        If (hv>0.0) Then
          ty = dt*tauy/(rho*hv)
          Ry(i, j) = dt*r*sqrt(vv**2+uv**2)/hv 
        End If
        dv(i, j) = ty+pgrdy+advy(i, j)

        ! 计算v
        If (mask(i, j)==1) Then
          If ((mask(i, j+1)==1) .or. (dv(i, j)>0.0)) v_nex(i, j) = (v_cur(i, j)+dv(i, j))/(1.0d0+Ry(i, j))
        Else
          If ((mask(i, j+1)==1) .and. (dv(i, j)<0.0)) v_nex(i, j) = (v_cur(i, j)+dv(i, j))/(1.0d0+Ry(i, j))
        End If
      End Do
    End Do
  End subroutine cal_v

  ! 计算eta子程序====================================================================
  Subroutine cal_eta
  Implicit None
    ! 定义cal_eta子程序用到的局部变量
    Real(8) :: CuP(0:nx+1, 0:ny+1), CuN(0:nx+1, 0:ny+1)
    Real(8) :: CvP(0:nx+1, 0:ny+1), CvN(0:nx+1, 0:ny+1)

    Do j = 0, ny+1
      Do i = 0, nx+1
        CuP(i, j) = 0.5d0*(u_nex(i, j)+abs(u_nex(i, j)))*dt/dx
        CuN(i, j) = 0.5d0*(u_nex(i, j)-abs(u_nex(i, j)))*dt/dx
        CvP(i, j) = 0.5d0*(v_nex(i, j)+abs(v_nex(i, j)))*dt/dy
        CvN(i, j) = 0.5d0*(v_nex(i, j)-abs(v_nex(i, j)))*dt/dy
      End Do
    End Do

    B_cur(:, :) = h(:, :)

    Call advect(CuP, CuN, CvP, CvN)

    Do j = 1, ny
      Do i = 1, nx
        ! 计算eta
        eta_nex(i, j) = eta_cur(i, j)+B_nex(i, j)
      End Do
    End Do
  End Subroutine cal_eta

! 计算对流项子程序===========================================================================
  Subroutine advect(CCuP, CCuN, CCvP, CCvN)
    Real(8), Intent(IN) :: CCuP(0:nx+1, 0:ny+1), CCuN(0:nx+1, 0:ny+1), CCvP(0:nx+1, 0:ny+1), CCvN(0:nx+1, 0:ny+1)
    Real(8)             :: RxP(0:nx+1, 0:ny+1), RxN(0:nx+1, 0:ny+1)
    Real(8)             :: RyP(0:nx+1, 0:ny+1), RyN(0:nx+1, 0:ny+1)
    Real(8)             :: dB, term1, term2, term3, term4
    Real(8)             :: BwP, BwN, BeP, BeN, BsP, BsN, BnP, BnN 

    ! 数组初始化
    RxP(:, :) = 0.0d0
    RxN(:, :) = 0.0d0
    RyP(:, :) = 0.0d0
    RyN(:, :) = 0.0d0

    ! 计算x方向，y方向的r+
    Do j = 1, ny
      Do i = 1, nx
        dB =  B_cur(i+1, j)-B_cur(i, j)
        If (abs(dB) > 0.0) RxP(i, j) = (B_cur(i, j)-B_cur(i-1, j))/dB
        dB =  B_cur(i, j+1)-B_cur(i, j)
        If (abs(dB) > 0.0) RyP(i, j) = (B_cur(i, j)-B_cur(i, j-1))/dB
      End Do
    End Do

    ! 计算x方向，y方向的r-
    Do j = 1, ny
      Do i = 0, nx-1
        dB = B_cur(i+1, j)-B_cur(i, j)
        If (abs(dB) > 0.0) RxN(i, j) = (B_cur(i+2, j)-B_cur(i+1, j))/dB
      End Do
    End Do

    Do j = 0, ny-1
      Do i = 1, nx
        dB = B_cur(i, j+1)-B_cur(i, j)
        If (abs(dB) > 0.0) RyN(i, j) = (B_cur(i, j+2)-B_cur(i, j+1))/dB
      End Do
    End Do
 
    ! 计算示踪物B浓度
    Do j = 1, ny
      Do i = 1, nx
        !x方向
        BwP = B_cur(i-1, j)+0.5d0*PSI(RxP(i-1, j), mode)*(1.0d0-CCuP(i-1, j))*(B_cur(i, j)-B_cur(i-1, j))
        BwN = B_cur(i, j)-0.5d0*PSI(RxN(i-1, j), mode)*(1.0d0+CCuN(i-1, j))*(B_cur(i, j)-B_cur(i-1, j))
        BeP = B_cur(i, j)+0.5d0*PSI(RxP(i, j), mode)*(1.0d0-CCuP(i, j))*(B_cur(i+1, j)-B_cur(i, j))
        BeN = B_cur(i+1, j)-0.5d0*PSI(RxN(i, j), mode)*(1.0d0+CCuN(i, j))*(B_cur(i+1, j)-B_cur(i, j))
        ! y方向
        BsP = B_cur(i, j-1)+0.5d0*PSI(RyP(i, j-1), mode)*(1.0d0-CCvP(i, j-1))*(B_cur(i, j)-B_cur(i, j-1))
        BsN = B_cur(i, j)-0.5d0*PSI(RyN(i, j-1), mode)*(1.0d0+CCvN(i, j-1))*(B_cur(i, j)-B_cur(i, j-1))
        BnP = B_cur(i, j)+0.5d0*PSI(RyP(i, j), mode)*(1.0d0-CCvP(i, j))*(B_cur(i, j+1)-B_cur(i, j))
        BnN = B_cur(i, j+1)-0.5d0*PSI(RyN(i, j), mode)*(1.0d0+CCvN(i, j))*(B_cur(i, j+1)-B_cur(i, j))

        term1 = CCuP(i-1, j)*BwP+CCuN(i-1, j)*BwN
        term2 = CCuP(i, j)*BeP+CCuN(i, j)*BeN
        term3 = CCvP(i, j-1)*BsP+CCvN(i, j-1)*BsN
        term4 = CCvP(i, j)*BnP+CCvN(i, j)*BnN

        B_nex(i, j) = term1-term2+term3-term4
      End Do
    End Do

    Contains
      Real(8) Function PSI(rr, mmode)
        Real(8), Intent(IN) :: rr 
        Integer, Intent(IN) :: mmode
        Real(8)             :: comp1, comp2, comp3
    
        If (mmode == 1) PSI = 0.0d0 
        If (mmode == 2) PSI = 1.0d0
        If (mmode == 3) Then
          comp1 = Min(2.0d0*rr, 1.0d0)
          comp2 = Min(rr, 2.0d0)
          comp3 = Max(comp1, comp2)
          PSI   = Max(comp3, 0.0d0)
        End If
      End Function
  End Subroutine advect

End Module

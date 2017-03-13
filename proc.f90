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
    Real(8) :: du(0:nx+1, 0:ny+1), upre(0:nx+1, 0:ny+1)
    Real(8) :: Rx(0:nx+1, 0:ny+1)
    Real(8) :: hu, uu, vu
    Real(8) :: pgrdx, tx, corx
    Real(8) :: advx(0:nx+1, 0:ny+1), div, div1, div2          ! x方向对流项相关变量
    Real(8) :: diffx, term1, term2, term3, term4, s1, s2, h1  ! x方向扩散项相关变量

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
        advx(i, j)= B_nex(i, j)+div  ! 计算出advx
      End Do
    End Do

    ! 分别计算tx，pgrdx及diffx
    Do j = 1, ny
      Do i = 1, nx
        ! 计算u的中间变量    
        pgrdx = -dt*g*(eta_cur(i+1, j)-eta_cur(i, j))/dx  ! 计算出pgrdx
        hu    = 0.5d0*(h(i, j)+h(i+1, j))
        uu    = u_cur(i, j)
        vu    = 0.25d0*(v_cur(i, j)+v_cur(i, j-1)+v_cur(i+1, j)+v_cur(i+1, j-1))
        diffx = 0.0d0
        If (hu>0.0) Then
          tx = dt*taux/(rho*hu)                 ! 计算出tx
          Rx(i, j) = dt*r*sqrt(uu**2+vu**2)/hu  ! 计算出Rx

          ! 根据滑移边界类型计算diffx
          term1 = h(i+1, j)*(u_cur(i+1, j)-u_cur(i, j))/dx
          term2 = h(i, j)*(u_cur(i, j)-u_cur(i-1, j))/dx

          s1 = 1.0d0
          s2 = 1.0d0
          h1 = h(i, j+1)+h(i+1, j+1)
          If (h1<hmin) Then
            s1 = 0.0d0
            s2 = slip
            h1 = h(i, j)+h(i+1, j)
          End If
          term3 = 0.25d0*(h(i, j)+h(i+1, j)+h1)*(s1*u_cur(i, j+1)-s2*u_cur(i, j))/dy

          s1 = 1.0d0
          s2 = 1.0d0
          h1 = h(i, j-1)+h(i+1, j-1)
          If (h1<hmin) Then
            s1 = 0.0d0
            s2 = slip
            h1 = h(i, j)+h(i+1, j)
          End If
          term4 = 0.25d0*(h(i, j)+h(i+1, j)+h1)*(s2*u_cur(i, j)-s1*u_cur(i, j-1))/dy

          diffx = dt*ah*((term1-term2)/dx+(term3-term4)/dy)/hu  ! 计算出diffx
        End If
        
        ! 第一步，得到不含科氏力项预估的upre
        upre(i, j) = u_cur(i, j)+tx+pgrdx+advx(i, j)+diffx

        ! 第二步，半隐式方法求解科氏力项
        corx     = dt*f*vu
        du(i, j) = (upre(i, j)-beta*u_cur(i, j)+corx)/(1.0d0+beta)-u_cur(i, j)  ! 计算出du

        ! 第三步，计算u
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
    Real(8) :: dv(0:nx+1, 0:ny+1), vpre(0:nx+1, 0:ny+1)
    Real(8) :: Ry(0:nx+1, 0:ny+1)
    Real(8) :: hv, vv, uv
    Real(8) :: pgrdy, ty, cory
    Real(8) :: advy(0:nx+1, 0:ny+1), div, div1, div2          ! y方向对流项相关变量
    Real(8) :: diffy, term1, term2, term3, term4, s1, s2, h1  ! y方向扩散项相关变量

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
        advy(i, j)= B_nex(i, j)+div  ! 计算出advy
      End Do
    End Do

    ! 分别计算ty，pgrdy及diffy
    Do j = 1, ny
      Do i = 1, nx
        ! 计算v的中间变量 
        pgrdy = -dt*g*(eta_cur(i, j+1)-eta_cur(i, j))/dy  ! 计算出pgrdy
        hv    = 0.5d0*(h(i, j)+h(i, j+1))
        vv    = v_cur(i, j)
        uv    = 0.25d0*(u_cur(i, j)+u_cur(i, j+1)+u_cur(i-1, j)+u_cur(i-1, j+1))
        diffy = 0.0d0
        If (hv>0.0) Then
          ty = dt*tauy/(rho*hv)  ! 计算出ty
          Ry(i, j) = dt*r*sqrt(vv**2+uv**2)/hv  ! 计算出Ry

          ! 根据滑移边界类型计算diffy
          s1 = 1.0d0
          s2 = 1.0d0
          h1 = h(i+1, j)+h(i+1, j+1)
          If (h1<hmin) Then
            s1 = 0.0d0
            s2 = slip
            h1 = h(i, j)+h(i, j+1)
          End If
          term1 = 0.25d0*(h(i, j)+h(i, j+1)+h1)*(s1*v_cur(i+1, j)-s2*v_cur(i, j))/dx

          s1 = 1.0d0
          s2 = 1.0d0
          h1 = h(i-1, j)+h(i-1, j+1)
          If (h1<hmin) Then
            s1 = 0.0d0
            s2 = slip
            h1 = h(i, j)+h(i, j+1)
          End If
          term2 = 0.25d0*(h(i, j)+h(i, j+1)+h1)*(s2*v_cur(i, j)-s1*v_cur(i-1, j))/dx

          term3 = h(i, j+1)*(v_cur(i, j+1)-v_cur(i, j))/dy
          term4 = h(i, j)*(v_cur(i, j)-v_cur(i-1, j-1))/dy
          
          diffy = dt*ah*((term1-term2)/dx+(term3-term4)/dy)/hv  ! 计算出diffy
        End If

        ! 第一步，得到不含科氏力项预估的vpre
        vpre(i, j) = v_cur(i, j)+ty+pgrdy+advy(i, j)+diffy

        ! 第二步，半隐式方法求解科氏力项
        cory     = -dt*f*uv
        dv(i, j) = (vpre(i, j)-beta*v_cur(i, j)+cory)/(1.0d0+beta)-v_cur(i, j)  ! 计算出dv

        ! 计算v
        If (mask(i, j)==1) Then
          If ((mask(i, j+1)==1) .or. (dv(i, j)>0.0)) v_nex(i, j) = (v_cur(i, j)+dv(i, j))/(1.0d0+Ry(i, j))
        Else
          If ((mask(i, j+1)==1) .and. (dv(i, j)<0.0)) v_nex(i, j) = (v_cur(i, j)+dv(i, j))/(1.0d0+Ry(i, j))
        End If

      End Do
    End Do
  End subroutine cal_v

  ! 计算eta子程序=========================================================================
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

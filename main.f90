Program advection
!******************************************************************
! 二维浅水方程，考虑了对流项，扩散项，科氏力，风应力，底摩擦，水平压强梯度力，干湿
! 网格算法。
! 本程序为主程序（main.f90），还有变量声明模块（var.f90）、
! 计算过程模块（proc.f90）以及水深文件topo1.dat
! rewritten by lihu, 2017.3
!******************************************************************

Use var
Use proc
Implicit None

  ! 变量初始化，设置初始条件*************************************************************
  dx   = 100.0d0     ! 水平空间步长
  dy   = 100.0d0     ! 水平空间步长
  dt   = 3.0d0       ! 时间步长

  nx   = 51
  ny   = 51
  nt   = Int(10*24.0*3600.0/dt)
  nout = Int(2.0*3600.0/dt)

  hmin  = 0.10d0      ! 最小截断水深（m）
  rho   = 1028.0d0    ! 海水密度
  theta = 35.0d0      ! 区域所在纬度
  f     = 2.0d0*omega*sin(theta*pi/180.0d0)  ! 科氏力参数
  beta  = (0.5d0*dt*f)**2  ! beta系数
  r     = 0.001d0     ! 底摩擦系数
  ah    = 1.0d0       ! 水平湍流粘性系数
!  eps  = 0.05d0      ! Shapiro滤波参数


  mode = 3         ! 对流项差分方法选项（1为迎风格式；2为Lax-Wendroff格式；3为Superbee格式）
  slip = 2.0d0     ! 滑移边界类型选项（0为full-slip类型；1为full-slip类型；2为no-slip类型）

  ! 分配动态数组
  Allocate(mask(0:nx+1, 0:ny+1))     ! 水陆掩膜，判断干湿格点（1为湿点，0为干点）
  Allocate(hs(0:nx+1, 0:ny+1))       ! 净水深
  Allocate(h(0:nx+1, 0:ny+1))        ! 实时水深
  Allocate(eta_cur(0:nx+1, 0:ny+1))  ! 第n步水位
  Allocate(eta_nex(0:nx+1, 0:ny+1))  ! 第n+1步水位
  Allocate(u_cur(0:nx+1, 0:ny+1))    ! 第n步u流速
  Allocate(u_nex(0:nx+1, 0:ny+1))    ! 第n+1步u流速
  Allocate(v_cur(0:nx+1, 0:ny+1))    ! 第n步v流速
  Allocate(v_nex(0:nx+1, 0:ny+1))    ! 第n+1步v流速

  ! 是否分配示踪物数组，这里表示u或v
  Allocate(B_cur(0:nx+1, 0:ny+1))    ! 第n步示踪物B浓度
  Allocate(B_nex(0:nx+1, 0:ny+1))    ! 第n+1步示踪物B浓度

  ! 读取净水深文件
  Open(10, file ='./topo1.dat', form = 'formatted')
  Do j = 1, ny
    Read(10, *) hs(1:nx, j)
  End Do
  Close(10)
  ! 设置外围边界的净水深
  hs(0, 0:ny+1)    = -10.0d0
  hs(nx+1, 0:ny+1) = -10.0d0
  hs(0:nx+1, 0)    = -10.0d0
  hs(0:nx+1, ny+1) = -10.0d0

  ! 设置eta初始条件
  Do j = 0, ny+1
    Do i = 0, nx+1
      eta_cur(i, j) = -Min(0.0d0, hs(i, j))
      eta_nex(i, j) = eta_cur(i, j)
    End Do
  End Do

  Do j = 0, ny+1
    Do i = 0, nx+1
      h(i, j)    = hs(i, j)+eta_cur(i, j)  ! 设置h初始条件
      mask(i, j) = 1
      If (h(i, j)<hmin) mask(i, j) = 0  ! 设置mask初始条件
    End Do
  End Do

  ! 设置u，v初始条件
  u_cur = 0.0d0
  u_nex = 0.0d0
  v_cur = 0.0d0
  v_nex = 0.0d0

  ! 新建并打开eta，u，v输出文件*********************************************************
  Open(10, file = './output/eta.dat', form = 'formatted')
  Open(20, file = './output/u.dat', form = 'formatted')
  Open(30, file = './output/v.dat', form = 'formatted')
  ! 计算过程*************************************************************************
  Do n = 1, nt
    time = Dble(n)*dt
    taux = 0.0d0
    tauy = 0.2d0*Min(time/(1.0d0*24.0d0*3600.0d0), 1.0d0)

    ! 调用计算过程子程序 
    Call cal_u
    Call cal_v
    Call cal_eta

    ! 调用Shapiro滤波子程序
    ! Call shapiro

    ! 迭代变量，以第n+1步的计算结果替换第n步的eta，h，mask，u，v
    Do j = 0, ny+1
      Do i = 0, nx+1
        eta_cur(i, j) = eta_nex(i, j)
        h(i, j)       = hs(i, j)+eta_cur(i, j)
        mask(i, j)    = 1
        If (h(i, j)<hmin) mask(i, j) = 0
        u_cur(i, j)   = u_nex(i, j)
        v_cur(i, j)   = v_nex(i, j)
      End Do
    End Do

    ! 输出计算结果到文本
    If (Mod(n, nout)==0) Then
      Do j = 0, ny+1
        Write(10, '(<nx+2>F12.6)') eta_cur(:, j)
        Write(20, '(<nx+2>F12.6)') u_cur(:, j)
        Write(30, '(<nx+2>F12.6)') v_cur(:, j)
      End Do
      Write(*, *) "Output data at time = ", time/(24.0*3600.0), "days"
    End If
  End Do

  ! 关闭输出文件***********************************************************************
  Close(10)
  Close(20)
  Close(30)
End Program

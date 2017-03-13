Module var
Implicit None
  ! 常量定义
  Integer :: nx, ny, nt, nout                            ! 水平空间步数, 时间歩数, 输出间隔（每隔多少步输出结果）
  Real(8), Parameter :: pi=3.14159265359d0, g=9.81d0     ! 常量参数
  Real(8), Parameter :: omega = 2.0d0*pi/(24.d0*3600.d0) ! 地球自转角速度，求f时使用
  Real(8) :: dx, dy, dt, hmin, rho                       ! 水平空间步长，时间歩长，最小截断水深，密度
  Real(8) :: theta, f, beta                              ! 区域所在纬度，科氏力参数，beta系数
  Real(8) :: r, taux, tauy, ah                           ! 底摩擦系数，风阻系数, 水平涡粘系数
  Real(8) :: time, slip                                  ! 时间变量，滑移边界类型
  Real(8) :: eps                                         ! Shapiro滤波参数

  ! 变量定义
  Integer              :: i, j, n                       ! x轴，y轴，时间轴的循环计数器
  Integer              :: mode                          ! 对流项差分格式选项
  Integer, Allocatable :: mask(:, :)                    ! 水陆掩膜，判断干湿格点（1为湿点，0为干点）
  Real(8), Allocatable :: hs(:, :), h(:, :)             ! 净水深，实时水深
  Real(8), Allocatable :: eta_cur(:, :), eta_nex(:, :)  ! 第n步水位，第n+1步水位
  Real(8), Allocatable :: u_cur(:, :), u_nex(:, :)      ! 第n步u流速，第n+1步u流速
  Real(8), Allocatable :: v_cur(:, :), v_nex(:, :)      ! 第n步v流速，第n+1步v流速
  Real(8), Allocatable :: B_cur(:, :), B_nex(:, :)      ! 第n步示踪物B的浓度，第n+1步示踪物B的浓度
End Module

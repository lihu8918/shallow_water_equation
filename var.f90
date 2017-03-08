Module var
Implicit None
! 常量定义
Integer :: nx, ny, nt, nout  ! 水平空间步数, 时间歩数, 输出间隔（每隔多少步输出结果）
Real(8) :: dx, dy, dt  ! 水平空间步长及时间歩长
Real(8) :: hmin, rho, g, r, taux, tauy  ! 最小截断水深，海水密度，重力加速度，底摩擦系数，风阻系数
Real(8) :: time  ! 时间变量

! 变量定义
Integer              :: i, j, n  ! x轴，y轴，时间轴的循环计数器
Integer              :: mode     ! 对流项差分格式选项
Integer, Allocatable :: mask(:, :)  ! 水陆掩膜，判断干湿格点（1为湿点，0为干点）
Real(8), Allocatable :: hs(:, :), h(:, :)  ! 净水深，实时水深
Real(8), Allocatable :: eta_cur(:, :), eta_nex(:, :)  ! 第n步水位，第n+1步水位
Real(8), Allocatable :: u_cur(:, :), u_nex(:, :)  ! 第n步u流速，第n+1步u流速
Real(8), Allocatable :: v_cur(:, :), v_nex(:, :)  ! 第n步v流速，第n+1步v流速
Real(8), Allocatable :: B_cur(:, :), B_nex(:, :)  ! 第n步示踪物B的浓度，第n+1步示踪物B的浓度

Real(8)              :: eps  ! Shapiro滤波参数
End Module

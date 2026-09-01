# SailModule

## 1. 模块用途

本目录实现一个独立、宿主无关的单风帆 Fortran 气动载荷模块。模块以一组一维 CFD 气动数据库为依据，输入船体坐标系中的风、船体运动、帆气动中心位置和帆角，输出单帆六自由度气动载荷。它面向风帆助推船波浪载荷研究中的时域风帆载荷计算接口，也可用于论文第三章的方法说明、后续维护复现以及单模块和端到端测试。

当前代码尚未接入 `cwbasic`、WALCS、WALCS-LE、`planning`、`MAINCAL` 或任何时域积分主循环；本目录只完成独立模块级计算链。

## 2. 坐标系与六自由度约定

船体坐标系为右手系：`+x_b` 指向船首，`+y_b` 指向左舷，`+z_b` 向上。六自由度载荷顺序为：

1. `Fx`：Surge；
2. `Fy`：Sway；
3. `Fz`：Heave；
4. `Mx`：Roll；
5. `My`：Pitch；
6. `Mz`：Yaw。

```text
LOAD_6DOF = (Fx, Fy, Fz, Mx, My, Mz)
```

## 3. 当前物理模型范围

当前模型仅以水平相对风计算攻角、动压及升阻力方向。完整三维相对风 `V_REL_BODY` 仍被计算和保留，但其垂向分量不进入攻角、动压或升阻力方向。因而当前 `Fz = 0`；当作用点 `R_SAIL_BODY` 具有非零坐标（包括非零高度）时，`R_SAIL_BODY × FORCE_BODY` 仍可产生 `Mx`、`My` 和 `Mz`。

模型目前为单风帆、单组 `alpha-CL-CD` 一维数据库及分段线性插值；不包含动态失速、三维垂向气动力、多帆相互干扰或帆结构变形。

## 4. 风帆与攻角约定

帆角 `DELTA_S_DEG = 0°` 时，弦线为：

```text
c = (0, -1, 0)
```

任意帆角 `delta` 的弦线为：

```text
c(delta) = (-sin(delta), -cos(delta), 0)
```

帆弦线按无向轴线处理，具有 180° 周期。帆气动中心相对船体重心的位置由宿主作为 `R_SAIL_BODY` 提供；物理上可写为 `r_sail = r_CE - r_CG`。代码中的相对风关系为：

```text
V_SAIL_POINT_BODY = V_CG_BODY + OMEGA_BODY × R_SAIL_BODY
V_REL_BODY        = V_WIND_BODY - V_SAIL_POINT_BODY
                  = V_WIND_BODY - V_CG_BODY - OMEGA_BODY × R_SAIL_BODY
```

令 `V_REL_H_BODY` 为 `V_REL_BODY` 的水平投影，则：

```text
e_upstream = -V_REL_H_BODY / |V_REL_H_BODY|

alpha_raw = atan2(c_x e_upstream,y - c_y e_upstream,x,
                  c_x e_upstream,x + c_y e_upstream,y)

alpha_db = MODULO(alpha_raw, 180°)
0° <= alpha_db < 180°
```

`DELTA_S_DEG`、数据库角和接口角输出均使用度（deg）；三角函数内部将角度转换为弧度（rad）。接近 0° 或 180°、且落入 `ANGLE_TOL_DEG` 的周期映射结果会归零。

## 5. 升阻力和载荷公式

```text
e_D = V_REL_H_BODY / |V_REL_H_BODY|
e_L = e_D × e_z = (e_D,y, -e_D,x, 0)

q = 0.5 RHO_AIR |V_REL_H_BODY|^2
F = q SAIL_AREA (CD e_D + CL e_L)
M = R_SAIL_BODY × F
```

`CL` 是有符号量；`CD` 应非负，低于 `-DATABASE_TOL` 的负值会以 `SAIL_ERR_INVALID_INPUT` 拒绝。当前常量为 `RHO_AIR = 1.225 kg/m^3`、`SAIL_CHORD = 2.5 m`、`SAIL_SPAN = 6.0 m`、`SAIL_AREA = 15.0 m^2`。

STAR 数据库基准约定为：当 `V_REL_H_BODY = (0,+V,0)` 时，`e_D=(0,1,0)`、`e_L=(1,0,0)`，因此 `Fx=qSCL`、`Fy=qSCD`；该符号约定与当前数据库定义一致。

## 6. 目录结构

```text
SailModule/
|-- SAILPARAM_MOD.F90              参数、坐标索引、容差和错误码
|-- SAILDATABASE_MOD.F90           数据库读取、验证、查询和清理
|-- SAILRELWIND_MOD.F90            帆作用点速度和相对风
|-- SAILANGLE_MOD.F90              帆弦线、上游方向和攻角映射
|-- SAILINTERP_MOD.F90             一维分段线性系数插值
|-- SAILFORCE_MOD.F90              升阻力、力矩和六自由度映射
|-- SAILMODULE_API_MOD.F90         初始化、总计算和清理 API
|-- sail_database.dat              alpha_deg、CL、CD 数据库
|-- test_saildatabase.f90          数据库单元测试
|-- test_sailrelwind.f90           相对风单元测试
|-- test_sailangle.f90             攻角单元测试
|-- test_sailinterp.f90            插值单元测试
|-- test_sailforce.f90             气动力单元测试
|-- test_sailmodule.f90            端到端与生命周期测试
|-- README.md                      使用、物理定义和构建说明
`-- INTERFACE_ACCEPTANCE_REPORT.md 独立接口验收记录
```

## 7. 模块调用关系

```text
InitializeSailModule
        ↓
ComputeSailModuleLoads
        |-- ComputeSailRelativeWind
        |-- ComputeSailAngle
        |-- GetSailCoeff
        `-- ComputeSailForce
        ↓
FinalizeSailModule
```

数据库只在初始化阶段读取。`ComputeSailModuleLoads` 不自动加载数据库，可在未来时域循环中重复调用。

## 8. PUBLIC 接口

以下参数顺序与源码一致。除特别注明外，`IERR` 为无量纲整数状态，`MESSAGE` 为诊断文本。

### 8.1 `SAILDATABASE_MOD`

| 过程 | 主要输入 | 主要输出 | 作用 |
|---|---|---|---|
| `ReadSailDatabase(FILE_NAME, IERR, MESSAGE)` | 可选路径 `FILE_NAME` | `IERR`, `MESSAGE` | 读取并验证数据库；缺省路径为 `sail_database.dat` |
| `ValidateSailDatabase(IERR, MESSAGE)` | 模块内数据库 | `IERR`, `MESSAGE` | 验证节点数、有限性、CD、顺序和端点 |
| `ClearSailDatabase()` | 无 | 无 | 释放数据库；可重复调用 |
| `IsSailDatabaseInitialized()` | 无 | 逻辑值 | 查询已初始化状态 |
| `GetSailDatabaseSize()` | 无 | 节点数 | 未初始化时返回 0 |
| `GetSailDatabaseNode(INDEX_NODE, ALPHA_DEG, CL, CD, IERR, MESSAGE)` | `INDEX_NODE` | `ALPHA_DEG`(deg), `CL`, `CD` | 复制指定节点 |
| `GetSailDatabaseBounds(ALPHA_MIN, ALPHA_MAX, IERR, MESSAGE)` | 无 | `ALPHA_MIN`, `ALPHA_MAX`(deg) | 返回首末攻角 |

### 8.2 计算模块

| 模块/过程 | 按源码顺序的输入 | 主要输出 | 作用 |
|---|---|---|---|
| `SAILRELWIND_MOD::ComputeSailRelativeWind` | `V_WIND_BODY`, `V_CG_BODY` (m/s), `OMEGA_BODY` (rad/s), `R_SAIL_BODY` (m) | `V_SAIL_POINT_BODY`, `V_REL_BODY`, `V_REL_H_BODY` (m/s), `V_REL_H_MAG` (m/s) | 计算帆点速度和相对风 |
| `SAILANGLE_MOD::ComputeSailAngle` | `V_REL_H_BODY` (m/s), `V_REL_H_MAG` (m/s), `DELTA_S_DEG` (deg) | `DELTA_NORMALIZED_DEG`, `C_CHORD_BODY`, `E_UPSTREAM_BODY`, `ALPHA_RAW_DEG`, `ALPHA_DB_DEG` | 计算弦线和数据库攻角 |
| `SAILINTERP_MOD::GetSailCoeff` | `ALPHA_QUERY_DEG` (deg) | `CL_OUT`, `CD_OUT` | 查询节点或分段线性插值；不外推 |
| `SAILFORCE_MOD::ComputeSailForce` | `V_REL_H_BODY`, `V_REL_H_MAG`, `CL`, `CD`, `R_SAIL_BODY` | `E_DRAG_BODY`, `E_LIFT_BODY`, `Q_DYNAMIC` (Pa), `FORCE_BODY` (N), `MOMENT_BODY` (N·m), `LOAD_6DOF` | 计算力、力矩和六自由度载荷 |

上述各计算过程的末尾均依次含 `IERR, MESSAGE`。

### 8.3 `SAILMODULE_API_MOD`

| 过程 | 参数与作用 |
|---|---|
| `InitializeSailModule(DATABASE_FILE, IERR, MESSAGE)` | 可选数据库路径、状态、消息；完成一次数据库初始化 |
| `ComputeSailModuleLoads(...)` | 依次接收 5 个物理输入，返回完整中间量、系数、力、力矩和 `LOAD_6DOF`，末尾为 `IERR, MESSAGE` |
| `FinalizeSailModule(IERR, MESSAGE)` | 清理数据库存储 |

`ComputeSailModuleLoads` 的准确参数顺序为：

```text
V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, DELTA_S_DEG,
V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, V_REL_H_MAG,
DELTA_NORMALIZED_DEG, C_CHORD_BODY, E_UPSTREAM_BODY,
ALPHA_RAW_DEG, ALPHA_DB_DEG, CL, CD, E_DRAG_BODY,
E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, LOAD_6DOF,
IERR, MESSAGE
```

## 9. API 使用生命周期

下面只展示调用结构，不假设任何 WALCS 变量名：

```fortran
CALL InitializeSailModule(IERR=IERR, MESSAGE=MESSAGE)
IF (IERR /= SAIL_OK) STOP

DO IT = 1, NT
  CALL ComputeSailModuleLoads( &
    V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, DELTA_S_DEG, &
    V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, V_REL_H_MAG, &
    DELTA_NORMALIZED_DEG, C_CHORD_BODY, E_UPSTREAM_BODY, &
    ALPHA_RAW_DEG, ALPHA_DB_DEG, CL, CD, E_DRAG_BODY, &
    E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, LOAD_6DOF, &
    IERR, MESSAGE)
END DO

CALL FinalizeSailModule(IERR, MESSAGE)
```

初始化通常只执行一次，计算可重复执行，清理在结束阶段执行。总计算接口不会自动加载数据库。

## 10. 输入输出单位表

| 变量 | 含义 | 单位 |
|---|---|---|
| `V_WIND_BODY` | 船体坐标系风速 | m/s |
| `V_CG_BODY` | 船体重心平动速度 | m/s |
| `OMEGA_BODY` | 船体角速度 | rad/s |
| `R_SAIL_BODY` | 重心至帆气动中心位置矢量 | m |
| `DELTA_S_DEG` | 输入帆角 | deg |
| `V_SAIL_POINT_BODY` | 帆气动中心速度 | m/s |
| `V_REL_BODY` | 完整三维相对风 | m/s |
| `V_REL_H_BODY` | 水平相对风 | m/s |
| `V_REL_H_MAG` | 水平相对风速标量 | m/s |
| `ALPHA_RAW_DEG` | 有符号原始攻角 | deg |
| `ALPHA_DB_DEG` | 180° 周期映射后的攻角 | deg |
| `CL` | 升力系数 | 无量纲 |
| `CD` | 阻力系数 | 无量纲 |
| `Q_DYNAMIC` | 动压 | Pa |
| `FORCE_BODY` | 船体轴气动力 `(Fx,Fy,Fz)` | N |
| `MOMENT_BODY` | 关于船体重心的力矩 `(Mx,My,Mz)` | N·m |
| `LOAD_6DOF` | `(Fx,Fy,Fz,Mx,My,Mz)` | N；N·m |

## 11. 错误码与异常处理

| 名称 | 值 | 含义 |
|---|---:|---|
| `SAIL_OK` | 0 | 成功 |
| `SAIL_ERR_FILE_NOT_FOUND` | 1 | 数据库文件不存在 |
| `SAIL_ERR_FILE_OPEN` | 2 | 文件打开失败 |
| `SAIL_ERR_FILE_READ` | 3 | 文件读取、解析、分配或关闭失败 |
| `SAIL_ERR_DATABASE_EMPTY` | 4 | 无有效节点或少于两个节点 |
| `SAIL_ERR_DATABASE_ORDER` | 5 | 攻角节点次序错误 |
| `SAIL_ERR_DATABASE_DUPLICATE` | 6 | 重复攻角节点 |
| `SAIL_ERR_DATABASE_RANGE` | 7 | 数据库端点范围错误 |
| `SAIL_ERR_DATABASE_NOT_INITIALIZED` | 8 | 数据库未初始化 |
| `SAIL_ERR_INVALID_INPUT` | 9 | 非法输入、非有限值、负 CD 等 |
| `SAIL_ERR_LOW_WIND_SPEED` | 10 | 水平相对风低于 `VREL_MIN` |
| `SAIL_ERR_INTERPOLATION` | 11 | 插值内部状态或结果错误 |

所有公开计算入口均拒绝相关输入中的 NaN/Inf。低风速状态是非致命、可由宿主识别的状态；总接口将本时间步下游载荷置零，不应因此终止整个时域计算。源码判定条件为 `V_REL_H_MAG < VREL_MIN`，其中 `VREL_MIN = 1.0E-8 m/s`；恰等于阈值时不是低风速错误。

## 12. 数据库格式

`sail_database.dat` 是无 BOM 的 ASCII 文本，严格使用三列 `alpha_deg CL CD`。空行以及首个非空字符为 `#` 或 `!` 的整行注释可被忽略；首个候选内容也可为同时包含 `alpha`、`cl`、`cd`（不区分 ASCII 大小写）的表头。行内注释、尾随第四字段或畸形行会被拒绝。

节点必须严格递增、不得重复，值必须有限，`CD` 不得显著为负，端点必须覆盖 0° 和 180°。当前数据库有 33 个节点，范围为 0° 至 180°：

| 位置 | `alpha_deg` | `CL` | `CD` |
|---|---:|---:|---:|
| 首节点 | 0 | -0.29395224 | 0.17179593 |
| 80° 节点 | 80 | -0.27873644 | 1.43531615 |
| 末节点 | 180 | 0.26555131 | 0.17291675 |

插值模块可直接查询 180° 端点；姿态映射模块的输出范围则为 `[0°,180°)`。不要修改数据库以改变这一区别。

## 13. 编译环境

已验证环境：Intel ifx 2024.2.0 Build 20240602，oneAPI 位于 `D:\intel\oneAPI`，Visual Studio 2022 位于 `D:\englishroad\vsstudio\product`。

在 `cmd.exe` 语法中设置环境时，`VS2022INSTALLDIR` 的赋值末尾与 `&&` 之间不要插入会进入变量值的空格，例如：

```bat
set "VS2022INSTALLDIR=D:\englishroad\vsstudio\product"&& call "D:\intel\oneAPI\setvars.bat" intel64 vs2022
```

## 14. 单模块编译示例

以下命令在 `SailModule` 目录执行，并遵循真实模块依赖顺序。

数据库测试：

```bat
ifx /stand:f18 /warn:all /check:all /traceback SAILPARAM_MOD.F90 SAILDATABASE_MOD.F90 test_saildatabase.f90 /exe:test_saildatabase.exe
```

相对风测试：

```bat
ifx /stand:f18 /warn:all /check:all /traceback SAILPARAM_MOD.F90 SAILRELWIND_MOD.F90 test_sailrelwind.f90 /exe:test_sailrelwind.exe
```

攻角测试：

```bat
ifx /stand:f18 /warn:all /check:all /traceback SAILPARAM_MOD.F90 SAILANGLE_MOD.F90 test_sailangle.f90 /exe:test_sailangle.exe
```

插值测试：

```bat
ifx /stand:f18 /warn:all /check:all /traceback SAILPARAM_MOD.F90 SAILDATABASE_MOD.F90 SAILINTERP_MOD.F90 test_sailinterp.f90 /exe:test_sailinterp.exe
```

气动力测试：

```bat
ifx /stand:f18 /warn:all /check:all /traceback SAILPARAM_MOD.F90 SAILFORCE_MOD.F90 test_sailforce.f90 /exe:test_sailforce.exe
```

端到端测试：

```bat
ifx /stand:f18 /warn:all /check:all /traceback SAILPARAM_MOD.F90 SAILDATABASE_MOD.F90 SAILRELWIND_MOD.F90 SAILANGLE_MOD.F90 SAILINTERP_MOD.F90 SAILFORCE_MOD.F90 SAILMODULE_API_MOD.F90 test_sailmodule.f90 /exe:test_sailmodule.exe
```

## 15. 测试结果摘要

现有测试程序及已完成的真实测试记录如下。数据库测试源码明确编号 1 至 15，因此计数可核验为 15。

| 测试程序 | PASS | FAIL | 状态 |
|---|---:|---:|---|
| `test_saildatabase` | 15 | 0 | PASS |
| `test_sailinterp` | 16 | 0 | PASS |
| `test_sailrelwind` | 16 | 0 | PASS |
| `test_sailangle` | 20 | 0 | PASS |
| `test_sailforce` | 27 | 0 | PASS |
| `test_sailmodule` | 23 | 0 | PASS |

端到端测试包括初始化/清理、STAR 基准、非零作用点力矩、平动和转动影响、非节点插值、低风、NaN/Inf、输入不修改、连续 1000 次一致调用和重新初始化。

## 16. 构建产物清理

可删除的本地编译产物为 `*.mod`、`*.obj`、`*.exe`、`*.pdb`：

```powershell
Remove-Item *.obj, *.mod, *.exe, *.pdb -ErrorAction SilentlyContinue
```

不要删除 `*.F90`、`*.f90`、`sail_database.dat`、`README.md` 或 `INTERFACE_ACCEPTANCE_REPORT.md`。

## 17. Git 建议

应提交源码和文档，不应提交编译产物。提交前可检查：

```powershell
git status --short
```

本说明不代表任何提交操作已经执行。

## 18. 当前限制与后续工作

下一阶段不应继续改变已经验收的独立模块物理公式，而应开展宿主适配设计：

1. 确定接入对象：`cwbasic`/WALCS、WALCS-LE 或 `planning`；
2. 从宿主取得 `V_WIND_BODY`、`V_CG_BODY`、`OMEGA_BODY`、`R_SAIL_BODY` 和 `DELTA_S_DEG`；
3. 将 `LOAD_6DOF` 写入宿主外载荷位置；
4. 核对坐标、单位和时间步一致性；
5. 先做单向接入验证，再进入时域风浪耦合。

当前限制仍是单风帆、一维数据库、分段线性插值和水平气动力；无动态失速、多帆干扰、宿主集成或时域耦合验证。

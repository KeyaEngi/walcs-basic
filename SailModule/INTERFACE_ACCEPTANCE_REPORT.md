# SailModule 独立接口验收报告

## 1. 验收范围

本次验收对象包括参数模块、数据库模块、相对风模块、攻角模块、插值模块、气动力模块、API 总接口、对应单元测试和端到端测试。

本次验收不包括 WALCS、WALCS-LE、`planning` 接入，不包括多风帆模型或时域积分；也不评价尚未实现的动态失速、帆结构变形或三维垂向气动力。

## 2. 文件清单

正式源码：

- `SAILPARAM_MOD.F90`
- `SAILDATABASE_MOD.F90`
- `SAILRELWIND_MOD.F90`
- `SAILANGLE_MOD.F90`
- `SAILINTERP_MOD.F90`
- `SAILFORCE_MOD.F90`
- `SAILMODULE_API_MOD.F90`

测试源码：

- `test_saildatabase.f90`
- `test_sailrelwind.f90`
- `test_sailangle.f90`
- `test_sailinterp.f90`
- `test_sailforce.f90`
- `test_sailmodule.f90`

数据库：`sail_database.dat`。

## 3. 接口验收

| 模块 | PUBLIC 接口 | 职责 | 依赖 | 状态 |
|---|---|---|---|---|
| `SAILPARAM_MOD` | `DP`、数学/物理常量、容差、索引、错误码、数据库范围和默认文件名 | 统一参数定义 | 无 | PASS |
| `SAILDATABASE_MOD` | `ReadSailDatabase`; `ValidateSailDatabase`; `ClearSailDatabase`; `IsSailDatabaseInitialized`; `GetSailDatabaseSize`; `GetSailDatabaseNode`; `GetSailDatabaseBounds` | 数据库生命周期、验证和只读查询 | `SAILPARAM_MOD` | PASS |
| `SAILRELWIND_MOD` | `ComputeSailRelativeWind` | 帆点速度、三维及水平相对风 | `SAILPARAM_MOD` | PASS |
| `SAILANGLE_MOD` | `ComputeSailAngle` | 弦线、上游方向、攻角及 180° 映射 | `SAILPARAM_MOD` | PASS |
| `SAILINTERP_MOD` | `GetSailCoeff` | 节点查询和分段线性插值 | 参数、数据库模块 | PASS |
| `SAILFORCE_MOD` | `ComputeSailForce` | 升阻力、力矩、六自由度载荷 | `SAILPARAM_MOD` | PASS |
| `SAILMODULE_API_MOD` | `InitializeSailModule`; `ComputeSailModuleLoads`; `FinalizeSailModule` | 编排完整生命周期和计算链 | 上述全部计算模块 | PASS |

接口名、参数顺序及 `INTENT` 已与当前源码逐项核对。源码与任务文字的一个需澄清差异是：数据库本身包含并允许直接查询 180° 端点，而 `ComputeSailAngle` 的周期映射输出为 `[0°,180°)`；两者并不冲突。

## 4. 物理定义验收

| 验收项 | 核验结果 | 状态 |
|---|---|---|
| 坐标系 | 右手系；`+x_b` 船首、`+y_b` 左舷、`+z_b` 向上 | PASS |
| 六自由度顺序 | `(Fx,Fy,Fz,Mx,My,Mz)` | PASS |
| 相对风 | `V_REL_BODY = V_WIND_BODY - V_CG_BODY - OMEGA_BODY × R_SAIL_BODY` | PASS |
| `Omega×r` 顺序 | 源码先算 `OMEGA_BODY × R_SAIL_BODY`，再加至帆点速度 | PASS |
| 帆弦线 | `c=(-sin(delta),-cos(delta),0)`；0° 时 `(0,-1,0)` | PASS |
| 上游方向 | `-V_REL_H_BODY / V_REL_H_MAG` | PASS |
| 攻角 `ATAN2` | `atan2(c_x e_up,y-c_y e_up,x, c_x e_up,x+c_y e_up,y)` | PASS |
| 180° 映射 | `MODULO(alpha_raw,180°)`，容差邻域归零 | PASS |
| 阻力方向 | `e_D=V_REL_H_BODY/V_REL_H_MAG` | PASS |
| 升力方向 | `e_L=e_D×e_z=(e_D,y,-e_D,x,0)` | PASS |
| 动压 | `q=0.5 RHO_AIR V_REL_H_MAG^2` | PASS |
| 气动力 | `F=q SAIL_AREA (CD e_D+CL e_L)`，并令 `Fz=0` | PASS |
| 力矩 | `M=R_SAIL_BODY×F` | PASS |
| 六自由度映射 | 力三分量后接力矩三分量 | PASS |

## 5. 数值与边界验收

| 验收项 | 核验结果 | 状态 |
|---|---|---|
| NaN/Inf | 相关公开入口拒绝非有限输入，失败输出按过程约定归零 | PASS |
| 未初始化数据库 | 查询/插值/总计算返回 `SAIL_ERR_DATABASE_NOT_INITIALIZED` | PASS |
| 重复初始化 | 同文件返回成功并跳过重载 | PASS |
| 清理 | 释放存储；重复清理安全 | PASS |
| 重新初始化 | 清理或失败后可重新加载有效数据库 | PASS |
| 低风速 | `< VREL_MIN` 返回可识别非致命状态；总接口载荷为零 | PASS |
| 等于 `VREL_MIN` | 源码使用严格小于号；阈值本身继续计算 | PASS |
| 越界攻角 | 超出含容差数据库范围时拒绝，不外推 | PASS |
| 数据库端点 | 0°、180° 均可精确查询 | PASS |
| 非节点插值 | 使用实际相邻节点做分段线性插值 | PASS |
| 负 `CL` | 允许，符号进入升力项 | PASS |
| 负 `CD` | 低于 `-DATABASE_TOL` 时拒绝 | PASS |
| 输入不修改 | 5 个总接口输入在测试后保持不变 | PASS |
| 重复调用 | 总接口连续 1000 次结果一致 | PASS |
| 输出有限性 | 计算后统一检查，溢出或非有限结果被拒绝 | PASS |

## 6. 编译验收

- 编译器：Intel ifx 2024.2.0 Build 20240602
- 编译选项：`/stand:f18 /warn:all /check:all /traceback`
- oneAPI：`D:\intel\oneAPI`
- Visual Studio 2022：`D:\englishroad\vsstudio\product`

已完成的端到端编译记录为：退出码 0、警告 0、错误 0、链接成功、运行时诊断 0。本文档收尾未重新编译，以避免生成或改动构建产物；结果取自已完成的真实测试记录及现有测试源码。

## 7. 测试验收

| 测试程序 | PASS | FAIL | 状态 |
|---|---:|---:|---|
| `test_saildatabase` | 15 | 0 | PASS |
| `test_sailinterp` | 16 | 0 | PASS |
| `test_sailrelwind` | 16 | 0 | PASS |
| `test_sailangle` | 20 | 0 | PASS |
| `test_sailforce` | 27 | 0 | PASS |
| `test_sailmodule` | 23 | 0 | PASS |

数据库测试计数由 `test_saildatabase.f90` 中编号 1 至 15 的 `ReportResult` 调用核验；其余计数与任务所给已完成测试记录及对应测试程序摘要一致。

## 8. 关键基准结果

### STAR 基准

输入 `V_REL_H_BODY=(0,10,0) m/s`，得到：

```text
alpha_db = 0°
CL = -0.29395224
CD = 0.17179593
q = 61.25 Pa
qS = 918.75 N
F = (-270.0686205, 157.8375106875, 0) N
```

### 非零作用点

对 `R_SAIL_BODY=(2,3,4) m`：

```text
M = (-631.35004275, -1080.274482, 1125.880882875) N·m
```

### 非节点插值

在 `alpha=2.5°`，由 0° 与 5° 节点线性插值得到：

```text
CL = -0.44574592
CD = 0.191519865
F = (-409.529064, 175.95887596875, 0) N
```

### 相对运动

船体沿风向平动 2 m/s 后，原 10 m/s 水平相对风变为 8 m/s。对 `OMEGA_BODY=(0,0,1) rad/s`、`R_SAIL_BODY=(2,0,0) m`：

```text
OMEGA_BODY × R_SAIL_BODY = (0,2,0) m/s
V_REL_H_BODY = (0,8,0) m/s
```

### 重复调用

总接口在一次初始化后连续 1000 次调用，所有被检查输出与首次结果一致。

## 9. 生命周期验收

以下完整流程已由端到端测试覆盖并通过：

```text
未初始化调用（拒绝且输出归零）
→ 初始化
→ 重复初始化（跳过同文件重载）
→ 多次计算
→ 清理
→ 清理后拒绝计算
→ 重新初始化
→ 再次成功计算
→ 最终清理
```

数据库仅在初始化阶段读取；`ComputeSailModuleLoads` 不自动加载数据库。

## 10. 已知限制

- 单风帆；
- 一维 `alpha-CL-CD` 数据库；
- 分段线性插值；
- 仅水平气动力，当前 `Fz=0`；
- 无动态失速；
- 无多帆干扰；
- 无宿主集成；
- 无时域耦合验证；
- 无帆结构变形模型。

## 11. 验收结论

当前 SailModule 已完成独立模块级开发、严格编译、单元测试及端到端测试，可进入宿主适配设计阶段。

下一阶段应确定 `cwbasic`/WALCS、WALCS-LE 或 `planning` 中的接入对象，建立输入映射与 `LOAD_6DOF` 外载荷写入位置，核对坐标、单位和时间步，并先完成单向接入验证，再开展时域风浪耦合。当前结论不表示已经完成任何宿主接入或时域耦合。

本项目为 Walcs 风帆助推船论文代码。

项目根目录固定为：

D:\GitHub\walcs-basic

所有新建、修改和生成的 Walcs、Fortran、Python、说明文档文件，必须保存在该项目目录内。不得在 C 盘用户目录、临时目录或桌面目录中生成项目代码。

Python 解释器固定使用：

D:\Englishroad\python\python.exe

涉及 Walcs 现有源代码时，修改前先分析调用关系。除非任务明确要求，不得大范围重构现有程序。

新增风帆功能优先放在：

D:\GitHub\walcs-basic\SailModule

禁止提交以下文件：

* Visual Studio 编译产物
* Debug、Release、x64、x86 目录
* .obj、.mod、.exe、.pdb 文件
* STAR-CCM+ 大型结果文件
* 临时文件和缓存文件

每次完成修改后，输出：

1. 修改了哪些文件
2. 每个文件修改了什么
3. 是否修改了现有 Walcs 主程序
4. 推荐的 Git commit 信息
5. 编译和验证步骤


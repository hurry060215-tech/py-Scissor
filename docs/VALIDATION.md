# 验证报告：scissor-py 0.1.0.dev0

日期：2026-10-09。范围：独立 Python Scissor APML1-Net 研究原型；本地构建，未发布。

## 结论

Gaussian、Binomial、Cox 的可运行 Python 内核、显式图、ID 对齐、联合分位数归一化、Pearson、固定折 CV、lambda.min 重拟合和 alpha 搜索已实现。运行时只需要 NumPy/SciPy/pandas。266 项 Python 测试全部通过。

固定 R 微型案例验证了官方**已定义且图索引安全的子集**。这不证明所有输入严格兼容。修正模式经过独立数学验收，并明确与官方模式分开。

## 参考基线

- 官方 Scissor SHA：`311560a1f160f665917e7da0fcf289914ee99773`
- 未修改计算语句；固定源码及其 SHA-256 保存在 `reference/upstream/` 与 `reference/manifest.json`
- 实际 R：4.5.0；Rcpp/RcppEigen 编译，Matrix/jsonlite，preprocessCore 1.75.1
- preprocessCore SHA：`5f4522fcbce0820b4983f453e3d56341da62da85`
- Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0、pandas 2.2.3；Linux
- 未验证用户另一个 R 4.5.2 环境

22 份主参考 fixture 的四个 R 生成脚本已从项目目录外重新运行，生成结果与已记录 SHA-256 逐字节一致。额外独立审查保存了 144 次随机路径尝试和 60 次不等大小折 CV；原始 JSON 与生成脚本在 `reference/review/`。

## 数值结果

### 三类显式路径与 CV

同一 24 样本、6 细胞安全非规则图，四个固定 lambda、四折，分别使用官方默认和收紧的公开 thresh。六组案例：

| family | 最大系数绝对误差 | 最大 CV 均值误差 | lambda.min |
|---|---:|---:|---|
| Gaussian | 2.22e-16 | 6.94e-17 | 全匹配 |
| Binomial | 1.78e-15 | 2.22e-16 | 全匹配 |
| Cox | 1.89e-15 | 1.33e-15 | 全匹配 |

这里是保存 fixture 上的实测值，不是对未知数据的精度保证。初始验收门槛仍为归一化/相关 atol=1e-8、rtol=1e-6；稳定系数 atol=1e-5、rtol=1e-4，没有靠放宽阈值取得通过。

### 自动路径、表达到选择

- Gaussian/Binomial 的 20 点自动 lambda 路径和固定折 CV 通过
- Cox 自动路径的实际缩放 λ、包装标签、CV/refit 行为被分别验证
- 三类 raw bulk/sc 表达→ID 重排→联合归一化→相关矩阵→CV→单 λ refit→alpha 停止均与实际 R harness 匹配
- 三类最终标签、alpha、lambda、cutoff 状态匹配；Gaussian/Binomial 有明确 cutoff 未达案例，Cox 达到 cutoff
- 表达 harness 使用显式图，跳过 Seurat；未声称执行完整 `Scissor()` 自动构图工作流
- 实际 preprocessCore 并列值输出验证 midpoint-rank 插值；不是整段 tie 目标均值

### 额外独立差分

- 120 个定义良好的官方随机路径：三个 family、8 个 R 种子、alpha=0/0.4/1，系数在 1e-7 内匹配
- 另外 24 个 p=1 Cox 原包装维度错误被记录，Python official 模式明确拒绝
- 60 个随机不等大小折 CV：系数、均值、SE、lambda.min 在 1e-7 内匹配
- 180 个 corrected 拟合：10 个种子、三类、alpha=0/0.4/1、两个 lambda、相关预测量；均收敛至 KKT≤1e-9，与独立 SciPy 约束优化目标值误差≤1e-7

## 上游问题与处理

1. **Gaussian/Cox 图节点偏移。** R 包装添加截距图节点，但这两类仍只有 p 个系数。实际 R 构造的 p=6 图含 loc=7，按 C++ 的 loc-1 访问会越界。危险内核访问未执行。Gaussian 另有 inactive-set loc 未减一。官方模式保守拒绝可能越界图；Gaussian 最后两个原始节点须孤立，Cox 最后一个须孤立，alpha=1 和 Gaussian p=1 有明确例外。已执行的安全图中，仅改为正确对齐就使系数最大改变 Gaussian 0.139213、Cox 0.076109。
2. **Gaussian 的停机不等于所需最优精度。** 单变量解析解 1.32807217734，原默认得到 1.31971056611。公开 thresh 收紧也不能消除内部固定 1e-5 停机阈值。实际 RSS 递推与返回系数重算 RSS 不一致。official 保留行为和 KKT 诊断；corrected 使用真实坐标二阶项与 KKT。
3. **Cox 自动 λ 标记不一致。** 已保存案例 lambda_max=0.2368093421544794。内核按此缩放，而包装用原比例向量做标签/CV。official 分别记录报告 λ 与真实求解 λ；corrected 使用真实 λ。
4. **单细胞 Cox 包装失败。** p=1 的 R 子集降维导致错误；official 拒绝，corrected 支持。
5. **失败处理更严格。** 候选若全数据拟合、任一折或最终 refit 未收敛则不得成功；所有候选失败抛错。maxit 对 active-set 重启也有效。这些是明确的安全差异。

完整边界与三个目标函数见 `docs/SEMANTICS.md`，实际 R 证据见 `reference/README.md`。未定义的 C++ 内存读取没有被模拟成数值结果。

## 边界测试

覆盖：重复/缺失/错序 ID、二分类正类编码、训练折单一类别、Cox 无事件训练/验证折、并列时间/删失、孤立节点、图非对称/负权重、自环二值化、全零选择、alpha 未达 cutoff、maxit 失败、active-set 重启、None/Generator RNG 跨 alpha 固定折、输入不被修改、AnnData 显式层/图、恒定表达、内存预算、细胞重排、近零诊断、多解预测稳定而系数不稳定。

## 内存与有界耗时

`scripts/validation_metrics.py` 生成可审计 JSON：

- 100,000 节点链图：199,998 非零边，CSR 占 2,799,980 bytes；相同密集矩阵需约 80 GB。测试禁止调用图的 toarray
- Cox 6,000 样本、3 预测量：likelihood score 的 tracemalloc 峰值 722,376 bytes；N×N dense risk matrix 单独需 288 MB。tracemalloc 不是整个进程 RSS
- 60 样本×80 细胞、六个 λ、corrected、thresh=1e-7、单线程 BLAS：Gaussian 0.186 s、Binomial 0.046 s、Cox 0.083 s；均收敛

这些小型数据耗时不外推到真实 atlas。表达归一化和 bulk×cell 相关仍为密集内存。默认 1 GB 参数是预估工作区预算，不是 OS 硬内存限额。未实现 Numba/GPU/并行。

## 安装与制品

最终安装检查及制品校验记录由 `docs/installation_check.json` 提供：已从源码构建 sdist，再从 sdist 构建 wheel；在不继承 system-site-packages 的独立 venv 中安装；离开源码目录、使用 isolated Python、清除 R 环境变量并禁止子进程，三类 smoke 全通过。pip check 无依赖问题。从 sdist 提取测试、关闭源码 pythonpath、针对安装后的 wheel 再跑 266 项测试全部通过。检查了 wheel 与已验收运行时代码的逐文件 SHA-256 一致性。

许可为 GPL-3.0-or-later。源码包含原作者、固定 SHA、原 LICENSE、修改说明、R fixture/生成器、测试、合成例子与本报告。wheel 只包含 Python 内核及 LICENSE/NOTICE，不捆绑 R 或开发编译产物。

## 尚未完成

- 公开 LUAD TP53/生存真实数据验证：官方教程五个数据端点的有界 HEAD 检查均超时，未下载完整数据、未运行或宣称生物学复现
- 完整 Seurat 自动构图结果的等价验证
- Omicverse 集成、技能、真实案例 notebook、远程 PR/PyPI 发布
- 跨 Python/操作系统/BLAS 版本矩阵测试、大规模性能优化、自动唯一性证明

该原型可用于核查实现与设计后续验证；不能据此把具体细胞标签当作已验证的生物学或临床结论。

# scissor-py · Python Scissor 原型

独立、可安装的 Python APML1-Net 内核。运行时只依赖 NumPy、SciPy、pandas；R 仅用于开发期生成参考数据。

当前核心365项测试、真实安装的 OmicVerse 18项接口测试通过；完整 LUAD 的部分高 alpha 收敛对照仍在进行，尚未最终验收。见 [当前验证状态](docs/VALIDATION_STATUS.md)，不要把接口测试通过等同于所有数值对照完成。

发行名暂定 `scissor-py`，导入名为 `scissor_py`。2026-10-09 查询 PyPI：`pyscissor` 已被 netCDF/shapefile 项目占用；`scissor-py` 和 `py-scissor` 的 JSON 接口返回 404。404 不代表名称已保留。本项目源码仓库为 https://github.com/hurry060215-tech/py-Scissor ，尚未发布到 PyPI。

## 先读：两种数值语义

固定参考为 [Scissor 311560a1](https://github.com/sunduanchen/Scissor/tree/311560a1f160f665917e7da0fcf289914ee99773)。实际 R 审计发现上游图索引、Gaussian 停机与 Cox λ 标记问题，因此不能把普通数学重实现直接称为完全兼容。

- `semantics="official"`：默认。复现已定义的官方行为，包括 Gaussian 的原停机逻辑、Cox 自动路径标记/CV 特性。对可能发生越界的 Gaussian/Cox 图明确报错。它通过固定小型 R fixture 对照，不代表任意输入或完整 Scissor 软件栈均兼容。
- `semantics="corrected"`：必须显式选择。图节点与细胞一一对应，Gaussian 使用真实坐标曲率，三类均采用目标函数的 KKT 收敛检查，Cox 使用实际 λ 路径。其结果可能与未修改官方版本不同。

一般真实细胞图建议先阅读 [兼容边界与数值定义](docs/SEMANTICS.md)，再决定要复现历史行为还是运行修正模型。原型不默默在两种模式间切换。

## 安装与测试

```bash
python -m pip install ./scissor_py-0.1.0.dev3-py3-none-any.whl
# 或在源码根目录
python -m pip install -e '.[test]'
python -m pytest -q
python examples/synthetic.py
```

需要 Python ≥3.10。实际验证环境为 Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0、pandas 2.2.3。普通测试读取保存的 R 结果，不需要安装或启动 R。

## 表达矩阵接口

```python
from scissor_py import scissor

result = scissor(
    bulk,                 # DataFrame：基因 × bulk 样本，唯一 ID
    single_cell,          # DataFrame：基因 × 细胞，唯一 ID
    phenotype,            # Series：样本 ID 索引；Cox 用 time/status DataFrame
    graph,                # 稀疏邻接矩阵或按细胞 ID 标记的 DataFrame
    graph_cell_ids=single_cell.columns,  # 数组/稀疏图必须明确给出
    family="binomial",
    positive_label="mutant",
    semantics="corrected",
    alpha=[0.05, 0.2, 0.6],
    foldid=foldid,         # 显式固定；跨语言对照必须使用相同折分
    nlambda=100,
)
print(result.positive_cells, result.negative_cells)
print(result.coefficients, result.metadata)
```

处理顺序：按 bulk 中的顺序取共同基因 → bulk/sc 联合分位数归一化 → bulk 样本 × 细胞 Pearson 相关 → 图去自环、二值化、度归一化 → 各 α 的 CV → `lambda.min` 从零重新拟合 → 按系数严格正/负号标记。最先满足选中比例 `< cutoff` 的 α 停止；若全未达到，返回最后成功结果并明确 `cutoff_reached=False`。

不自动 log 转换、填补缺失、选 HVG、批次校正或构图。输入必须是用户选择的可比较表达层。样本或细胞 ID 不匹配会报错；不会默默删除样本。表达列恒定时 Pearson 未定义，明确拒绝。

`family="binomial"` 要求 0/1 或显式 `positive_label`；Scissor+ 指向编码为 1 的类。Cox 的 `status=1` 表示事件，`0` 表示右删失，使用 Breslow ties；Scissor+ 指向更高 log hazard。所有这些是模型系数解释，不能作为临床判断。

## 低层与 AnnData

```python
from scissor_py import fit_path, cross_validate, scissor_from_anndata

path = fit_path(X, y, raw_graph, family="gaussian", alpha=.4,
                lambdas=[.12, .06, .03], semantics="corrected")
cv = cross_validate(X, y, raw_graph, family="gaussian", alpha=.4,
                     foldid=foldid, semantics="corrected")

result = scissor_from_anndata(
    adata, bulk, phenotype, layer="log1p", graph_key="connectivities",
    family="binomial", positive_label="mutant", semantics="corrected",
    foldid=foldid,
)
```

低层 X 为样本 × 细胞；图为原始非负无向邻接，不是已经归一化的 W。低层接口按位置解释输入；ID 安全对齐由高层完成。AnnData 适配器默认不修改对象，默认读路径不导入 AnnData/Scanpy/Omicverse；必须显式给层和图键，`layer="X"` 表示 adata.X。`writeback=True` 才显式写入真实内存 AnnData 的 obs/uns，冲突默认拒绝；见 [写回、JSON 与原子 H5AD 保存](docs/ANNDATA.md)。

默认 `engine="numpy"` 不导入或需要 Numba。可选安装 `.[acceleration]` 后，corrected Gaussian 可以显式选 `engine="numba"`，禁用 fastmath/并行 JIT；缺依赖或请求其他模型不会静默回退。JIT 编译成本与热运行耗时须分别记录，真实数据加速幅度正在验证。这个开关也可传给 `cross_validate`、`scissor` 和 AnnData 适配器。

## 结果与诊断

结果 dataclass 包含系数、标签、α、λ、截距、每轮 CV、foldid、中间相关矩阵、归一化表达、稀疏图及参考 SHA。路径记录实际求解 λ、迭代次数、收敛状态、KKT 残差。`official` 的 `converged` 表示原版停机条件成功，不保证数学最优；查看 KKT 诊断。

任一候选的全数据拟合或任一折未收敛，该候选不得赢得 CV；最终单 λ refit 也必须成功。全部候选失败时抛 `ConvergenceError`。近零非零系数单独记录，标签仍严格按符号；不靠改变阈值制造匹配。多解模型的单个系数可能不稳定，需要比较预测与目标函数。

## 范围与证据

- 三类固定 R 路径、CV、表达→选择的小型 reference harness，见 [验证报告](docs/VALIDATION.md)
- 原版疑点的实际证据、复现步骤及边界，见 [R 审计](reference/README.md)
- 100,000 节点稀疏图检查；Cox 不构造 N×N risk matrix
- 相似度/联合归一化仍需要密集内存；`max_dense_bytes` 是估算预算，默认 1 GB，不是整个进程 RSS 的硬限制
- 默认 NumPy 原型尚无大规模速度保证；可选 Gaussian Numba 引擎正在全量验证，未加入 GPU/并行
- 官方公开 LUAD TP53/生存完整输入与预处理已核验；模型对照尚在进行，未做生物学有效性结论
- OmicVerse `ov.single.Scissor` 本地薄适配候选与 AnnData 写回正在验证；skill/notebook、完整安装与远程 draft PR 尚未完成
- 仅 APML1-Net，未实现 APML1 的未使用 L0 后修剪、adaptive weights/signs 或其他 penalty

## 许可

GPL-3.0-or-later。见 [NOTICE](NOTICE)、[LICENSE](LICENSE) 和原始源快照。原作者/固定版本/修改记录已保留。论文引用不替代 GPL 义务。

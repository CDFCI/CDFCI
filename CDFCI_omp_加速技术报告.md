# CDFCI_omp 当前稳定版加速技术报告

## 1. 摘要

本文只比较两个版本：

- 原版：`source-full`，即 `CDFCI-source-full.zip` 中的原始实现。
- 当前稳定版：`CDFCI_omp` 工作目录中的当前实现。

当前稳定版没有改变 CDFCI 的数学迭代定义，优化目标集中在 Wave Function Update 阶段的执行方式。最终保留的有效策略是：

1. 在 `z_threshold = 0` 时，使用 callback/visitor 方式边生成 Hamiltonian entry 边更新 `z`，减少 `Column` 中间 vector 的构造和二次遍历。
2. 在 `num_coordinates` 已经提供足够外层并行度时，自适应减少内层 OpenMP task 数量，降低 task 调度开销。
3. 在 `z_threshold > 0` 时，保留接近原版的 column/update_z 压缩路径，因为该路径在压缩模式下更稳定；只保留低风险的小优化，例如 Hamiltonian column 预分配和避免压缩模式下过度 reserve。

需要特别说明：此前尝试过的 presence-filter/Bloom-filter 类优化，以及 sharded/direct/profile 类实验更新后端，已经从当前源码中撤回。它们虽然可以设计成数学结果正确，但实际测试使 Wave Function Update 变慢或收益不稳定，因此不属于当前稳定版成果。

## 2. 数学更新公式

CDFCI 维护两个主要量：

```text
x : 当前波函数系数
z : Hx 的稀疏近似
```

对被选中的 determinant 坐标 `j`，坐标下降先更新：

```text
x_j <- x_j + dx_j
```

因为

```text
z_i = sum_j H_ij x_j
```

所以当 `x_j` 增加 `dx_j` 时，对所有与列 `j` 相连的 `i` 有：

```text
dz_i = H_ij dx_j
z_i <- z_i + dz_i
```

同时需要维护能量相关内积 `x^T z`。若 `x` 先更新，随后 `z` 增加 `dz`，则：

```text
(x + dx)^T (z + dz)
  = x^T z + dx^T z + (x + dx)^T dz
```

代码中对应关系为：

- `update_x(det, dx)` 更新 `x_j`，同时维护 `x^T x` 和 `dx^T z_old`。
- `update_z(...)` 或 `update_z_entry(...)` 执行 `z_i += H_ij dx_j`，并累积 `(x + dx)^T dz`。
- `reinsert_new_z(det, new_z)` 对被选中的 determinant 重新写入精确的 `new_z_j`，保证被选坐标的 `z_j` 与当前 `x` 一致。

当前稳定版只改变遍历和并行组织方式，不改变上述数学关系。

## 3. `z_threshold` 语义

对每个候选更新：

```text
dz = H_ij dx_j
norm = || scale * dz ||_infinity
```

阈值规则为：

```text
if norm > z_threshold:
    插入或更新 z_i
else:
    只在 z_i 已经存在时更新，不新增很小的 z_i
```

因此 `z_threshold > 0` 时，压缩只减少新增的 `z` 分量，不代表可以跳过 Hamiltonian column 的生成。即使一个小 `dz` 不会新增 `z_i`，如果 `z_i` 已经存在，仍然必须更新它，否则数学结果会改变。

这也是当前稳定版在 `z_threshold > 0` 选择保守 column/update_z 路径的原因：安全剪枝空间很小，很多看似可以跳过的小项仍需要先判断目标 determinant 是否已经存在。

## 4. 原版实现

原版 Wave Function Update 的主要流程如下：

```text
function update_coordinate_original(det_picked):
    sub_xz.clear()

    for (det, dx) in det_picked:
        update_x(det, dx)

    parallel for (det, dx) in det_picked:
        new_z = 0

        taskgroup:
            spawn serial task:
                column = get_column_serial_part(det)
                new_z_part = update_z(column, dx, sub_xz_part)
                critical merge(sub_xz_part, new_z_part)

            for tid in 0 .. nelec / 2 - 1:
                spawn double-excitation task:
                    column = get_column_parallel_part(det, tid, nelec / 2)
                    new_z_part = update_z(column, dx, sub_xz_part)
                    critical merge(sub_xz_part, new_z_part)

        det_picked[det].z = new_z

    for (det, new_z) in det_picked:
        reinsert_new_z(det, new_z)

    update global xz and size_z
```

原版的主要开销来自：

1. 每个 task 先构造 `Column = vector<pair<determinant, H_ij>>`，随后 `update_z` 再遍历这个 vector。
2. `Column` 的内存分配、`push_back`、扩容和缓存访问成本较高。
3. 在 `num_coordinates = 8` 且 `OMP_NUM_THREADS = 8` 时，外层坐标已经可以占满线程，但原版仍按 `nelec / 2` 生成大量内层 task。
4. 每个 task 最后进入 `#pragma omp critical` 合并局部 `sub_xz` 和 `new_z`，task 越多，调度和合并成本越高。

## 5. 当前稳定版改动

### 5.1 `z_threshold = 0` 的 callback 快路径

当前稳定版增加了 `ColumnSink` 和 fast visitor 接口：

```text
for_each_column_serial_part_fast(det, sink)
for_each_column_parallel_part_fast(det, tid, inner_tasks, sink)
```

在 `z_threshold = 0` 时，不再先生成完整 `Column` vector，而是：

```text
Hamiltonian entry 生成出来
    -> 立即 callback
    -> update_z_entry(...)
```

伪代码：

```text
function visit_legacy_work_unit(det, dx, part):
    delta_xz = 0
    sub_local.reserve(estimated_size)
    new_z_part = 0

    sink(target, h_value):
        x_current = update_z_entry(target, h_value, dx,
                                   sub_local, delta_xz)
        new_z_part += x_current * h_value

    if part is serial:
        for_each_column_serial_part_fast(det, sink)
    else:
        for_each_column_parallel_part_fast(det, tid, inner_tasks, sink)

    sub_local.update_xz(delta_xz)
```

相对原版：

```text
原版:   generate entries -> Column vector -> update_z(Column)
当前版: generate entry   -> callback      -> update_z_entry
```

收益是减少中间 vector 的构造、写入和二次遍历。

### 5.2 自适应内层 task 数量

原版内层 task 数通常固定为：

```text
inner_tasks = nelec / 2
```

对 N2/cc-pVDZ 测试：

```text
nelec = 14
inner_tasks = 7
```

若每轮选 `num_coordinates = 8` 个坐标，则原版每轮近似产生：

```text
8 * (1 + 7) = 64 个 task
```

当前稳定版在 `z_threshold = 0` 时根据外层并行度自适应：

```text
max_threads = omp_get_max_threads()
active_coords = det_picked.size()
target_inner_tasks = ceil(max_threads / active_coords)
inner_tasks = min(nelec / 2, target_inner_tasks)
```

当：

```text
max_threads = 8
active_coords = 8
```

则：

```text
inner_tasks = 1
```

每轮 task 数近似降为：

```text
8 * (1 + 1) = 16 个 task
```

task 数减少约 75%，因此 OpenMP 调度和 critical 合并次数同步下降。

### 5.3 `z_threshold > 0` 的稳定压缩路径

当前稳定版在 `z_threshold > 0` 时不使用 callback 快路径，而使用接近原版的 column/update_z 路径：

```text
generate_column_parallel_task
    -> get_column_serial_part / get_column_parallel_part
    -> update_z(column, dx, sub_xz, z_threshold)
```

选择这个路径的原因是实际测试表明，在压缩模式中：

1. 许多 `dz` 不会新增到 `z`，hash insert/append 的成本下降。
2. 此时 Hamiltonian column 生成、阈值判断和已有 key 查询成为主要固定成本。
3. 原版式 column vector 路径虽然有中间存储，但访问模式更连续；callback 直接插入哈希表反而可能破坏局部性。
4. 任何试图跳过小 `dz` 的方案，都必须先保证目标 determinant 不存在，否则会改变数学结果。

因此当前稳定版对 `z_threshold > 0` 的定位是“结果正确优先、性能接近原版、保留低风险小优化”，不再追求激进剪枝。

### 5.4 Hamiltonian column 预分配

当前稳定版在 Hamiltonian column 构造处增加了 `reserve_size` 估计，例如：

```text
reserve_size = 1
reserve_size += single_excitation[i].size()
reserve_size += double_excitation[index(i, j)].size()
result.reserve(reserve_size)
```

这减少 `std::vector` 扩容和搬移。该优化不改变生成的 determinant 和 Hamiltonian value，只降低内存分配开销。

### 5.5 压缩模式避免过度 reserve

`update_z(column, dx, sub_xz, z_threshold)` 中，当前稳定版只在无压缩时按 `column.size()` 预留 `sub_xz`：

```text
if (z_threshold <= 0.0)
    sub_xz.reserve(column.size())
```

原因是 `z_threshold > 0` 时，最终真正进入 `sub_xz` 的条目可能远小于 column size。强行按 column size 预留会造成过多内存申请和缓存压力。

## 6. 正确性说明

当前稳定版保持正确性的关键点：

1. 所有实际 `z` 更新仍通过 `update_z(...)` 或 `update_z_entry(...)` 完成。
2. `z_threshold` 判断逻辑没有改变：大于阈值插入或更新，小于阈值只更新已有项。
3. 对被选中的 determinant，仍然在每轮末尾执行 `reinsert_new_z(det, new_z)`。
4. `x^T x`、`x^T z`、`size_x`、`size_z` 的维护逻辑没有改变数学定义。
5. callback 快路径只改变 entry 的传递方式，不改变 entry 内容和更新公式。

因此当前稳定版属于实现层重构，不属于算法近似修改。

## 7. 性能结论

### 7.1 `z_threshold = 0`

在无压缩模式下，当前稳定版的主要收益来自：

```text
减少 Column 中间 vector
+ 减少 OpenMP task 数
+ 减少 critical 合并次数
```

用户前序测试表明，当前稳定版相对 `source-full` 在 N2/cc-pVDZ、8 OpenMP 线程、`num_coordinates = 8`、100000 iterations 条件下可以获得约 10% 量级的加速。

这里的加速属于常数项优化，不改变算法复杂度。Wave Function Update 仍然需要遍历 Hamiltonian 连接并访问哈希表，因此不应预期出现数倍加速。

### 7.2 `z_threshold = 1e-5`

在典型压缩模式 `z_threshold = 1e-5` 下，当前稳定版已经回到保守 column/update_z 路径。该模式下压缩会显著减少 `|z|_0`，但不会同比例减少运行时间，原因是：

```text
必须仍然生成 Hamiltonian candidates
必须仍然判断 dz 是否过阈值
小 dz 若目标已存在仍必须 update_fn
```

## 8. 推荐输入

当前稳定版推荐使用：

```json
{
    "hamiltonian": {
        "type": "molecule",
        "molecule": {
            "fcidump_path": "/home/zyzhang/sourcecode/CDFCI/data/n2_ccpvdz_psi4.FCIDUMP",
            "verbose": 1
        }
    },
    "solver": {
        "type": "cdfci",
        "cdfci": {
            "verbose": 1,
            "num_iterations": 100000,
            "report_interval": 10000,
            "num_coordinates": 8,
            "coordinate_pick": "block_gcd_grad",
            "coordinate_update": "eig",
            "wavefunction_update_backend": "legacy_cuckoo",
            "z_threshold": 1e-5,
            "z_threshold_search": false
        }
    },
    "max_memory": 25.0,
    "max_load_factor": 0.79,
    "verbose": 1
}
```

其中最关键的是：

```text
wavefunction_update_backend = legacy_cuckoo
num_coordinates = 8
z_threshold = 1e-5
```

若要验证无压缩加速，将 `z_threshold` 改为 `0.0` 即可。

## 9. 后续优化空间判断

当前版本已经把低风险优化集中在可验证的位置。继续优化有两类可能：

1. 保持数学结果严格一致：空间较小，主要是继续减少 task/critical 常数开销，预期收益有限。
2. 改变压缩或剪枝策略：可能获得更大加速，但会改变数值路径，必须重新定义“数学结果正确”的容忍标准。

不建议继续投入的方向：

- presence-filter 或 Bloom-filter 负查找过滤：实测成本大于收益，已撤回。
- sharded/direct/profile 类实验更新后端：额外缓冲、锁、合并和统计成本抵消压缩收益，已从当前代码中删除。
- 在 `z_threshold > 0` 下直接按积分大小截断 Hamiltonian entry：会漏掉已有 `z_i` 的小更新，严格语义下不正确。
- 仅增大 `max_memory`：容量足够后，对 Wave Function Update 的主时间帮助有限。

## 10. 汇报用总结

当前稳定版的优化没有改变 CDFCI 的数学公式，而是重构了 Wave Function Update 的执行路径。在无压缩模式下，它将原版“先构造 Hamiltonian column vector，再遍历更新”的流程改为“边生成 Hamiltonian entry，边 callback 更新”，同时根据外层坐标并行度减少内层 OpenMP task 数，从而降低内存分配、重复遍历、task 调度和 critical 合并成本。在压缩模式 `z_threshold > 0` 下，当前版本选择保守的原版式 column/update_z 路径，以保证结果正确和性能稳定；因此该模式下预期接近原版，不承诺稳定 10% 加速。

# CDFCI

[![CDFCI\_Test](https://github.com/CDFCI/CDFCI/actions/workflows/cdfci_test.yml/badge.svg?branch=master)](https://github.com/CDFCI/CDFCI/actions/workflows/cdfci_test.yml)

CDFCI（坐标下降全组态相互作用方法）是一个用于量子化学电子结构基态计算的现代高效求解器，基于 C++17 实现，支持大规模稀疏波函数与自适应压缩策略。

English version: [README.md](../README.md)

---

## ✨ 概览

CDFCI 将全组态相互作用（FCI）本征值问题重写为一个大规模无约束非凸优化问题，并使用带确定性压缩的自适应坐标下降法进行求解。

核心特点：

* 基于 Slater 行列式的坐标方向更新
* 通过 `z_threshold` 进行自适应压缩
* 确定性且面向内存的波函数存储策略
* 可选支持 OpenMP（多线程）

---

## 📄 引用文献

本软件包的主要论文为：

* **CDFCI 软件论文**：
  Y. Zhang, Z. Wang, J. Lu, Y. Li, [*arXiv:2605.04483*, 2026](https://arxiv.org/abs/2605.04483)

具体算法相关论文包括：

* **CDFCI 算法提出**：
  Z. Wang, Y. Li, J. Lu, [*JCTC*, 15(6), 2019](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138)

* **最优轨道选择（OptOrbFCI）扩展**：
  Y. Li, J. Lu, [*JCTC*, 16(10), 2020](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00613)

* **低激发态版本（xCDFCI）**
  Z. Wang, Z. Zhang, J. Lu, Y.Li, [*JCTC*, 19(21), 2023](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.3c00452)

* **多坐标并行版本（mCDFCI）**：
  Y. Zhang, W. Gao, Y. Li, [*JCTC*, 21(5), 2025](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01530)

以及收敛性保障：

* **理论分析与收敛性证明**：
  Y. Li, J. Lu, Z. Wang, [*SIAM J. Sci. Comput.*, 41(4), 2019](https://doi.org/10.1137/18M1202505)

如果您使用了本软件，请引用相关论文。

---

## 🚀 快速开始

如需 Python 工作流或高级脚本化使用，请参见 [Python 接口用户手册](../PYTHON_INTERFACE_USER_MANUAL.md)。

### 克隆并编译

```bash
git clone --recursive https://github.com/CDFCI/CDFCI.git
cd CDFCI
make
```

将在 `bin` 目录下生成可执行文件：

* `cdfci`：单线程版本
* `cdfci_omp`：启用 OpenMP 的版本
* `xcdfci`：用于激发态
* `optorbfci`：用于轨道旋转与压缩
* `cdfci_tools`：辅助工具

或者，您也可以使用 `make cdfci`、`make optorbfci` 或 `make xcdfci` 生成指定程序。

发布相关常用目标：

* `make check`：测试门禁（`make test` 的别名）
* `make install PREFIX=/path`：安装程序到 `${PREFIX}/bin`（默认 `/usr/local/bin`）
* `make uninstall PREFIX=/path`：卸载已安装程序
* `make release-check`：执行 clean build + regression tests + examples 的发布前检查

### 编译要求

* 支持 C++17 的编译器（GCC >= 9、Clang >= 10 或 ICC）
* 依赖库：[Eigen3](https://eigen.tuxfamily.org)（通过递归克隆自动获取）
* 可选依赖：OpenMP

---

## 📥 输入格式（JSON）

求解器使用 JSON 输入文件。示例：

```json
{
    "hamiltonian": {
        "type": "molecule",
        "molecule": {
            "fcidump_path": "/path/to/FCIDUMP",
            "threshold": 0.0
        }
    },
    "solver":{
        "type": "cdfci",
        "cdfci": {
            "num_iterations": 150000,
            "report_interval": 10000,
            "z_threshold": 0,
            "z_threshold_search": false
        }
    },
    "max_memory": 0.005
}
```

更多示例（包括完整的输入输出说明）见 [examples 文档](../examples/README.md)。

快速运行：

```bash
bin/cdfci examples/demo_input_cdfci.json
```

---

## 🧪 测试

```bash
make test
```

示例体系与期望能量见 `regression_tests`。

---

## 📊 复现论文结果

`examples/` 与 `papers/` 目录包含复现实验结果所需的全部脚本与输入文件。生成 FCIDUMP 文件需要 Python 与 [PySCF](https://github.com/pyscf/pyscf)。

---

## 👥 开发者

* [张悦嘉](https://ninotreve.github.io/)（复旦大学数学科学学院）
* 王喆（杜克大学数学系）
* [鲁剑锋](https://services.math.duke.edu/~jianfeng/)（杜克大学数学系、物理系与化学系）
* [李颖洲](http://yingzhouli.com/)（复旦大学数学科学学院、上海市应用数学和力学研究所）

---

## 🪙 致谢

本工作部分受以下项目资助：国家自然科学基金项目 12271109 和 12526211；美国国家科学基金会（NSF）项目 DMS-1454939 和 DMS-2012286；上海市基础研究特区计划-复旦大学项目 21TQ1400100（22TQ017）；青年教师科研创新能力支持项目 SRICSPYF-ZY2025159；以及复旦大学学敏高等研究院。

---

## ⚖️ 开源协议

本项目遵循 BSD 3-Clause License，详见根目录下的 [LICENSE](../LICENSE) 文件。


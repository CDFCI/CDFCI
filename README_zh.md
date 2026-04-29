# CDFCI

[![CDFCI\_Test](https://github.com/CDFCI/CDFCI/actions/workflows/cdfci_test.yml/badge.svg?branch=master)](https://github.com/CDFCI/CDFCI/actions/workflows/cdfci_test.yml)

CDFCI（坐标下降全组态相互作用方法）是一个用于量子化学电子结构基态计算的现代高效求解器，基于 C++17 实现，支持大规模稀疏波函数压缩与高性能迭代求解。

English version: [README.md](../README.md)

---

## ✨ 算法简介

CDFCI 将全组态相互作用（FCI）本征值问题重写为一个无约束非凸优化问题，并使用带压缩策略的自适应坐标下降法进行求解。

核心特点：

* 基于 Slater 行列式的坐标方向更新
* 可控的压缩精度 (`z_threshold`)
* 面向内存的截断策略
* 支持 OpenMP 的并行实现（可选）

---

## 📄 引用文献

CDFCI 实现了以下论文中描述的算法：

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

如果您使用了本软件包或者CDFCI算法，请引用相关文章。

---

## 🚀 快速开始

### 克隆并编译

```bash
git clone --recursive https://github.com/CDFCI/CDFCI.git
cd CDFCI
make
```

将在 `bin` 目录下生成可执行文件：

* `cdfci`：单线程主程序
* `cdfci_omp`：开启 OpenMP 的多线程版本
* `xcdfci`：计算激发态
* `optorbfci`：轨道旋转与压缩
* `cdfci_tools`：工具程序

或者，您也可以用 `make cdfci`，`make optorbfci` 或 `make xcdfci` 生成相应的程序。

发布相关常用目标：

* `make check`：测试门禁（`make test` 的别名）
* `make install PREFIX=/path`：安装程序到 `${PREFIX}/bin`（默认 `/usr/local/bin`）
* `make uninstall PREFIX=/path`：卸载已安装程序
* `make release-check`：执行 clean + build + test + examples 的发布前检查

### 编译要求

* 支持 C++17 的编译器（推荐 GCC ≥ 9 或 Clang ≥ 10）
* 依赖库：Eigen3（通过`git clone --recursive`自动下载）
* 可选依赖：OpenMP

---

## 📥 输入格式（JSON）

CDFCI 使用 JSON 配置文件控制计算。

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

更多示例文件，以及完整的使用说明，见 `examples/` 文件夹。

---

## 🧪 测试用例

```bash
make test
```

将运行一组小分子体系下的基准测试，输出能量并进行误差检查。

---

## 📊 复现论文结果

`examples/` 与 `papers/` 目录中包含所有用于复现实验的输入文件与脚本。需安装 [PySCF](https://github.com/pyscf/pyscf) 来生成 FCIDUMP 文件。

---

## 👥 作者信息

* 王喆（杜克大学数学系）
* [鲁剑锋](https://services.math.duke.edu/~jianfeng/)（杜克大学数学系）
* [李颖洲](http://yingzhouli.com/)（复旦大学数学科学学院）
* [张悦嘉](https://ninotreve.github.io/)（复旦大学数学科学学院）

---

## 🪙 致谢

本工作部分受以下项目资助：美国国家科学基金会（NSF）项目 DMS-1454939 和 DMS-2012286；国家自然科学基金项目 12271109 和 12526211；上海市基础研究特区计划-复旦大学项目 21TQ1400100（22TQ017）；以及青年教师科研创新能力支持项目 SRICSPYF-ZY2025159。

---

## ⚖️ 开源协议

本项目遵循 BSD 3-Clause License，详见根目录下的 [LICENSE](./LICENSE) 文件。


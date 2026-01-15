# DMol3 - AtomGrid Toolkit

<details>
<summary><strong>🇬🇧 English README (Click to expand)</strong></summary>

<br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17089110.svg  )](https://doi.org/10.5281/zenodo.17089110) 

**A powerful post-processing tool for processing, analyzing, and converting 3D grid data from computational chemistry, with full support for formats like DMol³ .grd, VASP CHGCAR, and Gaussian .cube.**

The development context and application of this toolkit are detailed in our research paper:
> Wang, X., Zhang, R., et al. (2025). Theoretical Investigation of Pd and Sc Decorated PtS2 Monolayer with Enhanced and Strain-Tunable Sensing Properties for Toxic Gases from LIB Thermal Runaway. DOI: [].

`AtomGrid Toolkit` is designed to solve the interoperability challenges of volumetric data faced by users in the field of computational chemistry, particularly those using DMol³. By providing robust format conversion and quantitative analysis capabilities, it bridges the data gap between DMol³ and mainstream analysis and visualization software such as Bader, Critic2, VESTA, and VMD.

---

## Key Features

*   **Arithmetic on Grid Data**: Supports arithmetic operations (addition and subtraction) on volumetric data files, including `.grd`, VASP `CHGCAR`, and Gaussian `.cube` formats.
*   **1D Profile Analysis**: Performs a variety of quantitative analyses along any lattice axis (x, y, z), including the calculation of plane-averaged charge density (Δρ(z)) and charge displacement curves (ΔQ(z)).
*   **Universal Format Conversion**: Provides seamless and accurate conversion between the non-standard DMol³ `.grd` format, VASP `CHGCAR`, and Gaussian `.cube`.
*   **Structural Information Integration**: Allows for the integration of atomic structural information from an external `.cif` file into the output `CHGCAR` or `.cube` files during conversion.
*   **Geometric Robustness**: Offers full support for non-orthogonal cells, with correct handling of unit conversions (Å/Bohr) and coordinate system definitions across different formats.

---

## Installation

This tool requires Python 3 and the NumPy library.

1.  **Ensure you have Python 3 and NumPy installed**:
    ```bash
    pip install numpy
    ```

</details>

<details>
<summary><strong>🇨🇳 中文说明 (点击展开)</strong></summary>

<br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17089110.svg)](https://doi.org/10.5281/zenodo.17089110) 

**一个功能强大的后处理工具，用于处理、分析和转换计算化学中的三维网格数据，全面支持 DMol³ .grd、VASP CHGCAR 和 Gaussian .cube 等多种格式。**

该工具的开发背景和应用已在我们的研究论文中详细介绍：
> Wang, X., Zhang, R., et al. (2025). Theoretical Investigation of Pd and Sc Decorated PtS2 Monolayer with Enhanced and Strain-Tunable Sensing Properties for Toxic Gases from LIB Thermal Runaway. DOI: [].

`AtomGrid Toolkit` 旨在解决计算化学领域，特别是 DMol³ 用户面临的体数据（volumetric data）互操作性挑战。它通过提供稳健的格式转换和定量分析功能，打通了 DMol³ 与 Bader、Critic2、VESTA、VMD 等主流分析和可视化软件之间的数据壁垒。

---

## 主要功能

*   **网格数据计算**：支持对体积数据文件`.grd`、VASP `CHGCAR` 和 Gaussian `.cube` 格式进行加减运算。
*   **一维曲线分析**：可沿任意晶格轴（x, y, z）进行包含平面平均电荷密度 Δρ(z) 和电荷位移曲线 ΔQ(z)在内等多种定量分析。
*   **通用格式转换**：在 DMol³ 输出的非标准格式的`.grd`、VASP `CHGCAR` 和 Gaussian `.cube` 格式之间进行无缝、精确的相互转换。
*   **结构信息整合**：在转换过程中，可以从外部 `.cif` 文件读入原子结构信息，并将其整合到输出的 `CHGCAR` 或 `.cube` 文件中。
*   **几何鲁棒性**：完全支持非正交晶胞，并能正确处理不同格式间的单位换算（Å/Bohr）和坐标系定义。

---

## 安装

本工具依赖于 Python 3 和 NumPy 库。

1.  **确保已安装 Python 3 和 NumPy**：
    ```bash
    pip install numpy
    ```

</details>

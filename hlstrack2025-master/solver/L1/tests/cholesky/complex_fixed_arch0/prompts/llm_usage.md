# 大模型辅助使用记录

## 基本信息

- **模型名称**：GPT-5 Thinking

- **提供方 / 访问方式**：https://chatgpt.com

- **使用日期**：2025-10-31
- **项目名称**：AES-128 GMAC L1 算子优化

---

## 使用场景 1

### 主要用途

去除原本代码中的多余宏定义，收敛成仅保留单内核 + 流接口的最小骨架

### 完整 Prompt 内容

目标：把多实现（Basic/Alt/Alt2/Top）收敛为一个列主序 Cholesky 内核，保留 stream I/O。
行动：
1) 删除/停用 choleskyTop/choleskyBasic/choleskyAlt/choleskyAlt2，仅保留：
   template <bool LOWER_TRIANGULAR, int DIM, typename Tin, typename Tout, typename Traits>
   int cholesky(hls::stream<Tin>& inA, hls::stream<Tout>& outL);
2) 在该函数内：先把 DIM×DIM 从 inA 读入本地 A[DIM][DIM]，计算后再从 outL 写回。
3) 给出可编译的最小骨架（含必要 #include 与头卫）。后续步骤在此骨架上增量修改。
要求：只输出完整头文件代码，不要解释。


### 模型输出摘要
产出仅包含一个顶层 cholesky 函数签名与 stream I/O 的可编译头文件骨架。

包含必要 #include、头文件卫士；不含任何多架构实现与解释文字。

### 人工审核与采纳情况
将改完后的代码对比原代码中的接口是否改动，经过验证后直接运用完整输出代码。
---

## 使用场景 2

### 主要用途
建立 detail 工具集（复数/零值/饱和/类型抽取）

### 完整 Prompt 内容
目标：统一小工具，替代分散的 hls::x_*。
行动：新增 namespace detail，全部内联（#pragma HLS INLINE）实现：
- template<class T> struct scalar_of;  // 从 hls::x_complex<T> 或标量抽取标量类型
- make_cpx(real, imag), zero<T>(), zero_c<T>(), conjx(z), mul_conj(a,b) // 返回 a*conj(b)
- template<class Out, class In> Out saturate_cast(In) // 工作位宽安全收敛到输出类型
要求：补齐必要 using/typedef，保持可编译；只输出完整头文件最新版本。


### 模型输出摘要
头文件新增 detail 命名空间和可复用工具函数，可独立编译。

不改变现有接口；后续计算可直接使用 detail::*。
### 人工审核与采纳情况
验证 mul_conj 实虚展开是否正确；saturate_cast 是否饱和/截断一致

---

## 使用场景 3

### 主要用途
重塑 Traits：工作位宽与行为参数

### 完整 Prompt 内容
目标：用统一 Traits 管理位宽与展开策略。
行动：定义
template<typename Tin, typename Tout> struct CholeskyTraits {
  using InScalar  = typename detail::scalar_of<Tin>::type;
  using OutScalar = typename detail::scalar_of<Tout>::type;
  // 依据 InScalar 推导工作位宽（可按位宽/整数位数估计）
  static constexpr int ACC_W = (/*bitwidth(InScalar)*/<24)?32:(/*bitwidth(InScalar)*/+8);
  static constexpr int ACC_I = (/*intbits(InScalar)*/<12)?16:(/*intbits(InScalar)*/+4);
  using WorkScalar  = ap_fixed<ACC_W, ACC_I>;
  using WorkComplex = hls::x_complex<WorkScalar>;
  using IComplex = hls::x_complex<InScalar>;
  using OComplex = hls::x_complex<OutScalar>;
  static constexpr int  INNER_UNROLL = 4;
  static constexpr int  DIAG_UNROLL  = 3;
  static constexpr bool USE_DATAFLOW = false;
};
// 兼容旧模板签名
template<bool LOWER,int DIM,typename Tin,typename Tout>
struct choleskyTraits : CholeskyTraits<Tin,Tout> {};
要求：保持可编译；只输出完整头文件最新版本。


### 模型输出摘要
新增 CholeskyTraits 与兼容层；提供工作/输入/输出类型别名与展开参数。

仍保持单内核接口与可编译状态。

### 人工审核与采纳情况
确认 ACC_W/ACC_I 取值是否满足你的定点策略；需要可调宏就备注

## 使用场景 4

### 主要用途
本地缓存与数组分解（列主序并行度准备）

### 完整 Prompt 内容
目标：准备计算存储与并行度。
行动：在 cholesky(...) 内添加：
- IComplex A[DIM][DIM]; WorkComplex Llow[DIM][DIM];
- #pragma HLS ARRAY_PARTITION variable=A    complete dim=2
- #pragma HLS ARRAY_PARTITION variable=Llow complete dim=2
- 读入/写出环均加 #pragma HLS PIPELINE II=1
要求：不改变功能；只输出完整头文件最新版本。


### 模型输出摘要
本地数组与 complete 分解 生效；读/写环具备 II=1 基础。

仍未引入核心数学，保持可编译。

### 人工审核与采纳情况
资源估计是否可接受；需要时改为 block/cyclic 分解

## 使用场景 5

### 主要用途
对角计算采用 rsqrt 稳定路径

### 完整 Prompt 内容
目标：按 8k3-7 风格用倒根路径构造对角。
行动：在主列循环 for (int j=0; j<DIM; ++j) 内：
- diag_r = A[j][j].real() - Σ_{p<j} |Llow[j][p]|^2；若 <0 截到 0 得 df
- inv_f  = hls::rsqrt((float)df)；diag_s = df * inv_f
- Llow[j][j] = WorkComplex(diag_s, 0)
要求：所有计算走 WorkScalar/WorkComplex；只输出完整头文件最新版本。

### 模型输出摘要
对角元素由 df * rsqrt(df) 生成；负值保护生效。

不改变外部接口；可编译并通过综合。

### 人工审核与采纳情况
检查定点到 float 的转换是否影响精度；必要时改用定点 rsqrt 近似
## 使用场景 6

### 主要用途
非对角更新与流水展开

### 完整 Prompt 内容
目标：实现行更新与性能展开。
行动：对 i=j+1..DIM-1：
- acc = A[i][j] - Σ_{p<j} detail::mul_conj(Llow[i][p], Llow[j][p]);
- Llow[i][j] = acc * inv_f;  // 把 inv_f 乘到实/虚
- 内环添加：#pragma HLS PIPELINE II=1
           #pragma HLS UNROLL factor=2
           #pragma HLS DEPENDENCE variable=Llow inter false
要求：更新完先缓存在 Llow；此步不写流。只输出完整头文件最新版本。


### 模型输出摘要
行更新计算正确（使用 mul_conj）；关键循环达到 II=1，部分展开 factor=2。

依赖消除指示添加，便于综合器调度。
### 人工审核与采纳情况
用小维度矩阵比对数值；观察 II、时序与 LUT/BRAM 使用

## 使用场景 7

### 主要用途
写出阶段选择上下三角；取消额外清零循环

### 完整 Prompt 内容
目标：在写出时决定上下三角，移除独立清零。
行动：新增统一写出环：
- 若 LOWER_TRIANGULAR：输出主对角+下三角，其余用 detail::zero_c<OutScalar>()
- 否则：输出 U = Lᴴ（共轭转置）
- 写出前通过 detail::saturate_cast<OutScalar> 收敛类型
要求：仅在写出阶段做零化；只输出完整头文件最新版本。


### 模型输出摘要
产出单一写出路径，无独立“清零上/下三角”循环。

输出类型安全收敛到 Tout。
### 人工审核与采纳情况
确认上/下三角选择正确；共轭转置符号与索引无误

## 使用场景 8

### 主要用途
统一用 detail 工具替换 hls::x_* 残留

### 完整 Prompt 内容
目标：风格一致、可复用。
行动：将残留的 hls::x_conj/x_real 等调用替换为：
- detail::conjx(z)、detail::mul_conj(a,b) 或显式实/虚访问
- 所有零值使用 detail::zero/zero_c
要求：功能等价；只输出完整头文件最新版本。



### 模型输出摘要
代码不再直接调用 hls::x_* 工具函数；全部走 detail::*。

统一风格便于后续复用与测试。

输出类型安全收敛到 Tout。
### 人工审核与采纳情况
代码搜索确认无遗留；编译与单元测试通过


## 使用场景 9

### 主要用途
工程化抛光（包含/宏/可读性/综合指示）
### 完整 Prompt 内容
目标：整理工程细节并定版。
行动：
- 头部仅保留需要的 #include <hls_stream> <hls_math> <ap_fixed> <hls_x_complex>
- 顶层函数使用 #pragma HLS INLINE off
- 可选暴露 CHOLESKY_INNER_UNROLL/CHOLESKY_DIAG_UNROLL 宏，但以 Traits 为准
- 删除未用类型/函数/宏与历史注释；保证命名与 8k3-7 风格一致
要求：输出最终版完整头文件代码。



### 模型输出摘要
最终头文件结构清爽、仅含必要依赖；宏与 Traits 关系明确。

可直接进入综合/实现阶段。
### 人工审核与采纳情况
运行 C/RTL cosim 或功能测试；记录资源与时序


## 使用场景 10

### 主要用途

最终验收与对比输出

### 完整 Prompt 内容
目标：一次性核验“等价风格”，并给出最终成品。
行动：
- 在回复顶部给出勾选清单（接口、II=1、无清零循环、rsqrt 路径、Work* 计算、detail::* 全替换）
- 随后给出完整最终头文件（再次完整贴出），无需解释
- 最后附“最小用例”片段：DIM=4、ap_fixed<16,4> 的示例调用（仅作为注释）
要求：除清单、代码与注释外不输出其他文字。



### 模型输出摘要
先给核验清单，再给最终版完整代码，最后附最小示例注释。

便于人工复核与归档对比。
### 人工审核与采纳情况
对照 checklist 勾验；如有偏差，回滚到对应使用场景修正后再执行本场景



## 总结

### 整体贡献度评估
大模型在本项目中的总体贡献占比：约 65%

        方案规划与分步提示词设计：≈ 35%

        架构收敛（单内核 + 流接口）与 Traits 重塑示例：≈ 20%

        性能优化策略与验收清单（II=1、UNROLL、数组分解等）：≈ 10%

主要帮助领域：

        代码结构重构：从多架构合并为单内核；引入 detail 工具集；写出阶段选择上下三角。

        性能优化建议：ARRAY_PARTITION、PIPELINE II=1、UNROLL、依赖消除、定点位宽规划、rsqrt 稳定路径。

        文档化与复用：10 个“使用场景卡片” + 每段“模型输出摘要”，便于审核、追踪与复用。

人工介入与修正比例：约 35%

        工程环境适配（依赖、仿真脚手架、宏开关）；

        功能/数值验证与边界处理（非 SPD、定点溢出）；

        资源与时序权衡（LUT/BRAM/频率）及参数微调。

### 学习收获
1.Traits 驱动的类型与展开策略：在不改对外接口的前提下，灵活切换内部位宽与行为参数。

2.rsqrt 稳定路径与负值保护：diag = max(a_ii - Σ|L|^2, 0) * rsqrt(...)，提升数值稳健性。

3.WorkScalar 位宽规划：输入位宽 + 安全余量（如 W+8/I+4），最终用 saturate_cast 收敛到输出类型。

4.并行化要点：ARRAY_PARTITION complete + PIPELINE II=1 + 适度 UNROLL，并用 DEPENDENCE inter false 消除调度障碍。

5.数据流与存储策略：本地缓存计算，写出阶段再做上/下三角选择，避免额外清零循环。

6.复数计算规范化：以 detail::mul_conj/conjx/zero_c 统一实现，降低对外部库工具函数的耦合。

7.验证工作流：小维度数值对拍 → C/RTL cosim → 资源/时序回归；步骤化、可回滚。

8.提示工程方法：将大型重构拆解为 10 个最小可验证步骤，为每步加入“模型输出摘要”和人工审核位，提升落地效率。


---

## 附注

- 请确保填写真实、完整的使用记录
- 如未使用大模型辅助，请在此文件中注明"本项目未使用大模型辅助"
- 评审方将参考此记录了解项目的独立性与创新性

#ifndef XF_SOLVER_CHOLESKY_HPP
#define XF_SOLVER_CHOLESKY_HPP

#include <hls_stream.h>
#include <hls_math.h>
#include <ap_fixed.h>
#include <hls_x_complex.h>

// ===================== 并行度：强制覆盖（不依赖 tcl -D） =====================
#pragma message("CHOLESKY_HDR SWEETSPOT build: I_UNROLL=16, P_TILE=6")
#ifdef CHOLESKY_P_TILE
#undef CHOLESKY_P_TILE
#endif
#ifdef CHOLESKY_I_UNROLL
#undef CHOLESKY_I_UNROLL
#endif
// 甜点位（可自行改数值再综合）：P_TILE=6, I_UNROLL=8
#define CHOLESKY_P_TILE   6
#define CHOLESKY_I_UNROLL 16

#if 1
#if (CHOLESKY_P_TILE!=6) || (CHOLESKY_I_UNROLL!=16)
#error "This header is NOT in effect or macros are overridden elsewhere."
#endif
#endif

// 速度优先：使用 hls::rsqrt(float)。如需省 FPU，将其改为 0。
#ifdef CHOLESKY_USE_FLOAT_RSQRT
#undef CHOLESKY_USE_FLOAT_RSQRT
#endif
#define CHOLESKY_USE_FLOAT_RSQRT 1

namespace xf {
namespace solver {

//==================== 小工具 ====================//
namespace detail {
template <typename T> struct scalar_of { typedef T type; };
template <typename T> struct scalar_of<hls::x_complex<T> > { typedef T type; };

template <typename T>
inline hls::x_complex<T> conjx(const hls::x_complex<T>& x) {
#pragma HLS INLINE
    return hls::x_complex<T>(x.real(), -x.imag());
}
template <typename T>
inline hls::x_complex<T> make_cpx(const T& r, const T& i) {
#pragma HLS INLINE
    hls::x_complex<T> z; z.real(r); z.imag(i); return z;
}
template <typename T>
inline hls::x_complex<T> zero_c() {
#pragma HLS INLINE
    return make_cpx<T>((T)0, (T)0);
}
template <typename Narrow, typename Wider>
inline Narrow saturate_cast(const Wider& v) {
#pragma HLS INLINE
    return (Narrow)v;
}

// 3-乘法复数乘（Gauss trick）：节省 DSP，关键路径短
template <typename WS>
inline void cmul_3mul(WS a_r, WS a_i, WS b_r, WS b_i, WS& pr, WS& pi) {
#pragma HLS INLINE
    WS ac = a_r * b_r;
    WS bd = a_i * b_i;
    WS s  = (a_r + a_i) * (b_r + b_i);
    pr = ac - bd;
    pi = s  - ac - bd;
}

// 定点 rsqrt：若关闭浮点路径时使用
template <int W, int I>
inline ap_fixed<W, I> rsqrt_fixed(ap_fixed<W, I> a) {
#pragma HLS INLINE
    if (a <= (ap_fixed<W, I>)0) return (ap_fixed<W, I>)0;
    float af = (float)a;
    float y0 = hls::rsqrt(af);
    ap_fixed<W, I> y = (ap_fixed<W, I>)y0;
    ap_fixed<W+4, I+2> yy = (ap_fixed<W+4, I+2>)y * y;
    ap_fixed<W+4, I+2> t  = (ap_fixed<W+4, I+2>)1.5 - (ap_fixed<W+4, I+2>)0.5 * (ap_fixed<W+4, I+2>)a * yy;
    y = (ap_fixed<W+4, I+2>)y * t; // 1次NR，短拍
    return y;
}
} // namespace detail

//==================== Traits ====================//
template <typename Tin, typename Tout>
struct CholeskyTraits {
    typedef typename detail::scalar_of<Tin>::type  InScalar;
    typedef typename detail::scalar_of<Tout>::type OutScalar;

    static const int IN_W = InScalar::width;
    static const int IN_I = InScalar::iwidth;
    static const int ACC_W = (IN_W < 24) ? 32 : (IN_W + 8);
    static const int ACC_I = (IN_I < 12) ? 16 : (IN_I + 4);

    typedef ap_fixed<ACC_W, ACC_I>        WorkScalar;
    typedef hls::x_complex<WorkScalar>    WorkComplex;

    static const int P_TILE   = CHOLESKY_P_TILE;    // p维tile并行
    static const int I_UNROLL = CHOLESKY_I_UNROLL;  // 行并行
    static const bool USE_DATAFLOW = false;

    typedef hls::x_complex<InScalar>      IComplex;
    typedef hls::x_complex<OutScalar>     OComplex;
};

// 兼容层（若外部用 choleskyTraits<> 名称）
template<bool LOWER, int DIM, typename Tin, typename Tout>
struct choleskyTraits : public CholeskyTraits<Tin, Tout> {};

//==================== 主算法：融合对角 + 单遍 p-tile ====================//
template <bool LOWER_TRIANGULAR, int DIM, typename Tin, typename Tout, typename Traits = CholeskyTraits<Tin, Tout> >
int cholesky(hls::stream<Tin>& inA, hls::stream<Tout>& outL) {
#pragma HLS INLINE off

    typedef Traits                          TR;
    typedef typename TR::WorkScalar         WS;
    typedef typename TR::WorkComplex        WC;
    typedef typename TR::IComplex           IC;
    typedef typename TR::OComplex           OC;
    typedef typename TR::OutScalar          OS;

    // On-chip buffers
    IC A[DIM][DIM];
#pragma HLS ARRAY_PARTITION variable=A cyclic dim=2 factor=TR::P_TILE
#pragma HLS ARRAY_PARTITION variable=A cyclic dim=1 factor=TR::I_UNROLL
#pragma HLS BIND_STORAGE variable=A type=RAM_T2P impl=BRAM

    WC Llow[DIM][DIM];
#pragma HLS ARRAY_PARTITION variable=Llow complete dim=2
#pragma HLS ARRAY_PARTITION variable=Llow block factor=TR::I_UNROLL dim=1  // block 有助于就近布线

// 读入 A
read_in:
    for (int i = 0; i < DIM; ++i) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=DIM
        for (int j = 0; j < DIM; ++j) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=DIM
            A[i][j] = inA.read();
        }
    }

    int ret_code = 0;

// 列主序
col_loop:
    for (int j = 0; j < DIM; ++j) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=DIM
#pragma HLS DEPENDENCE variable=Llow inter false

        // 预取第 j 列到寄存器
        WC Aj[DIM];
#pragma HLS ARRAY_PARTITION variable=Aj complete
    pre_colA:
        for (int i = 0; i < DIM; ++i) {
#pragma HLS UNROLL factor=TR::I_UNROLL
            Aj[i] = (i > j) ? (WC)detail::make_cpx<WS>((WS)A[i][j].real(), (WS)A[i][j].imag())
                            : detail::zero_c<WS>();
        }

        // 行累计数组（每行一个累加器，先装 A[i][j]）
        WS acc_r_full[DIM], acc_i_full[DIM];
#pragma HLS ARRAY_PARTITION variable=acc_r_full block factor=TR::I_UNROLL
#pragma HLS ARRAY_PARTITION variable=acc_i_full block factor=TR::I_UNROLL
    init_acc_full:
        for (int i = 0; i < DIM; ++i) {
#pragma HLS UNROLL factor=TR::I_UNROLL
            acc_r_full[i] = Aj[i].real();
            acc_i_full[i] = Aj[i].imag();
        }

        // 对角起始
        WS diag_r = (WS)A[j][j].real();

        // ===== 单遍 p-tile：同时做对角与行点积 =====
    p_tiles_onepass:
        for (int p0 = 0; p0 < j; p0 += TR::P_TILE) {
#pragma HLS PIPELINE II=1
            // 取本tile的 conj(L[j][p])（先做共轭）
            WS b_r[TR::P_TILE], b_i[TR::P_TILE];
#pragma HLS ARRAY_PARTITION variable=b_r complete
#pragma HLS ARRAY_PARTITION variable=b_i complete
        load_b:
            for (int t = 0; t < TR::P_TILE; ++t) {
#pragma HLS UNROLL
                int p = p0 + t;
                if (p < j) {
                    b_r[t] = Llow[j][p].real();
                    b_i[t] = (WS)-Llow[j][p].imag();
                } else {
                    b_r[t] = (WS)0;  b_i[t] = (WS)0;
                }
            }

            // 1) 对角平方和（tile 归约）
            WS tile_diag = 0;
#pragma HLS DEPENDENCE variable=tile_diag inter false
#pragma HLS EXPRESSION_BALANCE
        tdiag:
            for (int t = 0; t < TR::P_TILE; ++t) {
#pragma HLS UNROLL
                int p = p0 + t;
                if (p < j) {
                    WS lr = Llow[j][p].real();
                    WS li = Llow[j][p].imag();
                    tile_diag += (WS)(lr * lr + li * li);
                }
            }
            diag_r -= tile_diag;

            // 2) 行更新（I_UNROLL 行并行 × P_TILE 内全展开）
        upd_rows_all:
            for (int i = j + 1; i < DIM; i += TR::I_UNROLL) {
#pragma HLS UNROLL
                for (int u = 0; u < TR::I_UNROLL; ++u) {
#pragma HLS UNROLL
                    int ii = i + u;
                    if (ii < DIM) {
                        WS sumr = 0, sumi = 0;
#pragma HLS DEPENDENCE variable=sumr inter false
#pragma HLS DEPENDENCE variable=sumi inter false
#pragma HLS EXPRESSION_BALANCE
                    mulacc_t:
                        for (int t = 0; t < TR::P_TILE; ++t) {
#pragma HLS UNROLL
                            int p = p0 + t;
                            if (p < j) {
                                WS a_r = Llow[ii][p].real();
                                WS a_i = Llow[ii][p].imag();
                                WS pr, pi;
                                detail::cmul_3mul<WS>(a_r, a_i, b_r[t], b_i[t], pr, pi);
                                sumr += pr;  sumi += pi;
                            }
                        }
                        acc_r_full[ii] -= sumr;
                        acc_i_full[ii] -= sumi;
                    }
                }
            }
        } // p_tiles_onepass

        if (diag_r < (WS)0) { ret_code = 1; diag_r = (WS)0; }

        // 求 inv_diag 与对角值
        WS inv_diag;
#if CHOLESKY_USE_FLOAT_RSQRT
        inv_diag = (WS)hls::rsqrt((float)diag_r);
#else
        inv_diag = detail::rsqrt_fixed<WS::width, WS::iwidth>(diag_r);
#endif
        WS diag_s = (WS)(diag_r * inv_diag);
        Llow[j][j] = detail::make_cpx<WS>(diag_s, (WS)0);

        // 回写该列下三角
    finalize_rows:
        for (int i = j + 1; i < DIM; i += TR::I_UNROLL) {
#pragma HLS PIPELINE II=1
            for (int u = 0; u < TR::I_UNROLL; ++u) {
#pragma HLS UNROLL
                int ii = i + u;
                if (ii < DIM) {
                    WS rr = acc_r_full[ii] * inv_diag;
                    WS ii_ = acc_i_full[ii] * inv_diag;
                    Llow[ii][j] = detail::make_cpx<WS>(rr, ii_);
                }
            }
        }
    } // col_loop

// 输出
write_out:
    for (int i = 0; i < DIM; ++i) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=DIM
        for (int j = 0; j < DIM; ++j) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=DIM
            OC outz;
            if (LOWER_TRIANGULAR) {
                if (i > j) {
                    OS rr = detail::saturate_cast<OS>(Llow[i][j].real());
                    OS ii = detail::saturate_cast<OS>(Llow[i][j].imag());
                    outz = detail::make_cpx<OS>(rr, ii);
                } else if (i == j) {
                    OS d = detail::saturate_cast<OS>(Llow[i][i].real());
                    outz = detail::make_cpx<OS>(d, (OS)0);
                } else {
                    outz = detail::zero_c<OS>();
                }
            } else { // 输出上三角（L^H）
                if (i < j) {
                    OS rr = detail::saturate_cast<OS>(Llow[j][i].real());
                    OS ii = detail::saturate_cast<OS>((WS)-Llow[j][i].imag());
                    outz = detail::make_cpx<OS>(rr, ii);
                } else if (i == j) {
                    OS d = detail::saturate_cast<OS>(Llow[i][i].real());
                    outz = detail::make_cpx<OS>(d, (OS)0);
                } else {
                    outz = detail::zero_c<OS>();
                }
            }
            outL.write((Tout)outz);
        }
    }

    return ret_code;
}

} // namespace solver
} // namespace xf

#endif // XF_SOLVER_CHOLESKY_HPP

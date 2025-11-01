#ifndef XF_SOLVER_CHOLESKY_HPP
#define XF_SOLVER_CHOLESKY_HPP

#include <hls_stream.h>
#include <hls_math.h>
#include <ap_fixed.h>
#include <hls_x_complex.h>

#ifndef CHOLESKY_INNER_UNROLL
#define CHOLESKY_INNER_UNROLL 4
#endif
#ifndef CHOLESKY_DIAG_UNROLL
#define CHOLESKY_DIAG_UNROLL 2
#endif

namespace xf {
namespace solver {

//==================== 工具：类型萃取与复数辅助 ====================//
namespace detail {

template <typename T> struct scalar_of { typedef T type; };
template <typename T> struct scalar_of<hls::x_complex<T> > { typedef T type; };

template <typename T>
inline hls::x_complex<T> conjx(const hls::x_complex<T>& x) {
#pragma HLS INLINE
    return hls::x_complex<T>(x.real(), -x.imag());
}

template <typename T>
inline T zero() {
#pragma HLS INLINE
    return T(0);
}

template <typename T>
inline hls::x_complex<T> zero_c() {
#pragma HLS INLINE
    return hls::x_complex<T>(T(0), T(0));
}

template <typename T>
inline hls::x_complex<T> make_cpx(const T& r, const T& i) {
#pragma HLS INLINE
    hls::x_complex<T> z;
    z.real(r);
    z.imag(i);
    return z;
}

template <typename Narrow, typename Wider>
inline Narrow saturate_cast(const Wider& v) {
#pragma HLS INLINE
    return (Narrow)v; // rely on ap_fixed rounding/saturation
}

template <typename T>
inline hls::x_complex<T> mul_conj(const hls::x_complex<T>& a, const hls::x_complex<T>& b) {
#pragma HLS INLINE
    T ar = a.real(), ai = a.imag();
    T br = b.real(), bi = b.imag();
    T rr = ar * br + ai * bi; // Re{a*conj(b)}
    T ii = ai * br - ar * bi; // Im{a*conj(b)}
    return make_cpx<T>(rr, ii);
}

} // namespace detail

//==================== Cholesky Traits（可按需调宽度） ====================//
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
        // === latency 调参开关（新增） ===
    static const int  INNER_UNROLL = 4;   // 行更新 sum_p 展开因子
    static const int  DIAG_UNROLL  = 3;   // 对角累加 acc_diag 展开因子
    static const bool USE_DATAFLOW = false; // 单帧测试先关
    typedef hls::x_complex<InScalar>      IComplex;
    typedef hls::x_complex<OutScalar>     OComplex;
};

// ======== 兼容层：保持旧签名 choleskyTraits<LOWER,DIM,Tin,Tout> 可用 ======== //
template<bool LOWER, int DIM, typename Tin, typename Tout>
struct choleskyTraits : public CholeskyTraits<Tin, Tout> {};

//==================== 主算法 ====================//
template <bool LOWER_TRIANGULAR, int DIM, typename Tin, typename Tout, typename Traits = CholeskyTraits<Tin, Tout> >
int cholesky(hls::stream<Tin>& inA, hls::stream<Tout>& outL) {
#pragma HLS INLINE off

    typedef Traits                          TR;
    typedef typename TR::WorkScalar         WS;
    typedef typename TR::WorkComplex        WC;
    typedef typename TR::IComplex           IC;
    typedef typename TR::OComplex           OC;
    typedef typename TR::OutScalar          OS;

    IC A[DIM][DIM];
#pragma HLS ARRAY_PARTITION variable=A dim=1 complete
#pragma HLS ARRAY_PARTITION variable=A complete dim=1
#pragma HLS ARRAY_PARTITION variable=A dim=2 complete

    WC Llow[DIM][DIM];
#pragma HLS ARRAY_PARTITION variable=Llow dim=1 complete
#pragma HLS ARRAY_PARTITION variable=Llow complete dim=1
#pragma HLS ARRAY_PARTITION variable=Llow dim=2 complete

// 读入
read_in:
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
#pragma HLS PIPELINE II=1
            A[i][j] = inA.read();
        }
    }

// (init_L removed: zeroing Llow eliminated to cut latency)
// 列主序 Cholesky：生成下三角 Llow
col_loop:
    for (int j = 0; j < DIM; ++j) {
#pragma HLS DEPENDENCE variable=Llow inter false
        WS diag_r = (WS)A[j][j].real();

acc_diag:
        for (int p = 0; p < j; ++p) {
                #pragma HLS PIPELINE II=1
            WS lr = Llow[j][p].real();
            WS li = Llow[j][p].imag();
            diag_r -= (WS)(lr * lr + li * li);
        }
        if (diag_r < (WS)0) diag_r = (WS)0;

        
float df = (float)((diag_r > (WS)0) ? diag_r : (WS)0);
float inv_f = (df > 0.0f) ? hls::rsqrt(df) : 0.0f;
WS diag_s = (WS)(df * inv_f);
WS inv_diag = (WS)inv_f;
Llow[j][j] = detail::make_cpx<WS>(diag_s, (WS)0);

row_update:
        for (int i = j + 1; i < DIM; ++i) {
#pragma HLS PIPELINE II=1
    #pragma HLS UNROLL factor=2
            WC acc = detail::make_cpx<WS>((WS)A[i][j].real(), (WS)A[i][j].imag());
sum_loop:
            for (int p = 0; p < j; ++p) {
    #pragma HLS PIPELINE II=1
                acc = acc - detail::mul_conj<WS>(Llow[i][p], Llow[j][p]);
            }
            WS rr = acc.real() * inv_diag;
            WS ii = acc.imag() * inv_diag;
            Llow[i][j] = detail::make_cpx<WS>(rr, ii);
        }
    }

// 写出（根据需求输出下三角 L 或上三角 U=Lᴴ）
write_out:
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
#pragma HLS PIPELINE II=1
            OC outz;
            if (LOWER_TRIANGULAR) {
                if (i == j) {
                    OS d = detail::saturate_cast<OS>(Llow[i][i].real());
                    outz = detail::make_cpx<OS>(d, (OS)0);
                } else if (i > j) {
                    OS r = detail::saturate_cast<OS>(Llow[i][j].real());
                    OS im = detail::saturate_cast<OS>(Llow[i][j].imag());
                    outz = detail::make_cpx<OS>(r, im);
                } else {
                    outz = detail::zero_c<OS>();
                }
            } else {
                if (j == i) {
                    OS d = detail::saturate_cast<OS>(Llow[i][i].real());
                    outz = detail::make_cpx<OS>(d, (OS)0);
                } else if (j > i) {
                    WS r = Llow[j][i].real();
                    WS im = Llow[j][i].imag();
                    OS rr = detail::saturate_cast<OS>(r);
                    OS ii = detail::saturate_cast<OS>((WS)-im);
                    outz = detail::make_cpx<OS>(rr, ii);
                } else {
                    outz = detail::zero_c<OS>();
                }
            }
            outL.write((Tout)outz);
        }
    }

    return 0;
}

} // namespace solver
} // namespace xf

#endif // XF_SOLVER_CHOL_OPT_HPP
/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) kinematic_solver_model_ ## ID
#endif

#include <math.h>
#include "include/kinematic_solver.h"

#ifndef casadi_real
#define casadi_real kinematic_solver_float
#endif

#ifndef casadi_int
#define casadi_int solver_int32_default
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_f1 CASADI_PREFIX(f1)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_s4 CASADI_PREFIX(s4)
#define casadi_s5 CASADI_PREFIX(s5)
#define casadi_s6 CASADI_PREFIX(s6)
#define casadi_sq CASADI_PREFIX(sq)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

casadi_real casadi_sq(casadi_real x) { return x*x;}

static const casadi_int casadi_s0[13] = {9, 1, 0, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8};
static const casadi_int casadi_s1[16] = {12, 1, 0, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
static const casadi_int casadi_s2[5] = {1, 1, 0, 1, 0};
static const casadi_int casadi_s3[18] = {1, 9, 0, 1, 2, 3, 4, 5, 5, 5, 5, 6, 0, 0, 0, 0, 0, 0};
static const casadi_int casadi_s4[15] = {1, 9, 0, 0, 0, 0, 1, 2, 2, 2, 2, 3, 0, 0, 0};
static const casadi_int casadi_s5[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
static const casadi_int casadi_s6[35] = {6, 9, 0, 4, 8, 9, 10, 11, 14, 18, 22, 23, 0, 1, 2, 3, 0, 1, 2, 4, 5, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 4, 5};

/* kinematic_solver_model_0_inner:(i0[9],i1[12])->(o0,o1[1x9,6nz],o2,o3[1x9,3nz],o4[6],o5[6x9,23nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a100, a101, a102, a103, a104, a105, a106, a107, a108, a109, a11, a110, a111, a112, a113, a114, a115, a116, a117, a118, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a6, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a7, a70, a71, a72, a73, a74, a75, a76, a77, a78, a79, a8, a80, a81, a82, a83, a84, a85, a86, a87, a88, a89, a9, a90, a91, a92, a93, a94, a95, a96, a97, a98, a99;
  a0=arg[1]? arg[1][3] : 0;
  a1=arg[1]? arg[1][0] : 0;
  a2=arg[1]? arg[1][4] : 0;
  a3=arg[0]? arg[0][8] : 0;
  a4=arg[1]? arg[1][5] : 0;
  a5=(a3-a4);
  a5=(a2*a5);
  a5=(a1+a5);
  a6=arg[0]? arg[0][3] : 0;
  a7=(a5-a6);
  a7=(a0*a7);
  a8=arg[1]? arg[1][1] : 0;
  a9=(a3-a4);
  a9=(a0*a9);
  a9=(a8+a9);
  a10=arg[0]? arg[0][4] : 0;
  a11=(a9-a10);
  a11=(a2*a11);
  a7=(a7-a11);
  a11=arg[1]? arg[1][6] : 0;
  a12=(a7*a11);
  a13=(a12*a7);
  a5=(a5-a6);
  a5=(a2*a5);
  a9=(a9-a10);
  a9=(a0*a9);
  a5=(a5+a9);
  a9=(a5*a11);
  a14=(a9*a5);
  a13=(a13+a14);
  a14=arg[1]? arg[1][8] : 0;
  a15=arg[0]? arg[0][2] : 0;
  a16=(a14*a15);
  a13=(a13-a16);
  a16=arg[0]? arg[0][0] : 0;
  a17=arg[1]? arg[1][9] : 0;
  a18=(a16*a17);
  a19=(a18*a16);
  a13=(a13+a19);
  a19=arg[0]? arg[0][1] : 0;
  a20=arg[1]? arg[1][10] : 0;
  a21=(a19*a20);
  a22=(a21*a19);
  a13=(a13+a22);
  if (res[0]!=0) res[0][0]=a13;
  a17=(a17*a16);
  a18=(a18+a17);
  if (res[1]!=0) res[1][0]=a18;
  a20=(a20*a19);
  a21=(a21+a20);
  if (res[1]!=0) res[1][1]=a21;
  a14=(-a14);
  if (res[1]!=0) res[1][2]=a14;
  a5=(a11*a5);
  a9=(a9+a5);
  a5=(a2*a9);
  a11=(a11*a7);
  a12=(a12+a11);
  a11=(a0*a12);
  a7=(a5+a11);
  a7=(-a7);
  if (res[1]!=0) res[1][3]=a7;
  a12=(a2*a12);
  a9=(a0*a9);
  a7=(a12-a9);
  if (res[1]!=0) res[1][4]=a7;
  a9=(a9-a12);
  a9=(a0*a9);
  a5=(a5+a11);
  a5=(a2*a5);
  a9=(a9+a5);
  if (res[1]!=0) res[1][5]=a9;
  a9=(a3-a4);
  a9=(a2*a9);
  a1=(a1+a9);
  a1=(a1-a6);
  a9=casadi_sq(a1);
  a4=(a3-a4);
  a4=(a0*a4);
  a8=(a8+a4);
  a8=(a8-a10);
  a4=casadi_sq(a8);
  a9=(a9+a4);
  a4=arg[1]? arg[1][11] : 0;
  a4=casadi_sq(a4);
  a9=(a9-a4);
  if (res[2]!=0) res[2][0]=a9;
  a1=(a1+a1);
  a9=(-a1);
  if (res[3]!=0) res[3][0]=a9;
  a8=(a8+a8);
  a9=(-a8);
  if (res[3]!=0) res[3][1]=a9;
  a0=(a0*a8);
  a2=(a2*a1);
  a0=(a0+a2);
  if (res[3]!=0) res[3][2]=a0;
  a0=2.7777777777777779e-03;
  a2=arg[0]? arg[0][6] : 0;
  a1=arg[0]? arg[0][5] : 0;
  a8=cos(a1);
  a9=(a2*a8);
  a4=2.;
  a5=8.3333333333333332e-03;
  a11=(a5*a16);
  a11=(a2+a11);
  a12=6.3000000000000000e-02;
  a7=(a2/a12);
  a14=arg[0]? arg[0][7] : 0;
  a21=tan(a14);
  a20=(a7*a21);
  a18=(a5*a20);
  a18=(a1+a18);
  a17=cos(a18);
  a13=(a11*a17);
  a13=(a4*a13);
  a9=(a9+a13);
  a13=(a5*a16);
  a13=(a2+a13);
  a22=(a11/a12);
  a23=(a5*a19);
  a23=(a14+a23);
  a24=tan(a23);
  a25=(a22*a24);
  a26=(a5*a25);
  a26=(a1+a26);
  a27=cos(a26);
  a28=(a13*a27);
  a28=(a4*a28);
  a9=(a9+a28);
  a28=1.6666666666666666e-02;
  a29=(a28*a16);
  a29=(a2+a29);
  a30=(a13/a12);
  a31=(a5*a19);
  a31=(a14+a31);
  a32=tan(a31);
  a33=(a30*a32);
  a34=(a28*a33);
  a34=(a1+a34);
  a35=cos(a34);
  a36=(a29*a35);
  a9=(a9+a36);
  a9=(a0*a9);
  a6=(a6+a9);
  a9=(a4*a16);
  a9=(a16+a9);
  a36=(a4*a16);
  a9=(a9+a36);
  a9=(a9+a16);
  a9=(a0*a9);
  a9=(a2+a9);
  a25=(a4*a25);
  a20=(a20+a25);
  a33=(a4*a33);
  a20=(a20+a33);
  a33=(a29/a12);
  a25=(a28*a19);
  a25=(a14+a25);
  a36=tan(a25);
  a37=(a33*a36);
  a20=(a20+a37);
  a20=(a0*a20);
  a20=(a1+a20);
  a37=cos(a20);
  a38=(a9*a37);
  a39=(a5*a16);
  a39=(a9+a39);
  a40=(a9/a12);
  a41=(a4*a19);
  a41=(a19+a41);
  a42=(a4*a19);
  a41=(a41+a42);
  a41=(a41+a19);
  a41=(a0*a41);
  a41=(a14+a41);
  a42=tan(a41);
  a43=(a40*a42);
  a44=(a5*a43);
  a44=(a20+a44);
  a45=cos(a44);
  a46=(a39*a45);
  a46=(a4*a46);
  a38=(a38+a46);
  a46=(a5*a16);
  a46=(a9+a46);
  a47=(a39/a12);
  a48=(a5*a19);
  a48=(a41+a48);
  a49=tan(a48);
  a50=(a47*a49);
  a51=(a5*a50);
  a51=(a20+a51);
  a52=cos(a51);
  a53=(a46*a52);
  a53=(a4*a53);
  a38=(a38+a53);
  a53=(a28*a16);
  a53=(a9+a53);
  a54=(a46/a12);
  a55=(a5*a19);
  a55=(a41+a55);
  a56=tan(a55);
  a57=(a54*a56);
  a58=(a28*a57);
  a58=(a20+a58);
  a59=cos(a58);
  a60=(a53*a59);
  a38=(a38+a60);
  a38=(a0*a38);
  a6=(a6+a38);
  a38=(a4*a16);
  a38=(a16+a38);
  a60=(a4*a16);
  a38=(a38+a60);
  a38=(a38+a16);
  a38=(a0*a38);
  a38=(a9+a38);
  a50=(a4*a50);
  a43=(a43+a50);
  a57=(a4*a57);
  a43=(a43+a57);
  a57=(a53/a12);
  a50=(a28*a19);
  a50=(a41+a50);
  a60=tan(a50);
  a61=(a57*a60);
  a43=(a43+a61);
  a43=(a0*a43);
  a43=(a20+a43);
  a61=cos(a43);
  a62=(a38*a61);
  a63=(a5*a16);
  a63=(a38+a63);
  a64=(a38/a12);
  a65=(a4*a19);
  a65=(a19+a65);
  a66=(a4*a19);
  a65=(a65+a66);
  a65=(a65+a19);
  a65=(a0*a65);
  a65=(a41+a65);
  a66=tan(a65);
  a67=(a64*a66);
  a68=(a5*a67);
  a68=(a43+a68);
  a69=cos(a68);
  a70=(a63*a69);
  a70=(a4*a70);
  a62=(a62+a70);
  a70=(a5*a16);
  a70=(a38+a70);
  a71=(a63/a12);
  a72=(a5*a19);
  a72=(a65+a72);
  a73=tan(a72);
  a74=(a71*a73);
  a75=(a5*a74);
  a75=(a43+a75);
  a76=cos(a75);
  a77=(a70*a76);
  a77=(a4*a77);
  a62=(a62+a77);
  a77=(a28*a16);
  a77=(a38+a77);
  a78=(a70/a12);
  a79=(a5*a19);
  a79=(a65+a79);
  a80=tan(a79);
  a81=(a78*a80);
  a82=(a28*a81);
  a82=(a43+a82);
  a83=cos(a82);
  a84=(a77*a83);
  a62=(a62+a84);
  a62=(a0*a62);
  a6=(a6+a62);
  if (res[4]!=0) res[4][0]=a6;
  a6=sin(a1);
  a62=(a2*a6);
  a84=sin(a18);
  a85=(a11*a84);
  a85=(a4*a85);
  a62=(a62+a85);
  a85=sin(a26);
  a86=(a13*a85);
  a86=(a4*a86);
  a62=(a62+a86);
  a86=sin(a34);
  a87=(a29*a86);
  a62=(a62+a87);
  a62=(a0*a62);
  a10=(a10+a62);
  a62=sin(a20);
  a87=(a9*a62);
  a88=sin(a44);
  a89=(a39*a88);
  a89=(a4*a89);
  a87=(a87+a89);
  a89=sin(a51);
  a90=(a46*a89);
  a90=(a4*a90);
  a87=(a87+a90);
  a90=sin(a58);
  a91=(a53*a90);
  a87=(a87+a91);
  a87=(a0*a87);
  a10=(a10+a87);
  a87=sin(a43);
  a91=(a38*a87);
  a92=sin(a68);
  a93=(a63*a92);
  a93=(a4*a93);
  a91=(a91+a93);
  a93=sin(a75);
  a94=(a70*a93);
  a94=(a4*a94);
  a91=(a91+a94);
  a94=sin(a82);
  a95=(a77*a94);
  a91=(a91+a95);
  a91=(a0*a91);
  a10=(a10+a91);
  if (res[4]!=0) res[4][1]=a10;
  a74=(a4*a74);
  a67=(a67+a74);
  a81=(a4*a81);
  a67=(a67+a81);
  a12=(a77/a12);
  a81=(a28*a19);
  a81=(a65+a81);
  a74=tan(a81);
  a10=(a12*a74);
  a67=(a67+a10);
  a67=(a0*a67);
  a67=(a43+a67);
  if (res[4]!=0) res[4][2]=a67;
  a67=(a4*a16);
  a67=(a16+a67);
  a10=(a4*a16);
  a67=(a67+a10);
  a67=(a67+a16);
  a67=(a0*a67);
  a67=(a38+a67);
  if (res[4]!=0) res[4][3]=a67;
  a67=(a4*a19);
  a67=(a19+a67);
  a16=(a4*a19);
  a67=(a67+a16);
  a67=(a67+a19);
  a67=(a0*a67);
  a67=(a65+a67);
  if (res[4]!=0) res[4][4]=a67;
  a67=(a4*a15);
  a67=(a15+a67);
  a19=(a4*a15);
  a67=(a67+a19);
  a67=(a67+a15);
  a67=(a0*a67);
  a3=(a3+a67);
  a67=(a4*a15);
  a67=(a15+a67);
  a19=(a4*a15);
  a67=(a67+a19);
  a67=(a67+a15);
  a67=(a0*a67);
  a3=(a3+a67);
  a67=(a4*a15);
  a67=(a15+a67);
  a19=(a4*a15);
  a67=(a67+a19);
  a67=(a67+a15);
  a67=(a0*a67);
  a3=(a3+a67);
  if (res[4]!=0) res[4][5]=a3;
  a3=(a5*a17);
  a3=(a4*a3);
  a67=(a5*a27);
  a15=sin(a26);
  a19=1.3227513227513227e-01;
  a16=(a19*a24);
  a10=(a5*a16);
  a91=(a15*a10);
  a91=(a13*a91);
  a67=(a67-a91);
  a67=(a4*a67);
  a3=(a3+a67);
  a67=(a28*a35);
  a91=sin(a34);
  a19=(a19*a32);
  a95=(a28*a19);
  a96=(a91*a95);
  a96=(a29*a96);
  a67=(a67-a96);
  a3=(a3+a67);
  a3=(a0*a3);
  a67=(a28*a37);
  a96=sin(a20);
  a16=(a4*a16);
  a19=(a4*a19);
  a16=(a16+a19);
  a19=2.6455026455026454e-01;
  a97=(a19*a36);
  a16=(a16+a97);
  a16=(a0*a16);
  a97=(a96*a16);
  a97=(a9*a97);
  a67=(a67-a97);
  a97=2.5000000000000001e-02;
  a98=(a97*a45);
  a99=sin(a44);
  a19=(a19*a42);
  a100=(a5*a19);
  a100=(a16+a100);
  a101=(a99*a100);
  a101=(a39*a101);
  a98=(a98-a101);
  a98=(a4*a98);
  a67=(a67+a98);
  a98=(a97*a52);
  a101=sin(a51);
  a102=3.9682539682539686e-01;
  a103=(a102*a49);
  a104=(a5*a103);
  a104=(a16+a104);
  a105=(a101*a104);
  a105=(a46*a105);
  a98=(a98-a105);
  a98=(a4*a98);
  a67=(a67+a98);
  a98=3.3333333333333333e-02;
  a105=(a98*a59);
  a106=sin(a58);
  a102=(a102*a56);
  a107=(a28*a102);
  a107=(a16+a107);
  a108=(a106*a107);
  a108=(a53*a108);
  a105=(a105-a108);
  a67=(a67+a105);
  a67=(a0*a67);
  a3=(a3+a67);
  a67=(a98*a61);
  a105=sin(a43);
  a103=(a4*a103);
  a19=(a19+a103);
  a102=(a4*a102);
  a19=(a19+a102);
  a102=5.2910052910052907e-01;
  a103=(a102*a60);
  a19=(a19+a103);
  a19=(a0*a19);
  a19=(a16+a19);
  a103=(a105*a19);
  a103=(a38*a103);
  a67=(a67-a103);
  a103=4.1666666666666664e-02;
  a108=(a103*a69);
  a109=sin(a68);
  a102=(a102*a66);
  a110=(a5*a102);
  a110=(a19+a110);
  a111=(a109*a110);
  a111=(a63*a111);
  a108=(a108-a111);
  a108=(a4*a108);
  a67=(a67+a108);
  a108=(a103*a76);
  a111=sin(a75);
  a112=6.6137566137566139e-01;
  a113=(a112*a73);
  a114=(a5*a113);
  a114=(a19+a114);
  a115=(a111*a114);
  a115=(a70*a115);
  a108=(a108-a115);
  a108=(a4*a108);
  a67=(a67+a108);
  a108=5.0000000000000003e-02;
  a115=(a108*a83);
  a116=sin(a82);
  a112=(a112*a80);
  a117=(a28*a112);
  a117=(a19+a117);
  a118=(a116*a117);
  a118=(a77*a118);
  a115=(a115-a118);
  a67=(a67+a115);
  a67=(a0*a67);
  a3=(a3+a67);
  if (res[5]!=0) res[5][0]=a3;
  a3=(a5*a84);
  a3=(a4*a3);
  a67=(a5*a85);
  a26=cos(a26);
  a10=(a26*a10);
  a10=(a13*a10);
  a67=(a67+a10);
  a67=(a4*a67);
  a3=(a3+a67);
  a67=(a28*a86);
  a34=cos(a34);
  a95=(a34*a95);
  a95=(a29*a95);
  a67=(a67+a95);
  a3=(a3+a67);
  a3=(a0*a3);
  a67=(a28*a62);
  a20=cos(a20);
  a16=(a20*a16);
  a16=(a9*a16);
  a67=(a67+a16);
  a16=(a97*a88);
  a44=cos(a44);
  a100=(a44*a100);
  a100=(a39*a100);
  a16=(a16+a100);
  a16=(a4*a16);
  a67=(a67+a16);
  a16=(a97*a89);
  a51=cos(a51);
  a104=(a51*a104);
  a104=(a46*a104);
  a16=(a16+a104);
  a16=(a4*a16);
  a67=(a67+a16);
  a16=(a98*a90);
  a58=cos(a58);
  a107=(a58*a107);
  a107=(a53*a107);
  a16=(a16+a107);
  a67=(a67+a16);
  a67=(a0*a67);
  a3=(a3+a67);
  a67=(a98*a87);
  a43=cos(a43);
  a16=(a43*a19);
  a16=(a38*a16);
  a67=(a67+a16);
  a16=(a103*a92);
  a68=cos(a68);
  a110=(a68*a110);
  a110=(a63*a110);
  a16=(a16+a110);
  a16=(a4*a16);
  a67=(a67+a16);
  a16=(a103*a93);
  a75=cos(a75);
  a114=(a75*a114);
  a114=(a70*a114);
  a16=(a16+a114);
  a16=(a4*a16);
  a67=(a67+a16);
  a16=(a108*a94);
  a82=cos(a82);
  a117=(a82*a117);
  a117=(a77*a117);
  a16=(a16+a117);
  a67=(a67+a16);
  a67=(a0*a67);
  a3=(a3+a67);
  if (res[5]!=0) res[5][1]=a3;
  a113=(a4*a113);
  a102=(a102+a113);
  a112=(a4*a112);
  a102=(a102+a112);
  a112=7.9365079365079372e-01;
  a112=(a112*a74);
  a102=(a102+a112);
  a102=(a0*a102);
  a19=(a19+a102);
  if (res[5]!=0) res[5][2]=a19;
  if (res[5]!=0) res[5][3]=a108;
  a23=cos(a23);
  a23=casadi_sq(a23);
  a19=(a5/a23);
  a19=(a22*a19);
  a102=(a5*a19);
  a112=(a15*a102);
  a112=(a13*a112);
  a112=(a4*a112);
  a31=cos(a31);
  a31=casadi_sq(a31);
  a113=(a5/a31);
  a113=(a30*a113);
  a3=(a28*a113);
  a67=(a91*a3);
  a67=(a29*a67);
  a112=(a112+a67);
  a112=(a0*a112);
  a19=(a4*a19);
  a113=(a4*a113);
  a19=(a19+a113);
  a25=cos(a25);
  a25=casadi_sq(a25);
  a113=(a28/a25);
  a113=(a33*a113);
  a19=(a19+a113);
  a19=(a0*a19);
  a113=(a96*a19);
  a113=(a9*a113);
  a41=cos(a41);
  a41=casadi_sq(a41);
  a67=(a28/a41);
  a67=(a40*a67);
  a16=(a5*a67);
  a16=(a19+a16);
  a117=(a99*a16);
  a117=(a39*a117);
  a117=(a4*a117);
  a113=(a113+a117);
  a48=cos(a48);
  a48=casadi_sq(a48);
  a117=(a97/a48);
  a117=(a47*a117);
  a114=(a5*a117);
  a114=(a19+a114);
  a110=(a101*a114);
  a110=(a46*a110);
  a110=(a4*a110);
  a113=(a113+a110);
  a55=cos(a55);
  a55=casadi_sq(a55);
  a97=(a97/a55);
  a97=(a54*a97);
  a110=(a28*a97);
  a110=(a19+a110);
  a107=(a106*a110);
  a107=(a53*a107);
  a113=(a113+a107);
  a113=(a0*a113);
  a112=(a112+a113);
  a117=(a4*a117);
  a67=(a67+a117);
  a97=(a4*a97);
  a67=(a67+a97);
  a50=cos(a50);
  a50=casadi_sq(a50);
  a97=(a98/a50);
  a97=(a57*a97);
  a67=(a67+a97);
  a67=(a0*a67);
  a67=(a19+a67);
  a97=(a105*a67);
  a97=(a38*a97);
  a65=cos(a65);
  a65=casadi_sq(a65);
  a98=(a98/a65);
  a98=(a64*a98);
  a117=(a5*a98);
  a117=(a67+a117);
  a113=(a109*a117);
  a113=(a63*a113);
  a113=(a4*a113);
  a97=(a97+a113);
  a72=cos(a72);
  a72=casadi_sq(a72);
  a113=(a103/a72);
  a113=(a71*a113);
  a107=(a5*a113);
  a107=(a67+a107);
  a104=(a111*a107);
  a104=(a70*a104);
  a104=(a4*a104);
  a97=(a97+a104);
  a79=cos(a79);
  a79=casadi_sq(a79);
  a103=(a103/a79);
  a103=(a78*a103);
  a104=(a28*a103);
  a104=(a67+a104);
  a100=(a116*a104);
  a100=(a77*a100);
  a97=(a97+a100);
  a97=(a0*a97);
  a112=(a112+a97);
  a112=(-a112);
  if (res[5]!=0) res[5][4]=a112;
  a102=(a26*a102);
  a102=(a13*a102);
  a102=(a4*a102);
  a3=(a34*a3);
  a3=(a29*a3);
  a102=(a102+a3);
  a102=(a0*a102);
  a19=(a20*a19);
  a19=(a9*a19);
  a16=(a44*a16);
  a16=(a39*a16);
  a16=(a4*a16);
  a19=(a19+a16);
  a114=(a51*a114);
  a114=(a46*a114);
  a114=(a4*a114);
  a19=(a19+a114);
  a110=(a58*a110);
  a110=(a53*a110);
  a19=(a19+a110);
  a19=(a0*a19);
  a102=(a102+a19);
  a19=(a43*a67);
  a19=(a38*a19);
  a117=(a68*a117);
  a117=(a63*a117);
  a117=(a4*a117);
  a19=(a19+a117);
  a107=(a75*a107);
  a107=(a70*a107);
  a107=(a4*a107);
  a19=(a19+a107);
  a104=(a82*a104);
  a104=(a77*a104);
  a19=(a19+a104);
  a19=(a0*a19);
  a102=(a102+a19);
  if (res[5]!=0) res[5][5]=a102;
  a113=(a4*a113);
  a98=(a98+a113);
  a103=(a4*a103);
  a98=(a98+a103);
  a81=cos(a81);
  a81=casadi_sq(a81);
  a103=(a108/a81);
  a103=(a12*a103);
  a98=(a98+a103);
  a98=(a0*a98);
  a67=(a67+a98);
  if (res[5]!=0) res[5][6]=a67;
  if (res[5]!=0) res[5][7]=a108;
  if (res[5]!=0) res[5][8]=a108;
  a108=1.;
  if (res[5]!=0) res[5][9]=a108;
  if (res[5]!=0) res[5][10]=a108;
  a67=sin(a1);
  a67=(a2*a67);
  a98=sin(a18);
  a103=(a11*a98);
  a103=(a4*a103);
  a67=(a67+a103);
  a103=(a13*a15);
  a103=(a4*a103);
  a67=(a67+a103);
  a103=(a29*a91);
  a67=(a67+a103);
  a67=(a0*a67);
  a103=(a9*a96);
  a113=(a39*a99);
  a113=(a4*a113);
  a103=(a103+a113);
  a113=(a46*a101);
  a113=(a4*a113);
  a103=(a103+a113);
  a113=(a53*a106);
  a103=(a103+a113);
  a103=(a0*a103);
  a67=(a67+a103);
  a103=(a38*a105);
  a113=(a63*a109);
  a113=(a4*a113);
  a103=(a103+a113);
  a113=(a70*a111);
  a113=(a4*a113);
  a103=(a103+a113);
  a113=(a77*a116);
  a103=(a103+a113);
  a103=(a0*a103);
  a67=(a67+a103);
  a67=(-a67);
  if (res[5]!=0) res[5][11]=a67;
  a1=cos(a1);
  a2=(a2*a1);
  a18=cos(a18);
  a1=(a11*a18);
  a1=(a4*a1);
  a2=(a2+a1);
  a1=(a13*a26);
  a1=(a4*a1);
  a2=(a2+a1);
  a1=(a29*a34);
  a2=(a2+a1);
  a2=(a0*a2);
  a1=(a9*a20);
  a67=(a39*a44);
  a67=(a4*a67);
  a1=(a1+a67);
  a67=(a46*a51);
  a67=(a4*a67);
  a1=(a1+a67);
  a67=(a53*a58);
  a1=(a1+a67);
  a1=(a0*a1);
  a2=(a2+a1);
  a1=(a38*a43);
  a67=(a63*a68);
  a67=(a4*a67);
  a1=(a1+a67);
  a67=(a70*a75);
  a67=(a4*a67);
  a1=(a1+a67);
  a67=(a77*a82);
  a1=(a1+a67);
  a1=(a0*a1);
  a2=(a2+a1);
  if (res[5]!=0) res[5][12]=a2;
  if (res[5]!=0) res[5][13]=a108;
  a2=1.5873015873015873e+01;
  a21=(a2*a21);
  a1=(a5*a21);
  a67=(a98*a1);
  a67=(a11*a67);
  a17=(a17-a67);
  a17=(a4*a17);
  a8=(a8+a17);
  a24=(a2*a24);
  a17=(a5*a24);
  a67=(a15*a17);
  a67=(a13*a67);
  a27=(a27-a67);
  a27=(a4*a27);
  a8=(a8+a27);
  a32=(a2*a32);
  a27=(a28*a32);
  a67=(a91*a27);
  a67=(a29*a67);
  a35=(a35-a67);
  a8=(a8+a35);
  a8=(a0*a8);
  a24=(a4*a24);
  a21=(a21+a24);
  a32=(a4*a32);
  a21=(a21+a32);
  a36=(a2*a36);
  a21=(a21+a36);
  a21=(a0*a21);
  a36=(a96*a21);
  a36=(a9*a36);
  a37=(a37-a36);
  a42=(a2*a42);
  a36=(a5*a42);
  a36=(a21+a36);
  a32=(a99*a36);
  a32=(a39*a32);
  a45=(a45-a32);
  a45=(a4*a45);
  a37=(a37+a45);
  a49=(a2*a49);
  a45=(a5*a49);
  a45=(a21+a45);
  a32=(a101*a45);
  a32=(a46*a32);
  a52=(a52-a32);
  a52=(a4*a52);
  a37=(a37+a52);
  a56=(a2*a56);
  a52=(a28*a56);
  a52=(a21+a52);
  a32=(a106*a52);
  a32=(a53*a32);
  a59=(a59-a32);
  a37=(a37+a59);
  a37=(a0*a37);
  a8=(a8+a37);
  a49=(a4*a49);
  a42=(a42+a49);
  a56=(a4*a56);
  a42=(a42+a56);
  a60=(a2*a60);
  a42=(a42+a60);
  a42=(a0*a42);
  a42=(a21+a42);
  a60=(a105*a42);
  a60=(a38*a60);
  a61=(a61-a60);
  a66=(a2*a66);
  a60=(a5*a66);
  a60=(a42+a60);
  a56=(a109*a60);
  a56=(a63*a56);
  a69=(a69-a56);
  a69=(a4*a69);
  a61=(a61+a69);
  a73=(a2*a73);
  a69=(a5*a73);
  a69=(a42+a69);
  a56=(a111*a69);
  a56=(a70*a56);
  a76=(a76-a56);
  a76=(a4*a76);
  a61=(a61+a76);
  a80=(a2*a80);
  a76=(a28*a80);
  a76=(a42+a76);
  a56=(a116*a76);
  a56=(a77*a56);
  a83=(a83-a56);
  a61=(a61+a83);
  a61=(a0*a61);
  a8=(a8+a61);
  if (res[5]!=0) res[5][14]=a8;
  a1=(a18*a1);
  a1=(a11*a1);
  a84=(a84+a1);
  a84=(a4*a84);
  a6=(a6+a84);
  a17=(a26*a17);
  a17=(a13*a17);
  a85=(a85+a17);
  a85=(a4*a85);
  a6=(a6+a85);
  a27=(a34*a27);
  a27=(a29*a27);
  a86=(a86+a27);
  a6=(a6+a86);
  a6=(a0*a6);
  a21=(a20*a21);
  a21=(a9*a21);
  a62=(a62+a21);
  a36=(a44*a36);
  a36=(a39*a36);
  a88=(a88+a36);
  a88=(a4*a88);
  a62=(a62+a88);
  a45=(a51*a45);
  a45=(a46*a45);
  a89=(a89+a45);
  a89=(a4*a89);
  a62=(a62+a89);
  a52=(a58*a52);
  a52=(a53*a52);
  a90=(a90+a52);
  a62=(a62+a90);
  a62=(a0*a62);
  a6=(a6+a62);
  a62=(a43*a42);
  a62=(a38*a62);
  a87=(a87+a62);
  a60=(a68*a60);
  a60=(a63*a60);
  a92=(a92+a60);
  a92=(a4*a92);
  a87=(a87+a92);
  a69=(a75*a69);
  a69=(a70*a69);
  a93=(a93+a69);
  a93=(a4*a93);
  a87=(a87+a93);
  a76=(a82*a76);
  a76=(a77*a76);
  a94=(a94+a76);
  a87=(a87+a94);
  a87=(a0*a87);
  a6=(a6+a87);
  if (res[5]!=0) res[5][15]=a6;
  a73=(a4*a73);
  a66=(a66+a73);
  a80=(a4*a80);
  a66=(a66+a80);
  a2=(a2*a74);
  a66=(a66+a2);
  a66=(a0*a66);
  a42=(a42+a66);
  if (res[5]!=0) res[5][16]=a42;
  if (res[5]!=0) res[5][17]=a108;
  a14=cos(a14);
  a14=casadi_sq(a14);
  a7=(a7/a14);
  a14=(a5*a7);
  a98=(a98*a14);
  a98=(a11*a98);
  a98=(a4*a98);
  a22=(a22/a23);
  a23=(a5*a22);
  a15=(a15*a23);
  a15=(a13*a15);
  a15=(a4*a15);
  a98=(a98+a15);
  a30=(a30/a31);
  a31=(a28*a30);
  a91=(a91*a31);
  a91=(a29*a91);
  a98=(a98+a91);
  a98=(a0*a98);
  a22=(a4*a22);
  a7=(a7+a22);
  a30=(a4*a30);
  a7=(a7+a30);
  a33=(a33/a25);
  a7=(a7+a33);
  a7=(a0*a7);
  a96=(a96*a7);
  a96=(a9*a96);
  a40=(a40/a41);
  a41=(a5*a40);
  a41=(a7+a41);
  a99=(a99*a41);
  a99=(a39*a99);
  a99=(a4*a99);
  a96=(a96+a99);
  a47=(a47/a48);
  a48=(a5*a47);
  a48=(a7+a48);
  a101=(a101*a48);
  a101=(a46*a101);
  a101=(a4*a101);
  a96=(a96+a101);
  a54=(a54/a55);
  a55=(a28*a54);
  a55=(a7+a55);
  a106=(a106*a55);
  a106=(a53*a106);
  a96=(a96+a106);
  a96=(a0*a96);
  a98=(a98+a96);
  a47=(a4*a47);
  a40=(a40+a47);
  a54=(a4*a54);
  a40=(a40+a54);
  a57=(a57/a50);
  a40=(a40+a57);
  a40=(a0*a40);
  a40=(a7+a40);
  a105=(a105*a40);
  a105=(a38*a105);
  a64=(a64/a65);
  a65=(a5*a64);
  a65=(a40+a65);
  a109=(a109*a65);
  a109=(a63*a109);
  a109=(a4*a109);
  a105=(a105+a109);
  a71=(a71/a72);
  a5=(a5*a71);
  a5=(a40+a5);
  a111=(a111*a5);
  a111=(a70*a111);
  a111=(a4*a111);
  a105=(a105+a111);
  a78=(a78/a79);
  a28=(a28*a78);
  a28=(a40+a28);
  a116=(a116*a28);
  a116=(a77*a116);
  a105=(a105+a116);
  a105=(a0*a105);
  a98=(a98+a105);
  a98=(-a98);
  if (res[5]!=0) res[5][18]=a98;
  a18=(a18*a14);
  a11=(a11*a18);
  a11=(a4*a11);
  a26=(a26*a23);
  a13=(a13*a26);
  a13=(a4*a13);
  a11=(a11+a13);
  a34=(a34*a31);
  a29=(a29*a34);
  a11=(a11+a29);
  a11=(a0*a11);
  a20=(a20*a7);
  a9=(a9*a20);
  a44=(a44*a41);
  a39=(a39*a44);
  a39=(a4*a39);
  a9=(a9+a39);
  a51=(a51*a48);
  a46=(a46*a51);
  a46=(a4*a46);
  a9=(a9+a46);
  a58=(a58*a55);
  a53=(a53*a58);
  a9=(a9+a53);
  a9=(a0*a9);
  a11=(a11+a9);
  a43=(a43*a40);
  a38=(a38*a43);
  a68=(a68*a65);
  a63=(a63*a68);
  a63=(a4*a63);
  a38=(a38+a63);
  a75=(a75*a5);
  a70=(a70*a75);
  a70=(a4*a70);
  a38=(a38+a70);
  a82=(a82*a28);
  a77=(a77*a82);
  a38=(a38+a77);
  a38=(a0*a38);
  a11=(a11+a38);
  if (res[5]!=0) res[5][19]=a11;
  a71=(a4*a71);
  a64=(a64+a71);
  a4=(a4*a78);
  a64=(a64+a4);
  a12=(a12/a81);
  a64=(a64+a12);
  a0=(a0*a64);
  a40=(a40+a0);
  if (res[5]!=0) res[5][20]=a40;
  if (res[5]!=0) res[5][21]=a108;
  if (res[5]!=0) res[5][22]=a108;
  return 0;
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_0_inner(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_0_inner_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_0_inner_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_0_inner_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_0_inner_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_0_inner_release(int mem) {
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_0_inner_incref(void) {
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_0_inner_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int kinematic_solver_model_0_inner_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int kinematic_solver_model_0_inner_n_out(void) { return 6;}

CASADI_SYMBOL_EXPORT casadi_real kinematic_solver_model_0_inner_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* kinematic_solver_model_0_inner_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* kinematic_solver_model_0_inner_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    case 2: return "o2";
    case 3: return "o3";
    case 4: return "o4";
    case 5: return "o5";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* kinematic_solver_model_0_inner_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* kinematic_solver_model_0_inner_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    case 1: return casadi_s3;
    case 2: return casadi_s2;
    case 3: return casadi_s4;
    case 4: return casadi_s5;
    case 5: return casadi_s6;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_0_inner_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 6;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* kinematic_solver_model_1_inner:(i0[9],i1[12])->(o0,o1[1x9,6nz],o2,o3[1x9,3nz]) */
static int casadi_f1(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[1]? arg[1][3] : 0;
  a1=arg[1]? arg[1][0] : 0;
  a2=arg[1]? arg[1][4] : 0;
  a3=arg[0]? arg[0][8] : 0;
  a4=arg[1]? arg[1][5] : 0;
  a5=(a3-a4);
  a5=(a2*a5);
  a5=(a1+a5);
  a6=arg[0]? arg[0][3] : 0;
  a7=(a5-a6);
  a7=(a0*a7);
  a8=arg[1]? arg[1][1] : 0;
  a9=(a3-a4);
  a9=(a0*a9);
  a9=(a8+a9);
  a10=arg[0]? arg[0][4] : 0;
  a11=(a9-a10);
  a11=(a2*a11);
  a7=(a7-a11);
  a11=arg[1]? arg[1][6] : 0;
  a12=(a7*a11);
  a13=(a12*a7);
  a5=(a5-a6);
  a5=(a2*a5);
  a9=(a9-a10);
  a9=(a0*a9);
  a5=(a5+a9);
  a9=(a5*a11);
  a14=(a9*a5);
  a13=(a13+a14);
  a14=arg[1]? arg[1][8] : 0;
  a15=arg[0]? arg[0][2] : 0;
  a15=(a14*a15);
  a13=(a13-a15);
  a15=arg[0]? arg[0][0] : 0;
  a16=arg[1]? arg[1][9] : 0;
  a17=(a15*a16);
  a18=(a17*a15);
  a13=(a13+a18);
  a18=arg[0]? arg[0][1] : 0;
  a19=arg[1]? arg[1][10] : 0;
  a20=(a18*a19);
  a21=(a20*a18);
  a13=(a13+a21);
  if (res[0]!=0) res[0][0]=a13;
  a16=(a16*a15);
  a17=(a17+a16);
  if (res[1]!=0) res[1][0]=a17;
  a19=(a19*a18);
  a20=(a20+a19);
  if (res[1]!=0) res[1][1]=a20;
  a14=(-a14);
  if (res[1]!=0) res[1][2]=a14;
  a5=(a11*a5);
  a9=(a9+a5);
  a5=(a2*a9);
  a11=(a11*a7);
  a12=(a12+a11);
  a11=(a0*a12);
  a7=(a5+a11);
  a7=(-a7);
  if (res[1]!=0) res[1][3]=a7;
  a12=(a2*a12);
  a9=(a0*a9);
  a7=(a12-a9);
  if (res[1]!=0) res[1][4]=a7;
  a9=(a9-a12);
  a9=(a0*a9);
  a5=(a5+a11);
  a5=(a2*a5);
  a9=(a9+a5);
  if (res[1]!=0) res[1][5]=a9;
  a9=(a3-a4);
  a9=(a2*a9);
  a1=(a1+a9);
  a1=(a1-a6);
  a6=casadi_sq(a1);
  a3=(a3-a4);
  a3=(a0*a3);
  a8=(a8+a3);
  a8=(a8-a10);
  a10=casadi_sq(a8);
  a6=(a6+a10);
  a10=arg[1]? arg[1][11] : 0;
  a10=casadi_sq(a10);
  a6=(a6-a10);
  if (res[2]!=0) res[2][0]=a6;
  a1=(a1+a1);
  a6=(-a1);
  if (res[3]!=0) res[3][0]=a6;
  a8=(a8+a8);
  a6=(-a8);
  if (res[3]!=0) res[3][1]=a6;
  a0=(a0*a8);
  a2=(a2*a1);
  a0=(a0+a2);
  if (res[3]!=0) res[3][2]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_1_inner(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f1(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_1_inner_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_1_inner_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_1_inner_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_1_inner_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_1_inner_release(int mem) {
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_1_inner_incref(void) {
}

CASADI_SYMBOL_EXPORT void kinematic_solver_model_1_inner_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int kinematic_solver_model_1_inner_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int kinematic_solver_model_1_inner_n_out(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_real kinematic_solver_model_1_inner_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* kinematic_solver_model_1_inner_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* kinematic_solver_model_1_inner_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    case 2: return "o2";
    case 3: return "o3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* kinematic_solver_model_1_inner_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* kinematic_solver_model_1_inner_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    case 1: return casadi_s3;
    case 2: return casadi_s2;
    case 3: return casadi_s4;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int kinematic_solver_model_1_inner_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 4;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif

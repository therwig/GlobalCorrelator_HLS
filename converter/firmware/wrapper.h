#ifndef WRAPPER_H
#define WRAPPER_H

#include <iostream>
#include <cmath>
#include "ap_int.h"
#include "ap_fixed.h"
#include "../../firmware/data.h"
#include "../../l1tk_types.h"
using namespace l1tk;
//#include "../../DiscretePFInputs.h"

#define PT_INV_TAB_SIZE 8
#define PT_INV_MAX_BITS 8
#define ETA_TAB_SIZE 8
#define DPHI_TAB_SIZE 8

typedef ap_uint<96> input_t;
typedef ap_uint<64> output_t;
void pf_input_track_conv_hw(input_t in, output_t& out, numlink_t nlink);

typedef ap_fixed<24,12, AP_RND_CONV, AP_SAT> bigfix_t; // helper type

typedef ap_fixed<17,10, AP_RND_CONV, AP_SAT> zdet_t; // helper type
typedef ap_fixed<18,6, AP_RND_CONV, AP_SAT> tanlam_help_t; // helper type

#define PF_PT_SCALE (4.0)
#define PF_ETAPHI_SCALE (4*180./3.1415)
#define PF_Z0_SCALE (20)
#define TT_NPHI_SECTORS 9

#define DETR 129
#define DETZ 320.9

void pack_L1T_track(ap_uint<kTrackWordSize> &tk,
                    rinv_t     rinv    ,
                    tkphi_t    tkphi   ,
                    tanlam_t   tanlam  ,
                    tkz0_t     tkz0    ,
                    tkd0_t     tkd0    ,
                    chi2rphi_t chi2rphi,
                    chi2rz_t   chi2rz  ,
                    bendChi2_t bendChi2,
                    hit_t      hit     ,
                    trackMVA_t trackMVA,
                    extraMVA_t extraMVA,
                    valid_t    valid   ); 

void unpack_L1T_track(ap_uint<kTrackWordSize> tk,
                      rinv_t     &rinv    ,
                      tkphi_t    &tkphi   ,
                      tanlam_t   &tanlam  ,
                      tkz0_t     &tkz0    ,
                      tkd0_t     &tkd0    ,
                      chi2rphi_t &chi2rphi,
                      chi2rz_t   &chi2rz  ,
                      bendChi2_t &bendChi2,
                      hit_t      &hit     ,
                      trackMVA_t &trackMVA,
                      extraMVA_t &extraMVA,
                      valid_t    &valid   );

void pack_pf_track(ap_uint<64> &tk,
                   pt_t     pf_pt   ,
                   pt_t     pf_pterr,
                   etaphi_t pf_eta  ,
                   etaphi_t pf_phi  ,
                   z0_t     pf_z0   ,
                   bool     pf_TightQuality);

void unpack_pf_track(ap_uint<64> tk,
                     pt_t     &pf_pt   ,
                     pt_t     &pf_pterr,
                     etaphi_t &pf_eta  ,
                     etaphi_t &pf_phi  ,
                     z0_t     &pf_z0   ,
                     bool     &pf_TightQuality);

template<class in_t, class out_t> void bit_copy(in_t in, out_t &out, int offset=0);

template<class pt_inv_T, class pt_T> void init_pt_inv_table(pt_T table_out[(1<<PT_INV_TAB_SIZE)]);
template<class pt_inv_T, class pt_T> void convert_pt(pt_inv_T inv, pt_T &pt);

template<class eta_T> void init_eta_table(eta_T table_out[(1<<ETA_TAB_SIZE)]);
template<class tanlam_T, class eta_T> void convert_eta(tanlam_T tanlam, eta_T &eta);

template<class phi_T> void init_dphi_table(phi_T table_out[(1<DPHI_TAB_SIZE)]);
template<class pt_inv_T, class phi_T> void convert_dphi(pt_inv_T inv, phi_T &dphi);

//void track_propagate(pt_t pt, etaphi_t eta, etaphi_t phi, z0_t z0,etaphi_t &eta_calo, etaphi_t &phi_calo);
void reso_calo(pt_t pt, etaphi_t eta_calo, pt_t& err);

void propagate_tanlam(tkz0_t z0, tanlam_t tanlam, tanlam_t &tanlam_at_det);

inline float tanlam_to_eta(float tanlam){return -log(tan((M_PI/2 - atan(tanlam))/2));}


#endif

/*
Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041


ROADMAP (triple_pendulum.sd)

Bodies        Inb
No  Name      body Joint type  Coords q
--- --------- ---- ----------- ----------------
 -1 $ground                                    
  0 body1      -1  Pin           0             
  1 body2       0  Pin           1             
  2 body3       1  Pin           2             

*/
#include <math.h>
#include <stdio.h>

typedef struct {
    int ground_,nbod_,ndof_,ncons_,nloop_,nldof_,nloopc_,nball_,nlball_,npres_,
      nuser_;
    int jtype_[3],inb_[3],outb_[3],njntdof_[3],njntc_[3],njntp_[3],firstq_[3],
      ballq_[3],firstm_[3],firstp_[3];
    int trans_[3];
} sdgtopo_t;
#define ground (sdgtopo.ground_)
#define nbod (sdgtopo.nbod_)
#define ndof (sdgtopo.ndof_)
#define ncons (sdgtopo.ncons_)
#define nloop (sdgtopo.nloop_)
#define nldof (sdgtopo.nldof_)
#define nloopc (sdgtopo.nloopc_)
#define nball (sdgtopo.nball_)
#define nlball (sdgtopo.nlball_)
#define npres (sdgtopo.npres_)
#define nuser (sdgtopo.nuser_)
#define jtype (sdgtopo.jtype_)
#define inb (sdgtopo.inb_)
#define outb (sdgtopo.outb_)
#define njntdof (sdgtopo.njntdof_)
#define njntc (sdgtopo.njntc_)
#define njntp (sdgtopo.njntp_)
#define firstq (sdgtopo.firstq_)
#define ballq (sdgtopo.ballq_)
#define firstm (sdgtopo.firstm_)
#define firstp (sdgtopo.firstp_)
#define trans (sdgtopo.trans_)

typedef struct {
    double grav_[3],mk_[3],ik_[3][3][3],pin_[3][3];
    double rk_[3][3],ri_[3][3],pres_[3],stabvel_,stabpos_;
    int mfrcflg_,roustate_,vpkflg_,inerflg_,mmflg_,mmlduflg_,wwflg_,ltauflg_,
      fs0flg_,ii_,mmap_[3];
    int gravq_[3],mkq_[3],ikq_[3][3][3],pinq_[3][3],rkq_[3][3],riq_[3][3],presq_
      [3],stabvelq_,stabposq_;
    double mtot_,psmkg_,rhead_[3][3],rcom_[3][3],mkrcomt_[3][3][3],psikg_[3][3],
      psrcomg_[3],psrkg_[3],psrig_[3],psmk_[3],psik_[3][3][3],psrcom_[3][3],
      psrk_[3][3],psri_[3][3];
} sdginput_t;
#define grav (sdginput.grav_)
#define mk (sdginput.mk_)
#define ik (sdginput.ik_)
#define pin (sdginput.pin_)
#define rk (sdginput.rk_)
#define ri (sdginput.ri_)
#define pres (sdginput.pres_)
#define stabvel (sdginput.stabvel_)
#define stabpos (sdginput.stabpos_)
#define rhead (sdginput.rhead_)
#define rcom (sdginput.rcom_)
#define psrcomg (sdginput.psrcomg_)
#define psrcom (sdginput.psrcom_)
#define mkrcomt (sdginput.mkrcomt_)
#define psmk (sdginput.psmk_)
#define psik (sdginput.psik_)
#define psrk (sdginput.psrk_)
#define psri (sdginput.psri_)
#define psmkg (sdginput.psmkg_)
#define psikg (sdginput.psikg_)
#define psrkg (sdginput.psrkg_)
#define psrig (sdginput.psrig_)
#define mtot (sdginput.mtot_)
#define mfrcflg (sdginput.mfrcflg_)
#define roustate (sdginput.roustate_)
#define vpkflg (sdginput.vpkflg_)
#define inerflg (sdginput.inerflg_)
#define mmflg (sdginput.mmflg_)
#define mmlduflg (sdginput.mmlduflg_)
#define wwflg (sdginput.wwflg_)
#define ltauflg (sdginput.ltauflg_)
#define fs0flg (sdginput.fs0flg_)
#define ii (sdginput.ii_)
#define mmap (sdginput.mmap_)
#define gravq (sdginput.gravq_)
#define mkq (sdginput.mkq_)
#define ikq (sdginput.ikq_)
#define pinq (sdginput.pinq_)
#define rkq (sdginput.rkq_)
#define riq (sdginput.riq_)
#define presq (sdginput.presq_)
#define stabvelq (sdginput.stabvelq_)
#define stabposq (sdginput.stabposq_)

typedef struct {
    double curtim_,q_[3],qn_[3],u_[3],cnk_[3][3][3],cnb_[3][3][3];
    double rnk_[3][3],vnk_[3][3],wk_[3][3],rnb_[3][3],vnb_[3][3],wb_[3][3],
      wbrcom_[3][3],com_[3],rnkg_[3];
    double Cik_[3][3][3],rikt_[3][3][3],Iko_[3][3][3],mkrk_[3][3][3],Cib_[3][3][
      3];
    double Wkk_[3][3],Vkk_[3][3],dik_[3][3],rpp_[3][3],rpk_[3][3],rik_[3][3],
      rik2_[3][3];
    double rpri_[3][3],Wik_[3][3],Vik_[3][3],Wirk_[3][3],rkWkk_[3][3],Wkrpk_[3][
      3],VikWkr_[3][3];
    double perr_[1],verr_[1],aerr_[1],mult_[1],ufk_[3][3],utk_[3][3],mfk_[3][3],
      mtk_[3][3];
    double utau_[3],mtau_[3],uacc_[3],uvel_[3],upos_[3];
    double s0_,c0_,s1_,c1_,s2_,c2_;
} sdgstate_t;
#define curtim (sdgstate.curtim_)
#define q (sdgstate.q_)
#define qn (sdgstate.qn_)
#define u (sdgstate.u_)
#define cnk (sdgstate.cnk_)
#define cnb (sdgstate.cnb_)
#define rnkg (sdgstate.rnkg_)
#define rnk (sdgstate.rnk_)
#define rnb (sdgstate.rnb_)
#define vnk (sdgstate.vnk_)
#define vnb (sdgstate.vnb_)
#define wk (sdgstate.wk_)
#define wb (sdgstate.wb_)
#define com (sdgstate.com_)
#define Cik (sdgstate.Cik_)
#define Cib (sdgstate.Cib_)
#define rikt (sdgstate.rikt_)
#define Iko (sdgstate.Iko_)
#define mkrk (sdgstate.mkrk_)
#define Wkk (sdgstate.Wkk_)
#define Vkk (sdgstate.Vkk_)
#define dik (sdgstate.dik_)
#define rpp (sdgstate.rpp_)
#define rpk (sdgstate.rpk_)
#define rik (sdgstate.rik_)
#define rik2 (sdgstate.rik2_)
#define rpri (sdgstate.rpri_)
#define Wik (sdgstate.Wik_)
#define Vik (sdgstate.Vik_)
#define Wirk (sdgstate.Wirk_)
#define rkWkk (sdgstate.rkWkk_)
#define Wkrpk (sdgstate.Wkrpk_)
#define VikWkr (sdgstate.VikWkr_)
#define wbrcom (sdgstate.wbrcom_)
#define perr (sdgstate.perr_)
#define verr (sdgstate.verr_)
#define aerr (sdgstate.aerr_)
#define mult (sdgstate.mult_)
#define ufk (sdgstate.ufk_)
#define utk (sdgstate.utk_)
#define utau (sdgstate.utau_)
#define mfk (sdgstate.mfk_)
#define mtk (sdgstate.mtk_)
#define mtau (sdgstate.mtau_)
#define uacc (sdgstate.uacc_)
#define uvel (sdgstate.uvel_)
#define upos (sdgstate.upos_)
#define s0 (sdgstate.s0_)
#define c0 (sdgstate.c0_)
#define s1 (sdgstate.s1_)
#define c1 (sdgstate.c1_)
#define s2 (sdgstate.s2_)
#define c2 (sdgstate.c2_)

typedef struct {
    double fs0_[3],qdot_[3],Otk_[3][3],Atk_[3][3],AiOiWi_[3][3],Fstar_[3][3];
    double Tstar_[3][3],Fstark_[3][3],Tstark_[3][3],IkWk_[3][3],WkIkWk_[3][3],
      gk_[3][3],IkbWk_[3][3],WkIkbWk_[3][3];
    double w0w0_[3],w1w1_[3],w2w2_[3],w0w1_[3],w0w2_[3],w1w2_[3];
    double w00w11_[3],w00w22_[3],w11w22_[3],ww_[1][1],qraux_[1];
    double mm_[3][3],mlo_[3][3],mdi_[3],IkWpk_[3][3][3],works_[3],workss_[3][3];
    double Wpk_[3][3][3],Vpk_[3][3][3],VWri_[3][3][3];
    int wmap_[1],multmap_[1],jpvt_[1],wsiz_,wrank_;
} sdglhs_t;
#define qdot (sdglhs.qdot_)
#define Otk (sdglhs.Otk_)
#define Atk (sdglhs.Atk_)
#define AiOiWi (sdglhs.AiOiWi_)
#define Fstar (sdglhs.Fstar_)
#define Tstar (sdglhs.Tstar_)
#define fs0 (sdglhs.fs0_)
#define Fstark (sdglhs.Fstark_)
#define Tstark (sdglhs.Tstark_)
#define IkWk (sdglhs.IkWk_)
#define IkbWk (sdglhs.IkbWk_)
#define WkIkWk (sdglhs.WkIkWk_)
#define WkIkbWk (sdglhs.WkIkbWk_)
#define gk (sdglhs.gk_)
#define w0w0 (sdglhs.w0w0_)
#define w1w1 (sdglhs.w1w1_)
#define w2w2 (sdglhs.w2w2_)
#define w0w1 (sdglhs.w0w1_)
#define w0w2 (sdglhs.w0w2_)
#define w1w2 (sdglhs.w1w2_)
#define w00w11 (sdglhs.w00w11_)
#define w00w22 (sdglhs.w00w22_)
#define w11w22 (sdglhs.w11w22_)
#define ww (sdglhs.ww_)
#define qraux (sdglhs.qraux_)
#define mm (sdglhs.mm_)
#define mlo (sdglhs.mlo_)
#define mdi (sdglhs.mdi_)
#define IkWpk (sdglhs.IkWpk_)
#define works (sdglhs.works_)
#define workss (sdglhs.workss_)
#define Wpk (sdglhs.Wpk_)
#define Vpk (sdglhs.Vpk_)
#define VWri (sdglhs.VWri_)
#define wmap (sdglhs.wmap_)
#define multmap (sdglhs.multmap_)
#define jpvt (sdglhs.jpvt_)
#define wsiz (sdglhs.wsiz_)
#define wrank (sdglhs.wrank_)

typedef struct {
    double fs_[3],udot_[3],tauc_[3],dyad_[3][3][3],fc_[3][3],tc_[3][3];
    double ank_[3][3],onk_[3][3],Onkb_[3][3],AOnkri_[3][3],Ankb_[3][3],AnkAtk_[3
      ][3],anb_[3][3],onb_[3][3],dyrcom_[3][3];
    double ffk_[3][3],ttk_[3][3],fccikt_[3][3],ffkb_[3][3],ttkb_[3][3];
} sdgrhs_t;
#define fs (sdgrhs.fs_)
#define udot (sdgrhs.udot_)
#define ank (sdgrhs.ank_)
#define anb (sdgrhs.anb_)
#define onk (sdgrhs.onk_)
#define onb (sdgrhs.onb_)
#define Onkb (sdgrhs.Onkb_)
#define AOnkri (sdgrhs.AOnkri_)
#define Ankb (sdgrhs.Ankb_)
#define AnkAtk (sdgrhs.AnkAtk_)
#define dyrcom (sdgrhs.dyrcom_)
#define ffk (sdgrhs.ffk_)
#define ttk (sdgrhs.ttk_)
#define fccikt (sdgrhs.fccikt_)
#define ffkb (sdgrhs.ffkb_)
#define ttkb (sdgrhs.ttkb_)
#define dyad (sdgrhs.dyad_)
#define fc (sdgrhs.fc_)
#define tc (sdgrhs.tc_)
#define tauc (sdgrhs.tauc_)

typedef struct {
    double temp_[3000],tmat1_[3][3],tmat2_[3][3],tvec1_[3],tvec2_[3],tvec3_[3],
      tvec4_[3],tvec5_[3];
    double tsc1_,tsc2_,tsc3_;
} sdgtemp_t;
#define temp (sdgtemp.temp_)
#define tmat1 (sdgtemp.tmat1_)
#define tmat2 (sdgtemp.tmat2_)
#define tvec1 (sdgtemp.tvec1_)
#define tvec2 (sdgtemp.tvec2_)
#define tvec3 (sdgtemp.tvec3_)
#define tvec4 (sdgtemp.tvec4_)
#define tvec5 (sdgtemp.tvec5_)
#define tsc1 (sdgtemp.tsc1_)
#define tsc2 (sdgtemp.tsc2_)
#define tsc3 (sdgtemp.tsc3_)

sdgtopo_t sdgtopo = {
/*  Topological information
*/
    /* ground */ 1,
    /* nbod */ 3,
    /* ndof */ 3,
    /* ncons */ 0,
    /* nloop */ 0,
    /* nldof */ 0,
    /* nloopc */ 0,
    /* nball */ 0,
    /* nlball */ 0,
    /* npres */ 0,
    /* nuser */ 0,
    /* jtype[0] */ 1,
    /* jtype[1] */ 1,
    /* jtype[2] */ 1,
    /* inb[0] */ -1,
    /* inb[1] */ 0,
    /* inb[2] */ 1,
    /* outb[0] */ 0,
    /* outb[1] */ 1,
    /* outb[2] */ 2,
    /* njntdof[0] */ 1,
    /* njntdof[1] */ 1,
    /* njntdof[2] */ 1,
    /* njntc[0] */ 0,
    /* njntc[1] */ 0,
    /* njntc[2] */ 0,
    /* njntp[0] */ 0,
    /* njntp[1] */ 0,
    /* njntp[2] */ 0,
    /* firstq[0] */ 0,
    /* firstq[1] */ 1,
    /* firstq[2] */ 2,
    /* ballq[0] */ -104,
    /* ballq[1] */ -104,
    /* ballq[2] */ -104,
    /* firstm[0] */ -1,
    /* firstm[1] */ -1,
    /* firstm[2] */ -1,
    /* firstp[0] */ -1,
    /* firstp[1] */ -1,
    /* firstp[2] */ -1,
    /* trans[0] */ 0,
    /* trans[1] */ 0,
    /* trans[2] */ 0,
};
sdginput_t sdginput = {
/* Model parameters from the input file */

/* gravity */
    /* grav[0] */ 0.,
    /* grav[1] */ -9.8,
    /* grav[2] */ 0.,

/* mass */
    /* mk[0] */ 10.,
    /* mk[1] */ 5.,
    /* mk[2] */ 2.,

/* inertia */
    /* ik[0][0][0] */ 0.,
    /* ik[0][0][1] */ 0.,
    /* ik[0][0][2] */ 0.,
    /* ik[0][1][0] */ 0.,
    /* ik[0][1][1] */ 0.,
    /* ik[0][1][2] */ 0.,
    /* ik[0][2][0] */ 0.,
    /* ik[0][2][1] */ 0.,
    /* ik[0][2][2] */ 0.,
    /* ik[1][0][0] */ 0.,
    /* ik[1][0][1] */ 0.,
    /* ik[1][0][2] */ 0.,
    /* ik[1][1][0] */ 0.,
    /* ik[1][1][1] */ 0.,
    /* ik[1][1][2] */ 0.,
    /* ik[1][2][0] */ 0.,
    /* ik[1][2][1] */ 0.,
    /* ik[1][2][2] */ 0.,
    /* ik[2][0][0] */ 0.,
    /* ik[2][0][1] */ 0.,
    /* ik[2][0][2] */ 0.,
    /* ik[2][1][0] */ 0.,
    /* ik[2][1][1] */ 0.,
    /* ik[2][1][2] */ 0.,
    /* ik[2][2][0] */ 0.,
    /* ik[2][2][1] */ 0.,
    /* ik[2][2][2] */ 0.,

/* tree hinge axis vectors */
    /* pin[0][0] */ 0.,
    /* pin[0][1] */ 0.,
    /* pin[0][2] */ 1.,
    /* pin[1][0] */ 0.,
    /* pin[1][1] */ 0.,
    /* pin[1][2] */ 1.,
    /* pin[2][0] */ 0.,
    /* pin[2][1] */ 0.,
    /* pin[2][2] */ 1.,

/* tree bodytojoint vectors */
    /* rk[0][0] */ -.5,
    /* rk[0][1] */ 0.,
    /* rk[0][2] */ 0.,
    /* rk[1][0] */ -.5,
    /* rk[1][1] */ 0.,
    /* rk[1][2] */ 0.,
    /* rk[2][0] */ -.25,
    /* rk[2][1] */ 0.,
    /* rk[2][2] */ 0.,

/* tree inbtojoint vectors */
    /* ri[0][0] */ 0.,
    /* ri[0][1] */ 0.,
    /* ri[0][2] */ 0.,
    /* ri[1][0] */ .5,
    /* ri[1][1] */ 0.,
    /* ri[1][2] */ 0.,
    /* ri[2][0] */ .5,
    /* ri[2][1] */ 0.,
    /* ri[2][2] */ 0.,

/* tree prescribed motion */
    /* pres[0] */ 0.,
    /* pres[1] */ 0.,
    /* pres[2] */ 0.,

/* stabilization parameters */
    /* stabvel */ 0.,
    /* stabpos */ 0.,

/* miscellaneous */
    /* mfrcflg */ 0,
    /* roustate */ 0,
    /* vpkflg */ 0,
    /* inerflg */ 0,
    /* mmflg */ 0,
    /* mmlduflg */ 0,
    /* wwflg */ 0,
    /* ltauflg */ 0,
    /* fs0flg */ 0,
    /* ii */ 0,
    /* mmap[0] */ 0,
    /* mmap[1] */ 1,
    /* mmap[2] */ 2,

/* Which parameters were "?" (1) or "<nominal>?" (3) */
    /* gravq[0] */ 3,
    /* gravq[1] */ 3,
    /* gravq[2] */ 3,
    /* mkq[0] */ 3,
    /* mkq[1] */ 3,
    /* mkq[2] */ 3,
    /* ikq[0][0][0] */ 3,
    /* ikq[0][0][1] */ 0,
    /* ikq[0][0][2] */ 0,
    /* ikq[0][1][0] */ 0,
    /* ikq[0][1][1] */ 3,
    /* ikq[0][1][2] */ 0,
    /* ikq[0][2][0] */ 0,
    /* ikq[0][2][1] */ 0,
    /* ikq[0][2][2] */ 3,
    /* ikq[1][0][0] */ 3,
    /* ikq[1][0][1] */ 0,
    /* ikq[1][0][2] */ 0,
    /* ikq[1][1][0] */ 0,
    /* ikq[1][1][1] */ 3,
    /* ikq[1][1][2] */ 0,
    /* ikq[1][2][0] */ 0,
    /* ikq[1][2][1] */ 0,
    /* ikq[1][2][2] */ 3,
    /* ikq[2][0][0] */ 3,
    /* ikq[2][0][1] */ 0,
    /* ikq[2][0][2] */ 0,
    /* ikq[2][1][0] */ 0,
    /* ikq[2][1][1] */ 3,
    /* ikq[2][1][2] */ 0,
    /* ikq[2][2][0] */ 0,
    /* ikq[2][2][1] */ 0,
    /* ikq[2][2][2] */ 3,
    /* pinq[0][0] */ 0,
    /* pinq[0][1] */ 0,
    /* pinq[0][2] */ 0,
    /* pinq[1][0] */ 0,
    /* pinq[1][1] */ 0,
    /* pinq[1][2] */ 0,
    /* pinq[2][0] */ 0,
    /* pinq[2][1] */ 0,
    /* pinq[2][2] */ 0,
    /* rkq[0][0] */ 3,
    /* rkq[0][1] */ 3,
    /* rkq[0][2] */ 3,
    /* rkq[1][0] */ 3,
    /* rkq[1][1] */ 3,
    /* rkq[1][2] */ 3,
    /* rkq[2][0] */ 3,
    /* rkq[2][1] */ 3,
    /* rkq[2][2] */ 3,
    /* riq[0][0] */ 0,
    /* riq[0][1] */ 0,
    /* riq[0][2] */ 0,
    /* riq[1][0] */ 3,
    /* riq[1][1] */ 3,
    /* riq[1][2] */ 3,
    /* riq[2][0] */ 3,
    /* riq[2][1] */ 3,
    /* riq[2][2] */ 3,
    /* presq[0] */ 0,
    /* presq[1] */ 0,
    /* presq[2] */ 0,
    /* stabvelq */ 3,
    /* stabposq */ 3,

/* End of values from input file */

};
sdgstate_t sdgstate;
sdglhs_t sdglhs;
sdgrhs_t sdgrhs;
sdgtemp_t sdgtemp;


void sdinit(void)
{
/*
Initialization routine


 This routine must be called before the first call to sdstate(), after
 supplying values for any `?' parameters in the input.
*/
    double sumsq,norminv;
    int i,j,k;


/* Check that all `?' parameters have been assigned values */

    for (k = 0; k < 3; k++) {
        if (gravq[k] == 1) {
            sdseterr(7,25);
        }
    }
    for (k = 0; k < 3; k++) {
        if (mkq[k] == 1) {
            sdseterr(7,26);
        }
        for (i = 0; i < 3; i++) {
            if (rkq[k][i] == 1) {
                sdseterr(7,29);
            }
            if (riq[k][i] == 1) {
                sdseterr(7,30);
            }
            for (j = 0; j < 3; j++) {
                if (ikq[k][i][j] == 1) {
                    sdseterr(7,27);
                }
            }
        }
    }
    for (k = 0; k < 3; k++) {
        for (i = 0; i < 3; i++) {
            if (pinq[k][i] == 1) {
                sdseterr(7,28);
            }
        }
    }

/* Normalize pin vectors if necessary */


/* Zero out Vpk and Wpk */

    for (i = 0; i < 3; i++) {
        for (j = i; j <= 2; j++) {
            for (k = 0; k < 3; k++) {
                Vpk[i][j][k] = 0.;
                Wpk[i][j][k] = 0.;
            }
        }
    }

/* Compute pseudobody-related constants */

    rcom[0][0] = 0.;
    rcom[0][1] = 0.;
    rcom[0][2] = 0.;
    rcom[1][0] = 0.;
    rcom[1][1] = 0.;
    rcom[1][2] = 0.;
    rcom[2][0] = 0.;
    rcom[2][1] = 0.;
    rcom[2][2] = 0.;
    dik[1][0] = (ri[1][0]-rk[0][0]);
    dik[1][1] = (ri[1][1]-rk[0][1]);
    dik[1][2] = (ri[1][2]-rk[0][2]);
    dik[2][0] = (ri[2][0]-rk[1][0]);
    dik[2][1] = (ri[2][1]-rk[1][1]);
    dik[2][2] = (ri[2][2]-rk[1][2]);

/* Compute mass properties-related constants */

    mtot = (mk[2]+(mk[0]+mk[1]));
    mkrk[0][0][1] = -(mk[0]*rk[0][2]);
    mkrk[0][0][2] = (mk[0]*rk[0][1]);
    mkrk[0][1][0] = (mk[0]*rk[0][2]);
    mkrk[0][1][2] = -(mk[0]*rk[0][0]);
    mkrk[0][2][0] = -(mk[0]*rk[0][1]);
    mkrk[0][2][1] = (mk[0]*rk[0][0]);
    mkrk[1][0][1] = -(mk[1]*rk[1][2]);
    mkrk[1][0][2] = (mk[1]*rk[1][1]);
    mkrk[1][1][0] = (mk[1]*rk[1][2]);
    mkrk[1][1][2] = -(mk[1]*rk[1][0]);
    mkrk[1][2][0] = -(mk[1]*rk[1][1]);
    mkrk[1][2][1] = (mk[1]*rk[1][0]);
    mkrk[2][0][1] = -(mk[2]*rk[2][2]);
    mkrk[2][0][2] = (mk[2]*rk[2][1]);
    mkrk[2][1][0] = (mk[2]*rk[2][2]);
    mkrk[2][1][2] = -(mk[2]*rk[2][0]);
    mkrk[2][2][0] = -(mk[2]*rk[2][1]);
    mkrk[2][2][1] = (mk[2]*rk[2][0]);
    Iko[0][0][0] = (ik[0][0][0]-((mkrk[0][0][1]*rk[0][2])-(mkrk[0][0][2]*
      rk[0][1])));
    Iko[0][0][1] = -(mkrk[0][0][2]*rk[0][0]);
    Iko[0][0][2] = (mkrk[0][0][1]*rk[0][0]);
    Iko[0][1][0] = (mkrk[0][1][2]*rk[0][1]);
    Iko[0][1][1] = (ik[0][1][1]-((mkrk[0][1][2]*rk[0][0])-(mkrk[0][1][0]*
      rk[0][2])));
    Iko[0][1][2] = -(mkrk[0][1][0]*rk[0][1]);
    Iko[0][2][0] = -(mkrk[0][2][1]*rk[0][2]);
    Iko[0][2][1] = (mkrk[0][2][0]*rk[0][2]);
    Iko[0][2][2] = (ik[0][2][2]-((mkrk[0][2][0]*rk[0][1])-(mkrk[0][2][1]*
      rk[0][0])));
    Iko[1][0][0] = (ik[1][0][0]-((mkrk[1][0][1]*rk[1][2])-(mkrk[1][0][2]*
      rk[1][1])));
    Iko[1][0][1] = -(mkrk[1][0][2]*rk[1][0]);
    Iko[1][0][2] = (mkrk[1][0][1]*rk[1][0]);
    Iko[1][1][0] = (mkrk[1][1][2]*rk[1][1]);
    Iko[1][1][1] = (ik[1][1][1]-((mkrk[1][1][2]*rk[1][0])-(mkrk[1][1][0]*
      rk[1][2])));
    Iko[1][1][2] = -(mkrk[1][1][0]*rk[1][1]);
    Iko[1][2][0] = -(mkrk[1][2][1]*rk[1][2]);
    Iko[1][2][1] = (mkrk[1][2][0]*rk[1][2]);
    Iko[1][2][2] = (ik[1][2][2]-((mkrk[1][2][0]*rk[1][1])-(mkrk[1][2][1]*
      rk[1][0])));
    Iko[2][0][0] = (ik[2][0][0]-((mkrk[2][0][1]*rk[2][2])-(mkrk[2][0][2]*
      rk[2][1])));
    Iko[2][0][1] = -(mkrk[2][0][2]*rk[2][0]);
    Iko[2][0][2] = (mkrk[2][0][1]*rk[2][0]);
    Iko[2][1][0] = (mkrk[2][1][2]*rk[2][1]);
    Iko[2][1][1] = (ik[2][1][1]-((mkrk[2][1][2]*rk[2][0])-(mkrk[2][1][0]*
      rk[2][2])));
    Iko[2][1][2] = -(mkrk[2][1][0]*rk[2][1]);
    Iko[2][2][0] = -(mkrk[2][2][1]*rk[2][2]);
    Iko[2][2][1] = (mkrk[2][2][0]*rk[2][2]);
    Iko[2][2][2] = (ik[2][2][2]-((mkrk[2][2][0]*rk[2][1])-(mkrk[2][2][1]*
      rk[2][0])));
    sdserialno(&i);
    if (i != 30123) {
        sdseterr(7,41);
    }
    roustate = 1;
}

/* Convert state to form using 1-2-3 Euler angles for ball joints. */

void sdst2ang(double st[3],
    double stang[3])
{
    int i;

    for (i = 0; i < 3; i++) {
        stang[i] = st[i];
    }
}

/* Convert 1-2-3 form of state back to Euler parameters for ball joints. */

void sdang2st(double stang[3],
    double st[3])
{
    int i;

    for (i = 0; i < 3; i++) {
        st[i] = stang[i];
    }
}

/* Normalize Euler parameters in state. */

void sdnrmsterr(double st[3],
    double normst[3],
    int routine)
{
    int i;

    for (i = 0; i < 3; i++) {
        normst[i] = st[i];
    }
}

void sdnormst(double st[3],
    double normst[3])
{

    sdnrmsterr(st,normst,0);
}

void sdstate(double timein,
    double qin[3],
    double uin[3])
{
/*
Compute kinematic information and store it in sdgstate.

Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041
*/
    int i,j,qchg,uchg;

    if ((roustate != 1) && (roustate != 2) && (roustate != 3)) {
        sdseterr(8,22);
        return;
    }
    if (roustate == 1) {
        for (i = 0; i < 3; i++) {
            if (presq[i] == 1) {
                sdseterr(8,31);
            }
        }
    }
/*
See if time or any qs have changed since last call
*/
    if ((roustate != 1) && (timein == curtim)) {
        qchg = 0;
        for (i = 0; i < 3; i++) {
            if (qin[i] != q[i]) {
                qchg = 1;
                break;
            }
        }
    } else {
        qchg = 1;
    }
/*
If time and qs are unchanged, check us
*/
    if (qchg == 0) {
        uchg = 0;
        for (i = 0; i < 3; i++) {
            if (uin[i] != u[i]) {
                uchg = 1;
                break;
            }
        }
    } else {
        uchg = 1;
    }
    curtim = timein;
    roustate = 2;
    if (qchg == 0) {
        goto skipqs;
    }
/*
Position-related variables need to be computed
*/
    vpkflg = 0;
    mmflg = 0;
    mmlduflg = 0;
    wwflg = 0;
    for (i = 0; i < 3; i++) {
        q[i] = qin[i];
    }
/*
Compute sines and cosines of q
*/
    s0 = sin(q[0]);
    c0 = cos(q[0]);
    s1 = sin(q[1]);
    c1 = cos(q[1]);
    s2 = sin(q[2]);
    c2 = cos(q[2]);
/*
Compute across-axis direction cosines Cik
*/
/*
Compute across-joint direction cosines Cib
*/
/*
Compute gravity
*/
    gk[0][0] = ((grav[0]*c0)+(grav[1]*s0));
    gk[0][1] = ((grav[1]*c0)-(grav[0]*s0));
    gk[1][0] = ((gk[0][0]*c1)+(gk[0][1]*s1));
    gk[1][1] = ((gk[0][1]*c1)-(gk[0][0]*s1));
    gk[2][0] = ((gk[1][0]*c2)+(gk[1][1]*s2));
    gk[2][1] = ((gk[1][1]*c2)-(gk[1][0]*s2));
/*
Compute cnk & cnb (direction cosines in N)
*/
    cnk[1][0][0] = ((c0*c1)-(s0*s1));
    cnk[1][0][1] = -((s0*c1)+(s1*c0));
    cnk[1][1][0] = ((s0*c1)+(s1*c0));
    cnk[1][1][1] = ((c0*c1)-(s0*s1));
    cnk[2][0][0] = ((cnk[1][0][0]*c2)+(cnk[1][0][1]*s2));
    cnk[2][0][1] = ((cnk[1][0][1]*c2)-(cnk[1][0][0]*s2));
    cnk[2][1][0] = ((cnk[1][1][0]*c2)+(cnk[1][1][1]*s2));
    cnk[2][1][1] = ((cnk[1][1][1]*c2)-(cnk[1][1][0]*s2));
    cnb[0][0][0] = c0;
    cnb[0][0][1] = -s0;
    cnb[0][0][2] = 0.;
    cnb[0][1][0] = s0;
    cnb[0][1][1] = c0;
    cnb[0][1][2] = 0.;
    cnb[0][2][0] = 0.;
    cnb[0][2][1] = 0.;
    cnb[0][2][2] = 1.;
    cnb[1][0][0] = cnk[1][0][0];
    cnb[1][0][1] = cnk[1][0][1];
    cnb[1][0][2] = 0.;
    cnb[1][1][0] = cnk[1][1][0];
    cnb[1][1][1] = cnk[1][1][1];
    cnb[1][1][2] = 0.;
    cnb[1][2][0] = 0.;
    cnb[1][2][1] = 0.;
    cnb[1][2][2] = 1.;
    cnb[2][0][0] = cnk[2][0][0];
    cnb[2][0][1] = cnk[2][0][1];
    cnb[2][0][2] = 0.;
    cnb[2][1][0] = cnk[2][1][0];
    cnb[2][1][1] = cnk[2][1][1];
    cnb[2][1][2] = 0.;
    cnb[2][2][0] = 0.;
    cnb[2][2][1] = 0.;
    cnb[2][2][2] = 1.;
/*
Compute q-related auxiliary variables
*/
    rik[1][0] = (((ri[1][0]*c1)+(ri[1][1]*s1))-rk[1][0]);
    rik[1][1] = (((ri[1][1]*c1)-(ri[1][0]*s1))-rk[1][1]);
    rik[1][2] = (ri[1][2]-rk[1][2]);
    rik[2][0] = (((ri[2][0]*c2)+(ri[2][1]*s2))-rk[2][0]);
    rik[2][1] = (((ri[2][1]*c2)-(ri[2][0]*s2))-rk[2][1]);
    rik[2][2] = (ri[2][2]-rk[2][2]);
/*
Compute rnk & rnb (mass center locations in N)
*/
    rnk[0][0] = ((rk[0][1]*s0)-(rk[0][0]*c0));
    rnk[0][1] = -((rk[0][0]*s0)+(rk[0][1]*c0));
    rnk[1][0] = ((rnk[0][0]+((ri[1][0]*c0)-(ri[1][1]*s0)))-((cnk[1][0][0]*
      rk[1][0])+(cnk[1][0][1]*rk[1][1])));
    rnk[1][1] = ((rnk[0][1]+((ri[1][0]*s0)+(ri[1][1]*c0)))-((cnk[1][1][0]*
      rk[1][0])+(cnk[1][1][1]*rk[1][1])));
    rnk[1][2] = ((ri[1][2]-rk[0][2])-rk[1][2]);
    rnk[2][0] = ((rnk[1][0]+((cnk[1][0][0]*ri[2][0])+(cnk[1][0][1]*ri[2][1])))-(
      (cnk[2][0][0]*rk[2][0])+(cnk[2][0][1]*rk[2][1])));
    rnk[2][1] = ((rnk[1][1]+((cnk[1][1][0]*ri[2][0])+(cnk[1][1][1]*ri[2][1])))-(
      (cnk[2][1][0]*rk[2][0])+(cnk[2][1][1]*rk[2][1])));
    rnk[2][2] = ((ri[2][2]+rnk[1][2])-rk[2][2]);
    rnb[0][0] = rnk[0][0];
    rnb[0][1] = rnk[0][1];
    rnb[0][2] = -rk[0][2];
    rnb[1][0] = rnk[1][0];
    rnb[1][1] = rnk[1][1];
    rnb[1][2] = rnk[1][2];
    rnb[2][0] = rnk[2][0];
    rnb[2][1] = rnk[2][1];
    rnb[2][2] = rnk[2][2];
/*
Compute com (system mass center location in N)
*/
    com[0] = ((1./mtot)*((mk[2]*rnk[2][0])+((mk[0]*rnk[0][0])+(mk[1]*rnk[1][0]))
      ));
    com[1] = ((1./mtot)*((mk[2]*rnk[2][1])+((mk[0]*rnk[0][1])+(mk[1]*rnk[1][1]))
      ));
    com[2] = ((1./mtot)*((mk[2]*rnk[2][2])+((mk[1]*rnk[1][2])-(mk[0]*rk[0][2])))
      );
    skipqs: ;
    if (uchg == 0) {
        goto skipus;
    }
/*
Velocity-related variables need to be computed
*/
    inerflg = 0;
    for (i = 0; i < 3; i++) {
        u[i] = uin[i];
    }
/*
Compute u-related auxiliary variables
*/
/*
Compute wk & wb (angular velocities)
*/
    wk[1][2] = (u[0]+u[1]);
    wk[2][2] = (u[2]+wk[1][2]);
    wb[0][0] = 0.;
    wb[0][1] = 0.;
    wb[0][2] = u[0];
    wb[1][0] = 0.;
    wb[1][1] = 0.;
    wb[1][2] = wk[1][2];
    wb[2][0] = 0.;
    wb[2][1] = 0.;
    wb[2][2] = wk[2][2];
/*
Compute auxiliary variables involving wk
*/
    Wirk[1][0] = -(ri[1][1]*u[0]);
    Wirk[1][1] = (ri[1][0]*u[0]);
    Wirk[2][0] = -(ri[2][1]*wk[1][2]);
    Wirk[2][1] = (ri[2][0]*wk[1][2]);
    Wkrpk[0][0] = (rk[0][1]*u[0]);
    Wkrpk[0][1] = -(rk[0][0]*u[0]);
    Wkrpk[1][0] = (rk[1][1]*wk[1][2]);
    Wkrpk[1][1] = -(rk[1][0]*wk[1][2]);
    Wkrpk[2][0] = (rk[2][1]*wk[2][2]);
    Wkrpk[2][1] = -(rk[2][0]*wk[2][2]);
    IkWk[0][2] = (ik[0][2][2]*u[0]);
    IkWk[1][2] = (ik[1][2][2]*wk[1][2]);
    IkWk[2][2] = (ik[2][2][2]*wk[2][2]);
/*
Compute temporaries for use in SDRHS
*/
    w2w2[0] = (u[0]*u[0]);
    w2w2[1] = (wk[1][2]*wk[1][2]);
    w2w2[2] = (wk[2][2]*wk[2][2]);
/*
Compute vnk & vnb (mass center linear velocities in N)
*/
    vnk[0][0] = ((Wkrpk[0][0]*c0)-(Wkrpk[0][1]*s0));
    vnk[0][1] = ((Wkrpk[0][0]*s0)+(Wkrpk[0][1]*c0));
    vnk[1][0] = ((vnk[0][0]+((Wirk[1][0]*c0)-(Wirk[1][1]*s0)))+((cnk[1][0][0]*
      Wkrpk[1][0])+(cnk[1][0][1]*Wkrpk[1][1])));
    vnk[1][1] = ((vnk[0][1]+((Wirk[1][0]*s0)+(Wirk[1][1]*c0)))+((cnk[1][1][0]*
      Wkrpk[1][0])+(cnk[1][1][1]*Wkrpk[1][1])));
    vnk[2][0] = ((vnk[1][0]+((cnk[1][0][0]*Wirk[2][0])+(cnk[1][0][1]*Wirk[2][1])
      ))+((cnk[2][0][0]*Wkrpk[2][0])+(cnk[2][0][1]*Wkrpk[2][1])));
    vnk[2][1] = ((vnk[1][1]+((cnk[1][1][0]*Wirk[2][0])+(cnk[1][1][1]*Wirk[2][1])
      ))+((cnk[2][1][0]*Wkrpk[2][0])+(cnk[2][1][1]*Wkrpk[2][1])));
    vnb[0][0] = vnk[0][0];
    vnb[0][1] = vnk[0][1];
    vnb[0][2] = 0.;
    vnb[1][0] = vnk[1][0];
    vnb[1][1] = vnk[1][1];
    vnb[1][2] = 0.;
    vnb[2][0] = vnk[2][0];
    vnb[2][1] = vnk[2][1];
    vnb[2][2] = 0.;
/*
Compute qdot (kinematical equations)
*/
    qdot[0] = u[0];
    qdot[1] = u[1];
    qdot[2] = u[2];
    skipus: ;
/*
Initialize applied forces and torques to zero
*/
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            ufk[i][j] = 0.;
            utk[i][j] = 0.;
        }
    }
    for (i = 0; i < 3; i++) {
        utau[i] = 0.;
    }
    ltauflg = 0;
    fs0flg = 0;
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   81 adds/subtracts/negates
                    104 multiplies
                      3 divides
                    139 assignments
*/
}

void sdqdot(double oqdot[3])
{
/*
Return position coordinate derivatives for tree joints.
*/
    int i;

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(63,23);
        return;
    }
    for (i = 0; i <= 2; i++) {
        oqdot[i] = qdot[i];
    }
}

void sdu2qdot(double uin[3],
    double oqdot[3])
{
/*
Convert velocities to qdots.
*/
    int i;

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(64,23);
        return;
    }
    for (i = 0; i <= 2; i++) {
        oqdot[i] = uin[i];
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdpsstate(double lqin[1])
{

    if (roustate != 2) {
        sdseterr(9,23);
        return;
    }
}

void sddovpk(void)
{

    if (vpkflg == 0) {
/*
Compute Wpk (partial angular velocities)
*/
        Wpk[0][0][2] = 1.;
        Wpk[0][1][2] = 1.;
        Wpk[0][2][2] = 1.;
        Wpk[1][1][2] = 1.;
        Wpk[1][2][2] = 1.;
        Wpk[2][2][2] = 1.;
/*
Compute Vpk (partial velocities)
*/
        Vpk[0][0][0] = rk[0][1];
        Vpk[0][0][1] = -rk[0][0];
        VWri[0][1][0] = (rk[0][1]-ri[1][1]);
        VWri[0][1][1] = (ri[1][0]-rk[0][0]);
        Vpk[0][1][0] = (rk[1][1]+((VWri[0][1][0]*c1)+(VWri[0][1][1]*s1)));
        Vpk[0][1][1] = (((VWri[0][1][1]*c1)-(VWri[0][1][0]*s1))-rk[1][0]);
        VWri[0][2][0] = (Vpk[0][1][0]-ri[2][1]);
        VWri[0][2][1] = (ri[2][0]+Vpk[0][1][1]);
        Vpk[0][2][0] = (rk[2][1]+((VWri[0][2][0]*c2)+(VWri[0][2][1]*s2)));
        Vpk[0][2][1] = (((VWri[0][2][1]*c2)-(VWri[0][2][0]*s2))-rk[2][0]);
        Vpk[1][1][0] = rk[1][1];
        Vpk[1][1][1] = -rk[1][0];
        VWri[1][2][0] = (rk[1][1]-ri[2][1]);
        VWri[1][2][1] = (ri[2][0]-rk[1][0]);
        Vpk[1][2][0] = (rk[2][1]+((VWri[1][2][0]*c2)+(VWri[1][2][1]*s2)));
        Vpk[1][2][1] = (((VWri[1][2][1]*c2)-(VWri[1][2][0]*s2))-rk[2][0]);
        Vpk[2][2][0] = rk[2][1];
        Vpk[2][2][1] = -rk[2][0];
        vpkflg = 1;
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   21 adds/subtracts/negates
                     12 multiplies
                      0 divides
                     24 assignments
*/
}

void sddoltau(void)
{

/*
Compute effect of loop hinge torques
*/
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      0 assignments
*/
}

void sddoiner(void)
{

/*
Compute inertial accelerations and related temps
*/
    if (inerflg == 0) {
/*
Compute Otk (inertial angular acceleration)
*/
/*
Compute Atk (inertial linear acceleration)
*/
        Atk[0][0] = -(u[0]*Wkrpk[0][1]);
        Atk[0][1] = (u[0]*Wkrpk[0][0]);
        AiOiWi[1][0] = (Atk[0][0]-(u[0]*Wirk[1][1]));
        AiOiWi[1][1] = (Atk[0][1]+(u[0]*Wirk[1][0]));
        Atk[1][0] = (((AiOiWi[1][0]*c1)+(AiOiWi[1][1]*s1))-(wk[1][2]*Wkrpk[1][1]
          ));
        Atk[1][1] = ((wk[1][2]*Wkrpk[1][0])+((AiOiWi[1][1]*c1)-(AiOiWi[1][0]*s1)
          ));
        AiOiWi[2][0] = (Atk[1][0]-(Wirk[2][1]*wk[1][2]));
        AiOiWi[2][1] = (Atk[1][1]+(Wirk[2][0]*wk[1][2]));
        Atk[2][0] = (((AiOiWi[2][0]*c2)+(AiOiWi[2][1]*s2))-(wk[2][2]*Wkrpk[2][1]
          ));
        Atk[2][1] = ((wk[2][2]*Wkrpk[2][0])+((AiOiWi[2][1]*c2)-(AiOiWi[2][0]*s2)
          ));
        inerflg = 1;
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   13 adds/subtracts/negates
                     18 multiplies
                      0 divides
                     10 assignments
*/
}

void sddofs0(void)
{

/*
Compute effect of all applied loads
*/
    if (fs0flg == 0) {
        sddoltau();
        sddoiner();
/*
Compute Fstar (forces)
*/
        Fstar[0][0] = ((mk[0]*(Atk[0][0]-gk[0][0]))-ufk[0][0]);
        Fstar[0][1] = ((mk[0]*(Atk[0][1]-gk[0][1]))-ufk[0][1]);
        Fstar[0][2] = -(ufk[0][2]+(grav[2]*mk[0]));
        Fstar[1][0] = ((mk[1]*(Atk[1][0]-gk[1][0]))-ufk[1][0]);
        Fstar[1][1] = ((mk[1]*(Atk[1][1]-gk[1][1]))-ufk[1][1]);
        Fstar[1][2] = -(ufk[1][2]+(grav[2]*mk[1]));
        Fstar[2][0] = ((mk[2]*(Atk[2][0]-gk[2][0]))-ufk[2][0]);
        Fstar[2][1] = ((mk[2]*(Atk[2][1]-gk[2][1]))-ufk[2][1]);
        Fstar[2][2] = -(ufk[2][2]+(grav[2]*mk[2]));
/*
Compute Tstar (torques)
*/
/*
Compute fs0 (RHS ignoring constraints)
*/
        sddovpk();
        fs0[0] = (utau[0]-((((Fstar[2][0]*Vpk[0][2][0])+(Fstar[2][1]*
          Vpk[0][2][1]))-utk[2][2])+((((Fstar[0][0]*rk[0][1])-(Fstar[0][1]*
          rk[0][0]))-utk[0][2])+(((Fstar[1][0]*Vpk[0][1][0])+(Fstar[1][1]*
          Vpk[0][1][1]))-utk[1][2]))));
        fs0[1] = (utau[1]-((((Fstar[1][0]*rk[1][1])-(Fstar[1][1]*rk[1][0]))-
          utk[1][2])+(((Fstar[2][0]*Vpk[1][2][0])+(Fstar[2][1]*Vpk[1][2][1]))-
          utk[2][2])));
        fs0[2] = (utau[2]-(((Fstar[2][0]*rk[2][1])-(Fstar[2][1]*rk[2][0]))-
          utk[2][2]));
        fs0flg = 1;
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   36 adds/subtracts/negates
                     21 multiplies
                      0 divides
                     12 assignments
*/
}

void sddomm(int routine)
{
    int dumroutine,errnum;
    int i;

    if (mmflg == 0) {
/*
Compute mass matrix (MM)
*/
        sddovpk();
        mm[0][0] = ((ik[2][2][2]+(mk[2]*((Vpk[0][2][0]*Vpk[0][2][0])+(
          Vpk[0][2][1]*Vpk[0][2][1]))))+((ik[0][2][2]+(mk[0]*((rk[0][0]*rk[0][0]
          )+(rk[0][1]*rk[0][1]))))+(ik[1][2][2]+(mk[1]*((Vpk[0][1][0]*
          Vpk[0][1][0])+(Vpk[0][1][1]*Vpk[0][1][1]))))));
        mm[0][1] = ((ik[1][2][2]+(mk[1]*((rk[1][1]*Vpk[0][1][0])-(rk[1][0]*
          Vpk[0][1][1]))))+(ik[2][2][2]+(mk[2]*((Vpk[0][2][0]*Vpk[1][2][0])+(
          Vpk[0][2][1]*Vpk[1][2][1])))));
        mm[0][2] = (ik[2][2][2]+(mk[2]*((rk[2][1]*Vpk[0][2][0])-(rk[2][0]*
          Vpk[0][2][1]))));
        mm[1][1] = ((ik[1][2][2]+(mk[1]*((rk[1][0]*rk[1][0])+(rk[1][1]*rk[1][1])
          )))+(ik[2][2][2]+(mk[2]*((Vpk[1][2][0]*Vpk[1][2][0])+(Vpk[1][2][1]*
          Vpk[1][2][1])))));
        mm[1][2] = (ik[2][2][2]+(mk[2]*((rk[2][1]*Vpk[1][2][0])-(rk[2][0]*
          Vpk[1][2][1]))));
        mm[2][2] = (ik[2][2][2]+(mk[2]*((rk[2][0]*rk[2][0])+(rk[2][1]*rk[2][1]))
          ));
/*
Check for singular mass matrix
*/
        for (i = 0; i < 3; i++) {
            if (mm[i][i] < 1e-13) {
                sdseterr(routine,47);
            }
        }
        sderror(&dumroutine,&errnum);
        if (errnum == 0) {
            mmflg = 1;
        }
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   24 adds/subtracts/negates
                     30 multiplies
                      0 divides
                      6 assignments
*/
}

void sddommldu(int routine)
{
    int i;
    int dumroutine,errnum;

    if (mmlduflg == 0) {
        sddomm(routine);
/*
Numerically decompose the mass matrix
*/
        sdldudcomp(3,3,mmap,1e-13,workss,works,mm,mlo,mdi);
/*
Check for singular mass matrix
*/
        for (i = 0; i < 3; i++) {
            if (mdi[i] <= 1e-13) {
                sdseterr(routine,47);
            }
        }
        sderror(&dumroutine,&errnum);
        if (errnum == 0) {
            mmlduflg = 1;
        }
    }
}

void sdlhs(int routine)
{
/* Compute all remaining state- and force-dependent quantities
*/

    roustate = 2;
    sddommldu(routine);
    sddofs0();
}

void sdmfrc(double imult[1])
{

}

void sdequivht(double tau[3])
{
/* Compute tree hinge torques to match effect of applied loads
*/
    double fstareq[3][3],tstareq[3][3];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(56,23);
        return;
    }
/*
Compute fstareq (forces)
*/
    fstareq[0][0] = -(ufk[0][0]+(gk[0][0]*mk[0]));
    fstareq[0][1] = -(ufk[0][1]+(gk[0][1]*mk[0]));
    fstareq[0][2] = -(ufk[0][2]+(grav[2]*mk[0]));
    fstareq[1][0] = -(ufk[1][0]+(gk[1][0]*mk[1]));
    fstareq[1][1] = -(ufk[1][1]+(gk[1][1]*mk[1]));
    fstareq[1][2] = -(ufk[1][2]+(grav[2]*mk[1]));
    fstareq[2][0] = -(ufk[2][0]+(gk[2][0]*mk[2]));
    fstareq[2][1] = -(ufk[2][1]+(gk[2][1]*mk[2]));
    fstareq[2][2] = -(ufk[2][2]+(grav[2]*mk[2]));
/*
Compute tstareq (torques)
*/
/*
Compute taus (RHS ignoring constraints and inertial forces)
*/
    sddovpk();
    tau[0] = (utau[0]-((((fstareq[2][0]*Vpk[0][2][0])+(fstareq[2][1]*
      Vpk[0][2][1]))-utk[2][2])+((((fstareq[0][0]*rk[0][1])-(fstareq[0][1]*
      rk[0][0]))-utk[0][2])+(((fstareq[1][0]*Vpk[0][1][0])+(fstareq[1][1]*
      Vpk[0][1][1]))-utk[1][2]))));
    tau[1] = (utau[1]-((((fstareq[1][0]*rk[1][1])-(fstareq[1][1]*rk[1][0]))-
      utk[1][2])+(((fstareq[2][0]*Vpk[1][2][0])+(fstareq[2][1]*Vpk[1][2][1]))-
      utk[2][2])));
    tau[2] = (utau[2]-(((fstareq[2][0]*rk[2][1])-(fstareq[2][1]*rk[2][0]))-
      utk[2][2]));
/*
Op counts below do not include called subroutines
*/
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   36 adds/subtracts/negates
                     21 multiplies
                      0 divides
                     12 assignments
*/
}

void sdfs0(void)
{

/*
Compute Fs (ignoring multiplier forces)
*/
    fs[0] = fs0[0];
    fs[1] = fs0[1];
    fs[2] = fs0[2];
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdfsmult(void)
{
    int i;

/*
Compute Fs (multiplier-generated forces only)
*/
    for (i = 0; i < 3; i++) {
        fs[i] = 0.;
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdfsfull(void)
{

/*
Compute Fs (including all forces)
*/
    sdfsmult();
    fs[0] = (fs[0]+fs0[0]);
    fs[1] = (fs[1]+fs0[1]);
    fs[2] = (fs[2]+fs0[2]);
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    3 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdfsgenmult(void)
{
    int i;

/*
Compute Fs (generic multiplier-generated forces)
*/
    for (i = 0; i < 3; i++) {
        fs[i] = 0.;
    }
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdfsgenfull(void)
{

/*
Compute Fs (incl generic mult & other forces)
*/
    sdfsgenmult();
    fs[0] = (fs[0]+fs0[0]);
    fs[1] = (fs[1]+fs0[1]);
    fs[2] = (fs[2]+fs0[2]);
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    3 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdfulltrq(double udotin[3],
    double multin[1],
    double trqout[3])
{
/* Compute hinge torques which would produce indicated udots
*/
    double fstarr[3][3],tstarr[3][3],Otkr[3][3],Atir[3][3],Atkr[3][3];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(61,23);
        return;
    }
/*
Account for inertial accelerations and supplied udots
*/
    Otkr[1][2] = (udotin[0]+udotin[1]);
    Otkr[2][2] = (Otkr[1][2]+udotin[2]);
    Atkr[0][0] = ((rk[0][1]*udotin[0])-(u[0]*Wkrpk[0][1]));
    Atkr[0][1] = ((u[0]*Wkrpk[0][0])-(rk[0][0]*udotin[0]));
    Atir[1][0] = (Atkr[0][0]-((ri[1][1]*udotin[0])+(u[0]*Wirk[1][1])));
    Atir[1][1] = (Atkr[0][1]+((ri[1][0]*udotin[0])+(u[0]*Wirk[1][0])));
    Atkr[1][0] = (((Atir[1][0]*c1)+(Atir[1][1]*s1))+((Otkr[1][2]*rk[1][1])-(
      wk[1][2]*Wkrpk[1][1])));
    Atkr[1][1] = (((Atir[1][1]*c1)-(Atir[1][0]*s1))+((wk[1][2]*Wkrpk[1][0])-(
      Otkr[1][2]*rk[1][0])));
    Atir[2][0] = (Atkr[1][0]-((Otkr[1][2]*ri[2][1])+(Wirk[2][1]*wk[1][2])));
    Atir[2][1] = (Atkr[1][1]+((Otkr[1][2]*ri[2][0])+(Wirk[2][0]*wk[1][2])));
    Atkr[2][0] = (((Atir[2][0]*c2)+(Atir[2][1]*s2))+((Otkr[2][2]*rk[2][1])-(
      wk[2][2]*Wkrpk[2][1])));
    Atkr[2][1] = (((Atir[2][1]*c2)-(Atir[2][0]*s2))+((wk[2][2]*Wkrpk[2][0])-(
      Otkr[2][2]*rk[2][0])));
/*
Accumulate all forces and torques
*/
    fstarr[0][0] = (ufk[0][0]+(mk[0]*(gk[0][0]-Atkr[0][0])));
    fstarr[0][1] = (ufk[0][1]+(mk[0]*(gk[0][1]-Atkr[0][1])));
    fstarr[0][2] = (ufk[0][2]+(grav[2]*mk[0]));
    fstarr[1][0] = (ufk[1][0]+(mk[1]*(gk[1][0]-Atkr[1][0])));
    fstarr[1][1] = (ufk[1][1]+(mk[1]*(gk[1][1]-Atkr[1][1])));
    fstarr[1][2] = (ufk[1][2]+(grav[2]*mk[1]));
    fstarr[2][0] = (ufk[2][0]+(mk[2]*(gk[2][0]-Atkr[2][0])));
    fstarr[2][1] = (ufk[2][1]+(mk[2]*(gk[2][1]-Atkr[2][1])));
    fstarr[2][2] = (ufk[2][2]+(grav[2]*mk[2]));
    tstarr[0][2] = (utk[0][2]-(ik[0][2][2]*udotin[0]));
    tstarr[1][2] = (utk[1][2]-(ik[1][2][2]*Otkr[1][2]));
    tstarr[2][2] = (utk[2][2]-(ik[2][2][2]*Otkr[2][2]));
/*
Now calculate the torques
*/
    sddovpk();
    trqout[0] = -(utau[0]+((tstarr[2][2]+((fstarr[2][0]*Vpk[0][2][0])+(
      fstarr[2][1]*Vpk[0][2][1])))+((tstarr[0][2]+((fstarr[0][0]*rk[0][1])-(
      fstarr[0][1]*rk[0][0])))+(tstarr[1][2]+((fstarr[1][0]*Vpk[0][1][0])+(
      fstarr[1][1]*Vpk[0][1][1]))))));
    trqout[1] = -(utau[1]+((tstarr[1][2]+((fstarr[1][0]*rk[1][1])-(fstarr[1][1]*
      rk[1][0])))+(tstarr[2][2]+((fstarr[2][0]*Vpk[1][2][0])+(fstarr[2][1]*
      Vpk[1][2][1])))));
    trqout[2] = -(utau[2]+(tstarr[2][2]+((fstarr[2][0]*rk[2][1])-(fstarr[2][1]*
      rk[2][0]))));
/*
Op counts below do not include called subroutines
*/
/*
 Used -0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   63 adds/subtracts/negates
                     52 multiplies
                      0 divides
                     27 assignments
*/
}

void sdcomptrq(double udotin[3],
    double trqout[3])
{
/* Compute hinge torques to produce these udots, ignoring constraints
*/
    double multin[1];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(60,23);
        return;
    }
    sdfulltrq(udotin,multin,trqout);
}

void sdmulttrq(double multin[1],
    double trqout[3])
{
/* Compute hinge trqs which would be produced by these mults.
*/
    int i;

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(65,23);
        return;
    }
    for (i = 0; i < 3; i++) {
        trqout[i] = 0.;
    }
}

void sdrhs(void)
{
/*
Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041
*/

/*
Compute hinge torques for tree hinges
*/
    tauc[0] = utau[0];
    tauc[1] = utau[1];
    tauc[2] = utau[2];
    sddoiner();
/*
Compute onk & onb (angular accels in N)
*/
    Onkb[1][2] = (udot[0]+udot[1]);
    Onkb[2][2] = (Onkb[1][2]+udot[2]);
    onb[0][0] = 0.;
    onb[0][1] = 0.;
    onb[0][2] = udot[0];
    onb[1][0] = 0.;
    onb[1][1] = 0.;
    onb[1][2] = Onkb[1][2];
    onb[2][0] = 0.;
    onb[2][1] = 0.;
    onb[2][2] = Onkb[2][2];
/*
Compute acceleration dyadics
*/
    dyad[0][0][0] = -w2w2[0];
    dyad[0][0][1] = -udot[0];
    dyad[0][0][2] = 0.;
    dyad[0][1][0] = udot[0];
    dyad[0][1][1] = -w2w2[0];
    dyad[0][1][2] = 0.;
    dyad[0][2][0] = 0.;
    dyad[0][2][1] = 0.;
    dyad[0][2][2] = 0.;
    dyad[1][0][0] = -w2w2[1];
    dyad[1][0][1] = -Onkb[1][2];
    dyad[1][0][2] = 0.;
    dyad[1][1][0] = Onkb[1][2];
    dyad[1][1][1] = -w2w2[1];
    dyad[1][1][2] = 0.;
    dyad[1][2][0] = 0.;
    dyad[1][2][1] = 0.;
    dyad[1][2][2] = 0.;
    dyad[2][0][0] = -w2w2[2];
    dyad[2][0][1] = -Onkb[2][2];
    dyad[2][0][2] = 0.;
    dyad[2][1][0] = Onkb[2][2];
    dyad[2][1][1] = -w2w2[2];
    dyad[2][1][2] = 0.;
    dyad[2][2][0] = 0.;
    dyad[2][2][1] = 0.;
    dyad[2][2][2] = 0.;
/*
Compute ank & anb (mass center linear accels in N)
*/
    Ankb[0][0] = (rk[0][1]*udot[0]);
    Ankb[0][1] = -(rk[0][0]*udot[0]);
    AOnkri[1][0] = (Ankb[0][0]-(ri[1][1]*udot[0]));
    AOnkri[1][1] = (Ankb[0][1]+(ri[1][0]*udot[0]));
    Ankb[1][0] = ((Onkb[1][2]*rk[1][1])+((AOnkri[1][0]*c1)+(AOnkri[1][1]*s1)));
    Ankb[1][1] = (((AOnkri[1][1]*c1)-(AOnkri[1][0]*s1))-(Onkb[1][2]*rk[1][0]));
    AOnkri[2][0] = (Ankb[1][0]-(Onkb[1][2]*ri[2][1]));
    AOnkri[2][1] = (Ankb[1][1]+(Onkb[1][2]*ri[2][0]));
    Ankb[2][0] = ((Onkb[2][2]*rk[2][1])+((AOnkri[2][0]*c2)+(AOnkri[2][1]*s2)));
    Ankb[2][1] = (((AOnkri[2][1]*c2)-(AOnkri[2][0]*s2))-(Onkb[2][2]*rk[2][0]));
    AnkAtk[0][0] = (Ankb[0][0]+Atk[0][0]);
    AnkAtk[0][1] = (Ankb[0][1]+Atk[0][1]);
    ank[0][0] = ((AnkAtk[0][0]*c0)-(AnkAtk[0][1]*s0));
    ank[0][1] = ((AnkAtk[0][0]*s0)+(AnkAtk[0][1]*c0));
    AnkAtk[1][0] = (Ankb[1][0]+Atk[1][0]);
    AnkAtk[1][1] = (Ankb[1][1]+Atk[1][1]);
    ank[1][0] = ((AnkAtk[1][0]*cnk[1][0][0])+(AnkAtk[1][1]*cnk[1][0][1]));
    ank[1][1] = ((AnkAtk[1][0]*cnk[1][1][0])+(AnkAtk[1][1]*cnk[1][1][1]));
    AnkAtk[2][0] = (Ankb[2][0]+Atk[2][0]);
    AnkAtk[2][1] = (Ankb[2][1]+Atk[2][1]);
    ank[2][0] = ((AnkAtk[2][0]*cnk[2][0][0])+(AnkAtk[2][1]*cnk[2][0][1]));
    ank[2][1] = ((AnkAtk[2][0]*cnk[2][1][0])+(AnkAtk[2][1]*cnk[2][1][1]));
    anb[0][0] = ank[0][0];
    anb[0][1] = ank[0][1];
    anb[0][2] = 0.;
    anb[1][0] = ank[1][0];
    anb[1][1] = ank[1][1];
    anb[1][2] = 0.;
    anb[2][0] = ank[2][0];
    anb[2][1] = ank[2][1];
    anb[2][2] = 0.;
    roustate = 3;
/*
 Used 0.02 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   36 adds/subtracts/negates
                     30 multiplies
                      0 divides
                     72 assignments
*/
}

void sdmassmat(double mmat[3][3])
{
/* Return the system mass matrix (LHS)
*/
    int i,j;

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(57,23);
        return;
    }
    sddomm(57);
    for (i = 0; i < 3; i++) {
        for (j = i; j <= 2; j++) {
            mmat[i][j] = mm[i][j];
            mmat[j][i] = mm[i][j];
        }
    }
}

void sdfrcmat(double fmat[3])
{
/* Return the system force matrix (RHS), excluding constraints
*/
    int i;

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(58,23);
        return;
    }
    sddofs0();
    for (i = 0; i < 3; i++) {
        fmat[i] = fs0[i];
    }
}

void sdpseudo(double lqout[1],
    double luout[1])
{
/*
Return pseudo-coordinates for loop joints.

*/
/*
There are no loop joints in this system.

*/
}

void sdpsqdot(double lqdout[1])
{
/*
Return pseudo-coordinate derivatives for loop joints.

*/
/*
There are no loop joints in this system.

*/
}

void sdpsudot(double ludout[1])
{
/*
Return pseudo-coordinate accelerations for loop joints.

*/
/*
There are no loop joints in this system.

*/
}

void sdperr(double errs[1])
{

}

void sdverr(double errs[1])
{

}

void sdaerr(double errs[1])
{

}
int 
sdchkbnum(int routine,
    int bnum)
{

    if ((bnum < -1) || (bnum > 2)) {
        sdseterr(routine,15);
        return 1;
    }
    return 0;
}
int 
sdchkjnum(int routine,
    int jnum)
{

    if ((jnum < 0) || (jnum > 2)) {
        sdseterr(routine,16);
        return 1;
    }
    return 0;
}
int 
sdchkucnum(int routine,
    int ucnum)
{

    if ((ucnum < 0) || (ucnum > -1)) {
        sdseterr(routine,21);
        return 1;
    }
    return 0;
}
int 
sdchkjaxis(int routine,
    int jnum,
    int axnum)
{
    int maxax;

    if (sdchkjnum(routine,jnum) != 0) {
        return 1;
    }
    if ((axnum < 0) || (axnum > 6)) {
        sdseterr(routine,17);
        return 1;
    }
    maxax = njntdof[jnum]-1;
    if ((jtype[jnum] == 4) || (jtype[jnum] == 6) || (jtype[jnum] == 21)) {
        maxax = maxax+1;
    }
    if (axnum > maxax) {
        sdseterr(routine,18);
        return 1;
    }
    return 0;
}
int 
sdchkjpin(int routine,
    int jnum,
    int pinno)
{
    int maxax,pinok;

    if (sdchkjnum(routine,jnum) != 0) {
        return 1;
    }
    if ((pinno < 0) || (pinno > 5)) {
        sdseterr(routine,17);
        return 1;
    }
    if (njntdof[jnum] >= 3) {
        maxax = 2;
    } else {
        maxax = njntdof[jnum]-1;
    }
    if (jtype[jnum] == 4) {
        maxax = -1;
    }
    if (jtype[jnum] == 7) {
        maxax = 0;
    }
    pinok = 0;
    if (pinno <= maxax) {
        pinok = 1;
    }
    if (pinok == 0) {
        sdseterr(routine,18);
        return 1;
    }
    return 0;
}
int 
sdindx(int joint,
    int axis)
{
    int offs,gotit;

    if (sdchkjaxis(36,joint,axis) != 0) {
        return 0;
    }
    gotit = 0;
    if (jtype[joint] == 4) {
        if (axis == 3) {
            offs = ballq[joint];
            gotit = 1;
        }
    } else {
        if ((jtype[joint] == 6) || (jtype[joint] == 21)) {
            if (axis == 6) {
                offs = ballq[joint];
                gotit = 1;
            }
        }
    }
    if (gotit == 0) {
        offs = firstq[joint]+axis;
    }
    return offs;
}

void sdpresacc(int joint,
    int axis,
    double prval)
{

}

void sdpresvel(int joint,
    int axis,
    double prval)
{

}

void sdprespos(int joint,
    int axis,
    double prval)
{

}

void sdgetht(int joint,
    int axis,
    double *torque)
{

    if (sdchkjaxis(30,joint,axis) != 0) {
        return;
    }
    if (roustate != 3) {
        sdseterr(30,24);
        return;
    }
    *torque = tauc[sdindx(joint,axis)];
}

void sdhinget(int joint,
    int axis,
    double torque)
{

    if (sdchkjaxis(10,joint,axis) != 0) {
        return;
    }
    if (roustate != 2) {
        sdseterr(10,23);
        return;
    }
    if (mfrcflg != 0) {
        mtau[sdindx(joint,axis)] = mtau[sdindx(joint,axis)]+torque;
    } else {
        fs0flg = 0;
        utau[sdindx(joint,axis)] = utau[sdindx(joint,axis)]+torque;
    }
}

void sdpointf(int body,
    double point[3],
    double force[3])
{
    double torque[3];

    if (sdchkbnum(11,body) != 0) {
        return;
    }
    if (roustate != 2) {
        sdseterr(11,23);
        return;
    }
    if (body == -1) {
        return;
    }
    torque[0] = point[1]*force[2]-point[2]*force[1];
    torque[1] = point[2]*force[0]-point[0]*force[2];
    torque[2] = point[0]*force[1]-point[1]*force[0];
    if (mfrcflg != 0) {
        mfk[body][0] = mfk[body][0]+force[0];
        mtk[body][0] = mtk[body][0]+torque[0];
        mfk[body][1] = mfk[body][1]+force[1];
        mtk[body][1] = mtk[body][1]+torque[1];
        mfk[body][2] = mfk[body][2]+force[2];
        mtk[body][2] = mtk[body][2]+torque[2];
    } else {
        fs0flg = 0;
        ufk[body][0] = ufk[body][0]+force[0];
        utk[body][0] = utk[body][0]+torque[0];
        ufk[body][1] = ufk[body][1]+force[1];
        utk[body][1] = utk[body][1]+torque[1];
        ufk[body][2] = ufk[body][2]+force[2];
        utk[body][2] = utk[body][2]+torque[2];
    }
}

void sdbodyt(int body,
    double torque[3])
{

    if (sdchkbnum(12,body) != 0) {
        return;
    }
    if (roustate != 2) {
        sdseterr(12,23);
        return;
    }
    if (body == -1) {
        return;
    }
    if (mfrcflg != 0) {
        mtk[body][0] = mtk[body][0]+torque[0];
        mtk[body][1] = mtk[body][1]+torque[1];
        mtk[body][2] = mtk[body][2]+torque[2];
    } else {
        fs0flg = 0;
        utk[body][0] = utk[body][0]+torque[0];
        utk[body][1] = utk[body][1]+torque[1];
        utk[body][2] = utk[body][2]+torque[2];
    }
}

void sddoww(int routine)
{

    roustate = 2;
    if (wwflg == 0) {
        wwflg = 1;
    }
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      0 assignments
*/
}

void sdxudot0(int routine,
    double oudot0[3])
{
/*
Compute unconstrained equations
*/
    int i;

    sdlhs(routine);
/*
Solve equations for udots
*/
    sdfs0();
    sdldubslv(3,3,mmap,works,mlo,mdi,fs,udot);
    for (i = 0; i <= 2; i++) {
        oudot0[i] = udot[i];
    }
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      3 assignments
*/
}

void sdudot0(double oudot0[3])
{

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(66,23);
        return;
    }
    sdxudot0(66,oudot0);
}

void sdsetudot(double iudot[3])
{
/*
Assign udots and advance to stage Dynamics Ready
*/
    int i;

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(68,23);
        return;
    }
    for (i = 0; i <= 2; i++) {
        udot[i] = iudot[i];
    }
    sdrhs();
}

void sdxudotm(int routine,
    double imult[1],
    double oudotm[3])
{
/*
Compute udots due only to multipliers
*/
    int i;

    sdlhs(routine);
    for (i = 0; i <= 2; i++) {
        udot[i] = 0.;
    }
    for (i = 0; i <= 2; i++) {
        oudotm[i] = udot[i];
    }
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      6 assignments
*/
}

void sdudotm(double imult[1],
    double oudotm[3])
{

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(67,23);
        return;
    }
    sdxudotm(67,imult,oudotm);
}

void sdderiv(double oqdot[3],
    double oudot[3])
{
/*
This is the derivative section for a 3-body ground-based
system with 3 hinge degree(s) of freedom.
*/
    int i;
    double udot0[3],udot1[3];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(17,23);
        return;
    }
    if (stabvelq == 1) {
        sdseterr(17,32);
    }
    if (stabposq == 1) {
        sdseterr(17,33);
    }
    wsiz = 0;
/*
Compute unconstrained equations
*/
    sdxudot0(17,udot0);
    sdrhs();
    wrank = 0;
    for (i = 0; i <= 2; i++) {
        oqdot[i] = qdot[i];
    }
    for (i = 0; i <= 2; i++) {
        oudot[i] = udot[i];
    }
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    0 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      6 assignments
*/
}
/*
Compute residuals for use with DAE integrator
*/

void sdresid(double eqdot[3],
    double eudot[3],
    double emults[1],
    double resid[6])
{
    int i;
    double uderrs[3];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(16,23);
        return;
    }
    if (stabposq == 1) {
        sdseterr(16,33);
    }
    sdfulltrq(eudot,emults,uderrs);
    for (i = 0; i < 3; i++) {
        resid[i] = eqdot[i]-qdot[i];
    }
    for (i = 0; i < 3; i++) {
        resid[3+i] = uderrs[i];
    }
    for (i = 0; i < 3; i++) {
        udot[i] = eudot[i];
    }
    sdrhs();
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain    3 adds/subtracts/negates
                      0 multiplies
                      0 divides
                      9 assignments
*/
}

void sdmult(double omults[1],
    int *owrank,
    int omultmap[1])
{

    if (roustate != 3) {
        sdseterr(34,24);
        return;
    }
    *owrank = wrank;
}

void sdreac(double force[3][3],
    double torque[3][3])
{
/*
Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041
*/

    if (roustate != 3) {
        sdseterr(31,24);
        return;
    }
/*
Compute reaction forces for non-weld tree joints
*/
    fc[2][0] = ((mk[2]*(AnkAtk[2][0]-gk[2][0]))-ufk[2][0]);
    fc[2][1] = ((mk[2]*(AnkAtk[2][1]-gk[2][1]))-ufk[2][1]);
    fc[2][2] = -(ufk[2][2]+(grav[2]*mk[2]));
    tc[2][0] = -(utk[2][0]+((fc[2][2]*rk[2][1])-(fc[2][1]*rk[2][2])));
    tc[2][1] = -(utk[2][1]+((fc[2][0]*rk[2][2])-(fc[2][2]*rk[2][0])));
    tc[2][2] = ((ik[2][2][2]*Onkb[2][2])-(utk[2][2]+((fc[2][1]*rk[2][0])-(
      fc[2][0]*rk[2][1]))));
    fccikt[2][0] = ((fc[2][0]*c2)-(fc[2][1]*s2));
    fccikt[2][1] = ((fc[2][0]*s2)+(fc[2][1]*c2));
    fccikt[2][2] = fc[2][2];
    ffk[1][0] = (ufk[1][0]-fccikt[2][0]);
    ffk[1][1] = (ufk[1][1]-fccikt[2][1]);
    ffk[1][2] = (ufk[1][2]-fccikt[2][2]);
    ttk[1][0] = (utk[1][0]-(((fccikt[2][2]*ri[2][1])-(fccikt[2][1]*ri[2][2]))+((
      tc[2][0]*c2)-(tc[2][1]*s2))));
    ttk[1][1] = (utk[1][1]-(((fccikt[2][0]*ri[2][2])-(fccikt[2][2]*ri[2][0]))+((
      tc[2][0]*s2)+(tc[2][1]*c2))));
    ttk[1][2] = (utk[1][2]-(tc[2][2]+((fccikt[2][1]*ri[2][0])-(fccikt[2][0]*
      ri[2][1]))));
    fc[1][0] = ((mk[1]*(AnkAtk[1][0]-gk[1][0]))-ffk[1][0]);
    fc[1][1] = ((mk[1]*(AnkAtk[1][1]-gk[1][1]))-ffk[1][1]);
    fc[1][2] = -(ffk[1][2]+(grav[2]*mk[1]));
    tc[1][0] = -(ttk[1][0]+((fc[1][2]*rk[1][1])-(fc[1][1]*rk[1][2])));
    tc[1][1] = -(ttk[1][1]+((fc[1][0]*rk[1][2])-(fc[1][2]*rk[1][0])));
    tc[1][2] = ((ik[1][2][2]*Onkb[1][2])-(ttk[1][2]+((fc[1][1]*rk[1][0])-(
      fc[1][0]*rk[1][1]))));
    fccikt[1][0] = ((fc[1][0]*c1)-(fc[1][1]*s1));
    fccikt[1][1] = ((fc[1][0]*s1)+(fc[1][1]*c1));
    fccikt[1][2] = fc[1][2];
    ffk[0][0] = (ufk[0][0]-fccikt[1][0]);
    ffk[0][1] = (ufk[0][1]-fccikt[1][1]);
    ffk[0][2] = (ufk[0][2]-fccikt[1][2]);
    ttk[0][0] = (utk[0][0]-(((fccikt[1][2]*ri[1][1])-(fccikt[1][1]*ri[1][2]))+((
      tc[1][0]*c1)-(tc[1][1]*s1))));
    ttk[0][1] = (utk[0][1]-(((fccikt[1][0]*ri[1][2])-(fccikt[1][2]*ri[1][0]))+((
      tc[1][0]*s1)+(tc[1][1]*c1))));
    ttk[0][2] = (utk[0][2]-(tc[1][2]+((fccikt[1][1]*ri[1][0])-(fccikt[1][0]*
      ri[1][1]))));
    fc[0][0] = ((mk[0]*(AnkAtk[0][0]-gk[0][0]))-ffk[0][0]);
    fc[0][1] = ((mk[0]*(AnkAtk[0][1]-gk[0][1]))-ffk[0][1]);
    fc[0][2] = -(ffk[0][2]+(grav[2]*mk[0]));
    tc[0][0] = -(ttk[0][0]+((fc[0][2]*rk[0][1])-(fc[0][1]*rk[0][2])));
    tc[0][1] = -(ttk[0][1]+((fc[0][0]*rk[0][2])-(fc[0][2]*rk[0][0])));
    tc[0][2] = ((ik[0][2][2]*udot[0])-(ttk[0][2]+((fc[0][1]*rk[0][0])-(fc[0][0]*
      rk[0][1]))));
    force[0][0] = fc[0][0];
    torque[0][0] = tc[0][0];
    force[0][1] = fc[0][1];
    torque[0][1] = tc[0][1];
    force[0][2] = fc[0][2];
    torque[0][2] = tc[0][2];
    force[1][0] = fc[1][0];
    torque[1][0] = tc[1][0];
    force[1][1] = fc[1][1];
    torque[1][1] = tc[1][1];
    force[1][2] = fc[1][2];
    torque[1][2] = tc[1][2];
    force[2][0] = fc[2][0];
    torque[2][0] = tc[2][0];
    force[2][1] = fc[2][1];
    torque[2][1] = tc[2][1];
    force[2][2] = fc[2][2];
    torque[2][2] = tc[2][2];
/*
Compute reaction forces for tree weld joints
*/
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   77 adds/subtracts/negates
                     58 multiplies
                      0 divides
                     54 assignments
*/
}

void sdmom(double lm[3],
    double am[3],
    double *ke)
{
/*
Compute system linear and angular momentum, and kinetic energy.

Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041
*/
    double lk[3][3],hnk[3][3];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(19,23);
        return;
    }
    lk[0][0] = (mk[0]*vnk[0][0]);
    lk[0][1] = (mk[0]*vnk[0][1]);
    lk[1][0] = (mk[1]*vnk[1][0]);
    lk[1][1] = (mk[1]*vnk[1][1]);
    lk[2][0] = (mk[2]*vnk[2][0]);
    lk[2][1] = (mk[2]*vnk[2][1]);
    hnk[0][2] = (ik[0][2][2]*u[0]);
    hnk[1][2] = (ik[1][2][2]*wk[1][2]);
    hnk[2][2] = (ik[2][2][2]*wk[2][2]);
    lm[0] = (lk[2][0]+(lk[0][0]+lk[1][0]));
    lm[1] = (lk[2][1]+(lk[0][1]+lk[1][1]));
    lm[2] = 0.;
    am[0] = ((com[2]*lm[1])+(((lk[0][1]*rk[0][2])-(lk[1][1]*rnk[1][2]))-(
      lk[2][1]*rnk[2][2])));
    am[1] = (((lk[2][0]*rnk[2][2])+((lk[1][0]*rnk[1][2])-(lk[0][0]*rk[0][2])))-(
      com[2]*lm[0]));
    am[2] = (((hnk[2][2]+((lk[2][1]*rnk[2][0])-(lk[2][0]*rnk[2][1])))+((
      hnk[0][2]+((lk[0][1]*rnk[0][0])-(lk[0][0]*rnk[0][1])))+(hnk[1][2]+((
      lk[1][1]*rnk[1][0])-(lk[1][0]*rnk[1][1])))))-((com[0]*lm[1])-(com[1]*lm[0]
      )));
    *ke = (.5*(((hnk[2][2]*wk[2][2])+((lk[2][0]*vnk[2][0])+(lk[2][1]*vnk[2][1]))
      )+(((hnk[0][2]*u[0])+((lk[0][0]*vnk[0][0])+(lk[0][1]*vnk[0][1])))+((
      hnk[1][2]*wk[1][2])+((lk[1][0]*vnk[1][0])+(lk[1][1]*vnk[1][1]))))));
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   28 adds/subtracts/negates
                     35 multiplies
                      0 divides
                     16 assignments
*/
}

void sdsys(double *mtoto,
    double cm[3],
    double icm[3][3])
{
/*
Compute system total mass, and instantaneous center of mass and
inertia matrix.

Generated 09-Aug-2008 11:36:26 by SD/FAST, Kane's formulation
(sdfast B.2.8 #30123) on machine ID unknown
Copyright (c) 1990-1997 Symbolic Dynamics, Inc.
Copyright (c) 1990-1997 Parametric Technology Corp.
RESTRICTED RIGHTS LEGEND: Use, duplication, or disclosure by the U.S.
Government is subject to restrictions as set forth in subparagraph
(c)(1)(ii) of the Rights in Technical Data and Computer Software
clause at DFARS 52.227-7013 and similar clauses in the FAR and NASA
FAR Supplement.  Symbolic Dynamics, Inc., Mountain View, CA 94041
*/
    double ikcnkt[3][3][3];

    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(20,23);
        return;
    }
    *mtoto = mtot;
    cm[0] = com[0];
    cm[1] = com[1];
    cm[2] = com[2];
    ikcnkt[0][0][0] = (ik[0][0][0]*c0);
    ikcnkt[0][0][1] = (ik[0][0][0]*s0);
    ikcnkt[0][1][0] = -(ik[0][1][1]*s0);
    ikcnkt[0][1][1] = (ik[0][1][1]*c0);
    ikcnkt[1][0][0] = (cnk[1][0][0]*ik[1][0][0]);
    ikcnkt[1][0][1] = (cnk[1][1][0]*ik[1][0][0]);
    ikcnkt[1][1][0] = (cnk[1][0][1]*ik[1][1][1]);
    ikcnkt[1][1][1] = (cnk[1][1][1]*ik[1][1][1]);
    ikcnkt[2][0][0] = (cnk[2][0][0]*ik[2][0][0]);
    ikcnkt[2][0][1] = (cnk[2][1][0]*ik[2][0][0]);
    ikcnkt[2][1][0] = (cnk[2][0][1]*ik[2][1][1]);
    ikcnkt[2][1][1] = (cnk[2][1][1]*ik[2][1][1]);
    icm[0][0] = ((((mk[2]*((rnk[2][1]*rnk[2][1])+(rnk[2][2]*rnk[2][2])))+((
      cnk[2][0][0]*ikcnkt[2][0][0])+(cnk[2][0][1]*ikcnkt[2][1][0])))+(((mk[0]*((
      rk[0][2]*rk[0][2])+(rnk[0][1]*rnk[0][1])))+((ikcnkt[0][0][0]*c0)-(
      ikcnkt[0][1][0]*s0)))+((mk[1]*((rnk[1][1]*rnk[1][1])+(rnk[1][2]*rnk[1][2])
      ))+((cnk[1][0][0]*ikcnkt[1][0][0])+(cnk[1][0][1]*ikcnkt[1][1][0])))))-(
      mtot*((com[1]*com[1])+(com[2]*com[2]))));
    icm[0][1] = ((mtot*(com[0]*com[1]))+((((cnk[2][0][0]*ikcnkt[2][0][1])+(
      cnk[2][0][1]*ikcnkt[2][1][1]))-(mk[2]*(rnk[2][0]*rnk[2][1])))+((((
      cnk[1][0][0]*ikcnkt[1][0][1])+(cnk[1][0][1]*ikcnkt[1][1][1]))-(mk[1]*(
      rnk[1][0]*rnk[1][1])))+(((ikcnkt[0][0][1]*c0)-(ikcnkt[0][1][1]*s0))-(mk[0]
      *(rnk[0][0]*rnk[0][1]))))));
    icm[0][2] = ((mtot*(com[0]*com[2]))+(((mk[0]*(rk[0][2]*rnk[0][0]))-(mk[1]*(
      rnk[1][0]*rnk[1][2])))-(mk[2]*(rnk[2][0]*rnk[2][2]))));
    icm[1][0] = icm[0][1];
    icm[1][1] = ((((mk[2]*((rnk[2][0]*rnk[2][0])+(rnk[2][2]*rnk[2][2])))+((
      cnk[2][1][0]*ikcnkt[2][0][1])+(cnk[2][1][1]*ikcnkt[2][1][1])))+(((mk[0]*((
      rk[0][2]*rk[0][2])+(rnk[0][0]*rnk[0][0])))+((ikcnkt[0][0][1]*s0)+(
      ikcnkt[0][1][1]*c0)))+((mk[1]*((rnk[1][0]*rnk[1][0])+(rnk[1][2]*rnk[1][2])
      ))+((cnk[1][1][0]*ikcnkt[1][0][1])+(cnk[1][1][1]*ikcnkt[1][1][1])))))-(
      mtot*((com[0]*com[0])+(com[2]*com[2]))));
    icm[1][2] = ((mtot*(com[1]*com[2]))+(((mk[0]*(rk[0][2]*rnk[0][1]))-(mk[1]*(
      rnk[1][1]*rnk[1][2])))-(mk[2]*(rnk[2][1]*rnk[2][2]))));
    icm[2][0] = icm[0][2];
    icm[2][1] = icm[1][2];
    icm[2][2] = (((ik[2][2][2]+(mk[2]*((rnk[2][0]*rnk[2][0])+(rnk[2][1]*
      rnk[2][1]))))+((ik[0][2][2]+(mk[0]*((rnk[0][0]*rnk[0][0])+(rnk[0][1]*
      rnk[0][1]))))+(ik[1][2][2]+(mk[1]*((rnk[1][0]*rnk[1][0])+(rnk[1][1]*
      rnk[1][1]))))))-(mtot*((com[0]*com[0])+(com[1]*com[1]))));
/*
 Used 0.00 seconds CPU time,
 0 additional bytes of memory.
 Equations contain   52 adds/subtracts/negates
                     90 multiplies
                      0 divides
                     25 assignments
*/
}

void sdpos(int body,
    double pt[3],
    double loc[3])
{
/*
Return inertial frame location of a point on a body.

*/
    double pv[3];

    if (sdchkbnum(21,body) != 0) {
        return;
    }
    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(21,23);
        return;
    }
    if (body  ==  -1) {
        loc[0] = pt[0];
        loc[1] = pt[1];
        loc[2] = pt[2];
    } else {
        pv[0] = rnb[body][0]+pt[0]*cnb[body][0][0]+pt[1]*cnb[body][0][1]+pt[2]*
          cnb[body][0][2];
        pv[1] = rnb[body][1]+pt[0]*cnb[body][1][0]+pt[1]*cnb[body][1][1]+pt[2]*
          cnb[body][1][2];
        pv[2] = rnb[body][2]+pt[0]*cnb[body][2][0]+pt[1]*cnb[body][2][1]+pt[2]*
          cnb[body][2][2];
        loc[0] = pv[0];
        loc[1] = pv[1];
        loc[2] = pv[2];
    }
}

void sdvel(int body,
    double pt[3],
    double velo[3])
{
/*
Return inertial frame velocity of a point on a body.

*/
    double pv[3];

    if (sdchkbnum(22,body) != 0) {
        return;
    }
    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(22,23);
        return;
    }
    if (body  ==  -1) {
        velo[0] = 0.;
        velo[1] = 0.;
        velo[2] = 0.;
    } else {
        pv[0] = wb[body][1]*pt[2]-wb[body][2]*pt[1];
        pv[1] = wb[body][2]*pt[0]-wb[body][0]*pt[2];
        pv[2] = wb[body][0]*pt[1]-wb[body][1]*pt[0];
        velo[0] = vnb[body][0]+pv[0]*cnb[body][0][0]+pv[1]*cnb[body][0][1]+pv[2]
          *cnb[body][0][2];
        velo[1] = vnb[body][1]+pv[0]*cnb[body][1][0]+pv[1]*cnb[body][1][1]+pv[2]
          *cnb[body][1][2];
        velo[2] = vnb[body][2]+pv[0]*cnb[body][2][0]+pv[1]*cnb[body][2][1]+pv[2]
          *cnb[body][2][2];
    }
}

void sdorient(int body,
    double dircos[3][3])
{
/*
Return orientation of body w.r.t. ground frame.

*/

    if (sdchkbnum(23,body) != 0) {
        return;
    }
    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(23,23);
        return;
    }
    if (body == -1) {
        dircos[0][0] = 1.;
        dircos[0][1] = 0.;
        dircos[0][2] = 0.;
        dircos[1][0] = 0.;
        dircos[1][1] = 1.;
        dircos[1][2] = 0.;
        dircos[2][0] = 0.;
        dircos[2][1] = 0.;
        dircos[2][2] = 1.;
    } else {
        dircos[0][0] = cnb[body][0][0];
        dircos[0][1] = cnb[body][0][1];
        dircos[0][2] = cnb[body][0][2];
        dircos[1][0] = cnb[body][1][0];
        dircos[1][1] = cnb[body][1][1];
        dircos[1][2] = cnb[body][1][2];
        dircos[2][0] = cnb[body][2][0];
        dircos[2][1] = cnb[body][2][1];
        dircos[2][2] = cnb[body][2][2];
    }
}

void sdangvel(int body,
    double avel[3])
{
/*
Return angular velocity of the body.

*/

    if (sdchkbnum(24,body) != 0) {
        return;
    }
    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(24,23);
        return;
    }
    if (body == -1) {
        avel[0] = 0.;
        avel[1] = 0.;
        avel[2] = 0.;
    } else {
        avel[0] = wb[body][0];
        avel[1] = wb[body][1];
        avel[2] = wb[body][2];
    }
}

void sdtrans(int frbod,
    double ivec[3],
    int tobod,
    double ovec[3])
{
/*
Transform ivec from frbod frame to tobod frame.

*/
    double pv[3];

    if (sdchkbnum(25,frbod) != 0) {
        return;
    }
    if (sdchkbnum(25,tobod) != 0) {
        return;
    }
    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(25,23);
        return;
    }
    if (frbod == tobod) {
        sdvcopy(ivec,ovec);
        return;
    }
    if (frbod == -1) {
        sdvcopy(ivec,pv);
        ovec[0] = pv[0]*cnb[tobod][0][0]+pv[1]*cnb[tobod][1][0]+pv[2]*cnb[tobod
          ][2][0];
        ovec[1] = pv[0]*cnb[tobod][0][1]+pv[1]*cnb[tobod][1][1]+pv[2]*cnb[tobod
          ][2][1];
        ovec[2] = pv[0]*cnb[tobod][0][2]+pv[1]*cnb[tobod][1][2]+pv[2]*cnb[tobod
          ][2][2];
        return;
    }
    if (tobod == -1) {
        sdvcopy(ivec,pv);
        ovec[0] = pv[0]*cnb[frbod][0][0]+pv[1]*cnb[frbod][0][1]+pv[2]*cnb[frbod
          ][0][2];
        ovec[1] = pv[0]*cnb[frbod][1][0]+pv[1]*cnb[frbod][1][1]+pv[2]*cnb[frbod
          ][1][2];
        ovec[2] = pv[0]*cnb[frbod][2][0]+pv[1]*cnb[frbod][2][1]+pv[2]*cnb[frbod
          ][2][2];
        return;
    }
    pv[0] = ivec[0]*cnb[frbod][0][0]+ivec[1]*cnb[frbod][0][1]+ivec[2]*cnb[frbod
      ][0][2];
    pv[1] = ivec[0]*cnb[frbod][1][0]+ivec[1]*cnb[frbod][1][1]+ivec[2]*cnb[frbod
      ][1][2];
    pv[2] = ivec[0]*cnb[frbod][2][0]+ivec[1]*cnb[frbod][2][1]+ivec[2]*cnb[frbod
      ][2][2];
    ovec[0] = pv[0]*cnb[tobod][0][0]+pv[1]*cnb[tobod][1][0]+pv[2]*cnb[tobod][2][
      0];
    ovec[1] = pv[0]*cnb[tobod][0][1]+pv[1]*cnb[tobod][1][1]+pv[2]*cnb[tobod][2][
      1];
    ovec[2] = pv[0]*cnb[tobod][0][2]+pv[1]*cnb[tobod][1][2]+pv[2]*cnb[tobod][2][
      2];
}

void sdrel2cart(int coord,
    int body,
    double point[3],
    double linchg[3],
    double rotchg[3])
{
/* Return derivative of pt loc and body orient w.r.t. hinge rate
*/
    int x,i,gnd;
    double lin[3],pv[3];

    if ((coord < 0) || (coord > 2)) {
        sdseterr(59,45);
        return;
    }
    if (sdchkbnum(59,body) != 0) {
        return;
    }
    if ((roustate != 2) && (roustate != 3)) {
        sdseterr(59,23);
        return;
    }
    gnd = -1;
    if (body == gnd) {
        x = -1;
    } else {
        x = firstq[body]+njntdof[body]-1;
    }
    if (x < coord) {
        sdvset(0.,0.,0.,linchg);
        sdvset(0.,0.,0.,rotchg);
        return;
    }
    sddovpk();
    for (i = 0; i < 3; i++) {
        rotchg[i] = Wpk[coord][x][i];
        lin[i] = Vpk[coord][x][i];
    }
    if (body == gnd) {
        sdvcopy(point,pv);
    } else {
        pv[0] = rcom[body][0]+point[0];
        pv[1] = rcom[body][1]+point[1];
        pv[2] = rcom[body][2]+point[2];
    }
    sdvcross(rotchg,pv,linchg);
    sdvadd(linchg,lin,linchg);
}

void sdacc(int body,
    double pt[3],
    double accel[3])
{
/*
Return linear acceleration a point of the specified body.

*/
    double pv[3];

    if (sdchkbnum(32,body) != 0) {
        return;
    }
    if (roustate != 3) {
        sdseterr(32,24);
        return;
    }
    if (body  ==  -1) {
        accel[0] = 0.;
        accel[1] = 0.;
        accel[2] = 0.;
    } else {
        pv[0] = pt[0]*dyad[body][0][0]+pt[1]*dyad[body][0][1]+pt[2]*dyad[body][0
          ][2];
        pv[1] = pt[0]*dyad[body][1][0]+pt[1]*dyad[body][1][1]+pt[2]*dyad[body][1
          ][2];
        pv[2] = pt[0]*dyad[body][2][0]+pt[1]*dyad[body][2][1]+pt[2]*dyad[body][2
          ][2];
        accel[0] = anb[body][0]+pv[0]*cnb[body][0][0]+pv[1]*cnb[body][0][1]+pv[2
          ]*cnb[body][0][2];
        accel[1] = anb[body][1]+pv[0]*cnb[body][1][0]+pv[1]*cnb[body][1][1]+pv[2
          ]*cnb[body][1][2];
        accel[2] = anb[body][2]+pv[0]*cnb[body][2][0]+pv[1]*cnb[body][2][1]+pv[2
          ]*cnb[body][2][2];
    }
}

void sdangacc(int body,
    double aacc[3])
{
/*
Return angular acceleration of the body.

*/

    if (sdchkbnum(33,body) != 0) {
        return;
    }
    if (roustate != 3) {
        sdseterr(33,24);
        return;
    }
    if (body == -1) {
        aacc[0] = 0.;
        aacc[1] = 0.;
        aacc[2] = 0.;
    } else {
        aacc[0] = onb[body][0];
        aacc[1] = onb[body][1];
        aacc[2] = onb[body][2];
    }
}

void sdgrav(double gravin[3])
{

    grav[0] = gravin[0];
    gravq[0] = 3;
    grav[1] = gravin[1];
    gravq[1] = 3;
    grav[2] = gravin[2];
    gravq[2] = 3;
    roustate = 0;
}

void sdmass(int body,
    double massin)
{

    if (sdchkbnum(2,body) != 0) {
        return;
    }
    if (body == -1) {
        sdseterr(2,15);
        return;
    }
    if (mkq[body] != 0) {
        mk[body] = massin;
        mkq[body] = 3;
    } else {
        sdseterr(2,19);
    }
    roustate = 0;
}

void sdiner(int body,
    double inerin[3][3])
{
    int anyques;

    if (sdchkbnum(3,body) != 0) {
        return;
    }
    if (body == -1) {
        sdseterr(3,15);
        return;
    }
    anyques = 0;
    if (ikq[body][0][0]  !=  0) {
        ik[body][0][0] = inerin[0][0];
        ikq[body][0][0] = 3;
        anyques = 1;
    }
    if (ikq[body][0][1]  !=  0) {
        ik[body][0][1] = inerin[0][1];
        ikq[body][0][1] = 3;
        ik[body][1][0] = inerin[0][1];
        ikq[body][1][0] = 3;
        anyques = 1;
    }
    if (ikq[body][0][2]  !=  0) {
        ik[body][0][2] = inerin[0][2];
        ikq[body][0][2] = 3;
        ik[body][2][0] = inerin[0][2];
        ikq[body][2][0] = 3;
        anyques = 1;
    }
    if (ikq[body][1][1]  !=  0) {
        ik[body][1][1] = inerin[1][1];
        ikq[body][1][1] = 3;
        anyques = 1;
    }
    if (ikq[body][1][2]  !=  0) {
        ik[body][1][2] = inerin[1][2];
        ikq[body][1][2] = 3;
        ik[body][2][1] = inerin[1][2];
        ikq[body][2][1] = 3;
        anyques = 1;
    }
    if (ikq[body][2][2]  !=  0) {
        ik[body][2][2] = inerin[2][2];
        ikq[body][2][2] = 3;
        anyques = 1;
    }
    if (anyques == 0) {
        sdseterr(3,19);
    }
    roustate = 0;
}

void sdbtj(int joint,
    double btjin[3])
{
    int anyques;

    if (sdchkjnum(4,joint) != 0) {
        return;
    }
    anyques = 0;
    if (rkq[joint][0]  !=  0) {
        rk[joint][0] = btjin[0];
        rkq[joint][0] = 3;
        anyques = 1;
    }
    if (rkq[joint][1]  !=  0) {
        rk[joint][1] = btjin[1];
        rkq[joint][1] = 3;
        anyques = 1;
    }
    if (rkq[joint][2]  !=  0) {
        rk[joint][2] = btjin[2];
        rkq[joint][2] = 3;
        anyques = 1;
    }
    if (anyques == 0) {
        sdseterr(4,19);
    }
    roustate = 0;
}

void sditj(int joint,
    double itjin[3])
{
    int anyques;

    if (sdchkjnum(5,joint) != 0) {
        return;
    }
    anyques = 0;
    if (riq[joint][0]  !=  0) {
        ri[joint][0] = itjin[0];
        riq[joint][0] = 3;
        anyques = 1;
    }
    if (riq[joint][1]  !=  0) {
        ri[joint][1] = itjin[1];
        riq[joint][1] = 3;
        anyques = 1;
    }
    if (riq[joint][2]  !=  0) {
        ri[joint][2] = itjin[2];
        riq[joint][2] = 3;
        anyques = 1;
    }
    if (anyques == 0) {
        sdseterr(5,19);
    }
    roustate = 0;
}

void sdpin(int joint,
    int pinno,
    double pinin[3])
{
    int anyques,offs;

    if (sdchkjpin(6,joint,pinno) != 0) {
        return;
    }
    anyques = 0;
    offs = firstq[joint]+pinno;
    if (jtype[joint] == 21) {
        offs = offs+3;
    }
    if (jtype[joint] == 11) {
        offs = offs+1;
    }
    if (pinq[offs][0]  !=  0) {
        pin[offs][0] = pinin[0];
        pinq[offs][0] = 3;
        anyques = 1;
    }
    if (pinq[offs][1]  !=  0) {
        pin[offs][1] = pinin[1];
        pinq[offs][1] = 3;
        anyques = 1;
    }
    if (pinq[offs][2]  !=  0) {
        pin[offs][2] = pinin[2];
        pinq[offs][2] = 3;
        anyques = 1;
    }
    if (anyques == 0) {
        sdseterr(6,19);
    }
    roustate = 0;
}

void sdpres(int joint,
    int axis,
    int presin)
{
    int anyques;

    if (sdchkjaxis(37,joint,axis) != 0) {
        return;
    }
    if ((presin != 0) && (presin != 1)) {
        sdseterr(37,20);
    }
    anyques = 0;
    if (presq[sdindx(joint,axis)]  !=  0) {
        if (presin  !=  0) {
            pres[sdindx(joint,axis)] = 1.;
        } else {
            pres[sdindx(joint,axis)] = 0.;
        }
        presq[sdindx(joint,axis)] = 3;
        anyques = 1;
    }
    if (anyques == 0) {
        sdseterr(37,19);
    }
    wwflg = 0;
}

void sdconschg(void)
{

    wwflg = 0;
}

void sdstab(double velin,
    double posin)
{

    stabvel = velin;
    stabvelq = 3;
    stabpos = posin;
    stabposq = 3;
}

void sdgetgrav(double gravout[3])
{

    gravout[0] = grav[0];
    gravout[1] = grav[1];
    gravout[2] = grav[2];
}

void sdgetmass(int body,
    double *massout)
{

    if (sdchkbnum(40,body) != 0) {
        return;
    }
    if (body == -1) {
        sdseterr(40,15);
        return;
    }
    *massout = mk[body];
}

void sdgetiner(int body,
    double inerout[3][3])
{

    if (sdchkbnum(41,body) != 0) {
        return;
    }
    if (body == -1) {
        sdseterr(41,15);
        return;
    }
    inerout[0][0] = ik[body][0][0];
    inerout[0][1] = ik[body][0][1];
    inerout[0][2] = ik[body][0][2];
    inerout[1][0] = ik[body][1][0];
    inerout[1][1] = ik[body][1][1];
    inerout[1][2] = ik[body][1][2];
    inerout[2][0] = ik[body][2][0];
    inerout[2][1] = ik[body][2][1];
    inerout[2][2] = ik[body][2][2];
}

void sdgetbtj(int joint,
    double btjout[3])
{

    if (sdchkjnum(42,joint) != 0) {
        return;
    }
    btjout[0] = rk[joint][0];
    btjout[1] = rk[joint][1];
    btjout[2] = rk[joint][2];
}

void sdgetitj(int joint,
    double itjout[3])
{

    if (sdchkjnum(43,joint) != 0) {
        return;
    }
    itjout[0] = ri[joint][0];
    itjout[1] = ri[joint][1];
    itjout[2] = ri[joint][2];
}

void sdgetpin(int joint,
    int pinno,
    double pinout[3])
{
    int offs;

    if (sdchkjpin(44,joint,pinno) != 0) {
        return;
    }
    offs = firstq[joint]+pinno;
    if (jtype[joint] == 21) {
        offs = offs+3;
    }
    if (jtype[joint] == 11) {
        offs = offs+1;
    }
    pinout[0] = pin[offs][0];
    pinout[1] = pin[offs][1];
    pinout[2] = pin[offs][2];
}

void sdgetpres(int joint,
    int axis,
    int *presout)
{

    if (sdchkjaxis(45,joint,axis) != 0) {
        return;
    }
    if (pres[sdindx(joint,axis)]  !=  0.) {
        *presout = 1;
    } else {
        *presout = 0;
    }
}

void sdgetstab(double *velout,
    double *posout)
{

    *velout = stabvel;
    *posout = stabpos;
}

void sdinfo(int info[50])
{

    info[0] = ground;
    info[1] = nbod;
    info[2] = ndof;
    info[3] = ncons;
    info[4] = nloop;
    info[5] = nldof;
    info[6] = nloopc;
    info[7] = nball;
    info[8] = nlball;
    info[9] = npres;
    info[10] = nuser;
    info[11] = 0;
/* info entries from 12-49 are reserved */
}

void sdjnt(int joint,
    int info[50],
    int tran[6])
{
    int i,offs;

    if (sdchkjnum(48,joint) != 0) {
        return;
    }
    info[0] = jtype[joint];
    info[1] = 0;
    offs = 0;
    info[2] = inb[joint];
    info[3] = outb[joint];
    info[4] = njntdof[joint];
    info[5] = njntc[joint];
    info[6] = njntp[joint];
    info[7] = firstq[joint];
    info[8] = ballq[joint];
    info[9] = firstm[joint];
    info[10] = firstp[joint];
/* info entries from 11-49 are reserved */

    for (i = 0; i <= 5; i++) {
        if (i  <  njntdof[joint]) {
            tran[i] = trans[offs+firstq[joint]+i];
        } else {
            tran[i] = -1;
        }
    }
}

void sdcons(int consno,
    int info[50])
{

    if (sdchkucnum(49,consno) != 0) {
        return;
    }
/* There are no user constraints in this problem. */
}

void sdgentime(int *gentm)
{

    *gentm = 113626;
}
/*
Done. CPU seconds used: 0.05  Memory used: 1753088 bytes.
Equation complexity:
  sdstate:    81 adds   104 multiplies     3 divides   139 assignments
  sdderiv:   143 adds   130 multiplies     3 divides   201 assignments
  sdresid:   123 adds    94 multiplies     0 divides   132 assignments
  sdreac:     77 adds    58 multiplies     0 divides    54 assignments
  sdmom:      28 adds    35 multiplies     0 divides    16 assignments
  sdsys:      52 adds    90 multiplies     0 divides    25 assignments
*/

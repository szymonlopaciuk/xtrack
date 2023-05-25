# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

from cpymad.madx import Madx

# hllhc15 can be found at git@github.com:lhcopt/hllhc15.git

mad = Madx()

mad.input("""
call,file="../../../hllhc15/util/lhc.seq";
call,file="../../../hllhc15/hllhc_sequence.madx";
call,file="../../../hllhc15/toolkit/macro.madx";
seqedit,sequence=lhcb1;flatten;cycle,start=IP3;flatten;endedit;
seqedit,sequence=lhcb2;flatten;cycle,start=IP3;flatten;endedit;
exec,mk_beam(7000);
call,file="../../../hllhc15/round/opt_round_150_1500.madx";
on_x1 = 1;
on_x5 = 1;
on_disp = 0;
exec,check_ip(b1);
exec,check_ip(b2);
""")

mad.input('''

brho:=NRJ*1e9/clight;
! Octupoles, added by hand
i_oct_b1=0;
  kof.a12b1:=kmax_mo*i_oct_b1/imax_mo/brho; kof.a23b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kof.a34b1:=kmax_mo*i_oct_b1/imax_mo/brho; kof.a45b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kof.a56b1:=kmax_mo*i_oct_b1/imax_mo/brho; kof.a67b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kof.a78b1:=kmax_mo*i_oct_b1/imax_mo/brho; kof.a81b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kod.a12b1:=kmax_mo*i_oct_b1/imax_mo/brho; kod.a23b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kod.a34b1:=kmax_mo*i_oct_b1/imax_mo/brho; kod.a45b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kod.a56b1:=kmax_mo*i_oct_b1/imax_mo/brho; kod.a67b1:=kmax_mo*i_oct_b1/imax_mo/brho;
  kod.a78b1:=kmax_mo*i_oct_b1/imax_mo/brho; kod.a81b1:=kmax_mo*i_oct_b1/imax_mo/brho;
''')


mad.use(sequence="lhcb1")
mad.globals['vrf400'] = 16
mad.globals['lagrf400.b1'] = 0.5
mad.globals['lagrf400.b2'] = 0.
mad.twiss()
mad.save(sequence=['lhcb1', 'lhcb2'], beam=True, file="sequence_thick.madx")

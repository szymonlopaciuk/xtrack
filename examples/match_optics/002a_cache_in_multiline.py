import time

import numpy as np

import xtrack as xt
import xdeps as xd

# xt._print.suppress = True

# Load the line
collider = xt.Multiline.from_json(
    '../../test_data/hllhc15_collider/collider_00_from_mad.json')
collider.build_trackers()

collider.lhcb1.twiss_default['method'] = '4d'
collider.lhcb2.twiss_default['method'] = '4d'
collider.lhcb2.twiss_default['reverse'] = True

collider.vars['on_x2'] = 123
collider.vars.cache = True

collider.vars['on_x1'] = 2
collider.vars['on_x5'] = 3

collider._var_sharing = None



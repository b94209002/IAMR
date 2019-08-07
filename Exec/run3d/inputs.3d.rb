
#*******************************************************************************
# INPUTS.3D.EULER
#*******************************************************************************

max_step 		= 2000

stop_time 		= 10.0

amr.n_cell              = 128 128 32
amr.max_grid_size       = 32

amr.max_level           = 0

ns.do_liquid_ref        = 0
ns.do_vorticity_ref     = 0

amr.regrid_int		= -1
#amr.regrid_file         = fixed_grids_rb_z32

ns.v                    = 1
amr.v                   = 1

amr.check_int		= 100 
#amr.restart		= chk03000

amr.plot_int		= 5
#amr.plot_per            = 10.

ns.cfl                  = 0.5  # CFL number used to set dt

ns.init_shrink          = 1.0  # factor which multiplies the very first time step
ns.init_iter            = 0
ns.do_init_proj         = 0

ns.be_cn_theta          = 0.5

ns.vel_visc_coef        = 0.001
ns.scal_diff_coefs      = 0.0014 0.0014

ns.do_cons_trac         = 1
ns.do_trac2             = 1
ns.do_cons_trac2        = 1


amr.probin_file 	= probin.3d.rb

geometry.coord_sys   =  0

geometry.prob_lo     =  -2. -2. 0.
geometry.prob_hi     =  2. 2. 1.

geometry.is_periodic =  1 1 0

ns.gravity = 0.0

ns.lo_bc             = 0 0 5
ns.hi_bc             = 0 0 4

# 0 = Interior/Periodic  3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall

amr.blocking_factor     = 8
amr.ref_ratio           = 2 2 2 2

#ns.sum_interval         = 1
#ns.turb_interval        = 10
#ns.liquid_interval      = 1

amr.derive_plot_vars    = energy liquid_water mag_vort humidity

proj.v     = 4
proj.Pcode = 2

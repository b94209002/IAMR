
#*******************************************************************************
# INPUTS.3D.RT for regression testing
#*******************************************************************************

ns.do_mom_diff = 1
ns.do_cons_trac = 1
godunov.use_ppm = 1
godunov.use_force_in_trans = 1 

max_step 		= 7

amr.n_cell              = 128 128 32

amr.max_level           = 0

# Refinement criterion, use vorticity
amr.refinement_indicators = vorticity
amr.vorticity.vorticity_greater =  1000.0

amr.regrid_int		= 2 

ns.v                    = 1
amr.v                   = 1

amr.checkpoint_files_output = 0
amr.check_int		= 10 

amr.plot_int		= 10

ns.cfl                  = 0.7  # CFL number used to set dt

ns.init_shrink          = 1.0  # factor which multiplies the very first time step
ns.init_iter            = 0 

ns.vel_visc_coef        = 0.001
ns.scal_diff_coefs      = 0.0014 0.0014

geometry.coord_sys   =  0

geometry.prob_lo     =  -2. -2. 0.
geometry.prob_hi     =  2. 2. 1.

geometry.is_periodic =  1 1 0

ns.gravity = 0.
ns.do_trac2 = 1

ns.lo_bc             = 0 0 5
ns.hi_bc             = 0 0 4

# 0 = Interior/Periodic  3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall

# Problem parameters
prob.probtype = 12
prob.D0 = 0.
prob.dD = 1.
prob.M0 = 0.
prob.dM = -3.
prob.N2 = 4.
prob.omega = 0.
  
amr.blocking_factor     = 8
amr.ref_ratio           = 2 2 2 2



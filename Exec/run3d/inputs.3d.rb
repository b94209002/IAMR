
#*******************************************************************************
# INPUTS.3D.RT for regression testing
#*******************************************************************************

ns.do_mom_diff = 1
ns.do_cons_trac = 1
ns.do_cons_trac2 = 1

#ns.advection_scheme = BDS
#ns.advection_scheme = godunov
#godunov.use_ppm = 1
godunov.use_force_in_trans = 1

#ns.getForceVerbose = 2

max_step 		= 4000

amr.n_cell              = 384 384 32

amr.max_level           = 0

# Refinement criterion, use vorticity
amr.refinement_indicators = liquid_water
amr.liquid_water.liquid_water_greater =  0.0001

amr.regrid_int		= 2

ns.v                    = 1
amr.v                   = 1

#amr.checkpoint_files_output = 0
#amr.restart             = chk02200
amr.check_int		= 200

amr.plot_int		= 10

ns.cfl                  = 0.7  # CFL number used to set dt

ns.init_shrink          = 1.0  # factor which multiplies the very first time step
ns.init_iter            = 0

ns.vel_visc_coef        = 0.001
ns.scal_diff_coefs      = 0.0014 0.0014

geometry.coord_sys   =  0

geometry.prob_lo     =  0. 0. 0.
geometry.prob_hi     =  12. 12. 1.

geometry.is_periodic =  1 1 0

ns.gravity = 0.
ns.do_trac2 = 1

ns.lo_bc             = 0 0 5
ns.hi_bc             = 0 0 4

# 0 = Interior/Periodic  3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall

amr.plot_vars = x_velocity y_velocity z_velocity density tracer tracer2
amr.derive_plot_vars    = mag_vort  liquid_water

# Problem parameters
prob.probtype = 12
prob.D0 = 0.
prob.dDz = 1.
prob.M0 = 0.
prob.dMz = -3.
prob.N2 = 4.
prob.U0 = 0.0
prob.qrad = 0.
prob.omega = 0.
prob.perturbation_amplitude = 20

amr.blocking_factor     = 8
amr.ref_ratio           = 2 2 2 2

mac_proj.verbose     = 0
mac_proj.semicoarsening = 1
mac_proj.max_semicoarsening_level = 2

diffuse.v = 0
diffuse.semicoarsening = 1
diffuse.max_semicoarsening_level = 2

nodal_proj.verbose = 0
nodal_proj.max_coarsening_level = 5
nodal_proj.semicoarsening = 0
nodal_proj.max_semicoarsening_level = 1
nodal_proj.proj_tol = 1.e-12
#nodal_proj.verbose = 5
#nodal_proj.bottom_verbose = 2
#nodal_proj.bottom_solver = hypre
#nodal_proj.mg_max_coarsening_level=0
#nodal_proj.bottom_rtol       = 1.e-12
#nodal_proj.bottom_maxiter = 100
#hypre.hypre_solver = BiCGSTAB
#hypre.hypre_preconditioner = BoomerAMG


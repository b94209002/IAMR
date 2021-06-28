
#*******************************************************************************
# INPUTS.3D.RT for regression testing
#*******************************************************************************

ns.do_mom_diff = 1
ns.do_cons_trac = 1
godunov.use_ppm = 1
godunov.use_force_in_trans = 1 

max_step 		= 10

amr.n_cell              = 256 256 32

amr.max_level           = 0

# Refinement criterion, use vorticity
amr.refinement_indicators = vorticity
amr.vorticity.vorticity_greater =  1000.0

amr.regrid_int		= 2 

ns.v                    = 1
amr.v                   = 1

#amr.checkpoint_files_output = 0
amr.check_int		= 80 

amr.plot_int		= 10

ns.cfl                  = 0.7  # CFL number used to set dt

ns.init_shrink          = 1.0  # factor which multiplies the very first time step
ns.init_iter            = 0 

ns.vel_visc_coef        = 0.001
ns.scal_diff_coefs      = 0.0014 0.0014

geometry.coord_sys   =  0

geometry.prob_lo     =  -4. -4. 0.
geometry.prob_hi     =  4. 4. 1.

geometry.is_periodic =  1 1 0

ns.gravity = 0.
ns.do_trac2 = 1

ns.lo_bc             = 0 0 5
ns.hi_bc             = 0 0 4

# 0 = Interior/Periodic  3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall

amr.plot_vars = x_velocity y_velocity z_velocity tracer tracer2
amr.derive_plot_vars    = mag_vort  liquid_water

# Problem parameters
prob.probtype = 12
prob.D0 = 0.
prob.dD = 8.
prob.M0 = 0.
prob.dM = -16.
prob.N2 = 24.
prob.qrad = 0.
prob.omega = 0.
prob.perturbation_amplitude = 200

amr.blocking_factor     = 8
amr.ref_ratio           = 2 2 2 2

mac.v     = 0
mac.semicoarsening = 1
mac.max_semicoarsening_level = 2

diffuse.v = 0
diffuse.semicoarsening = 1
diffuse.max_semicoarsening_level = 2

#proj.v = 1
proj.max_coarsening_level = 5
proj.semicoarsening = 1 
proj.max_semicoarsening_level = 1
proj.proj_tol = 1.e-12
#nodal_proj.verbose = 5
#nodal_proj.bottom_verbose = 2
nodal_proj.bottom_solver = hypre
#nodal_proj.mg_max_coarsening_level=0
#nodal_proj.bottom_rtol       = 1.e-12
#nodal_proj.bottom_maxiter = 100
hypre.hypre_solver = BiCGSTAB
hypre.hypre_preconditioner = BoomerAMG


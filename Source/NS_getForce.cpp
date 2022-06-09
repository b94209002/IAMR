
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>

#include <NavierStokes.H>

using namespace amrex;

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// NOTE: This function returns a rho weighted source term.
//
// For conservative (i.e. do_mom_diff=1, do_cons_trac=1), velocities
// are integrated according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q = Scal) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// For non-conservative, this rho-weighted source term will get divided
// by rho in the predict_velocity, velocity_advection, scalar_advection,
// and advection_update routines.
//
// For temperature (which is always non-conservative), we evolve
//
//     dT/dt - U dot grad T = [del dot lambda grad T + S_T] / (rho*c_p)
//     ===> tforces =  S_T/c_p
//
//
// For user-defined forcing, this means
//   - For conservative variables, the force term computed here gets used
//     as-is
//   - For non-conservative variables, the force term computed here is
//     divided by rho before use
//

void
NavierStokesBase::getForce (FArrayBox&       force,
                            const Box&       bx,
                            int              scomp,
                            int              ncomp,
                            const Real       time,
                            const FArrayBox& State,
                            const FArrayBox& Aux,
                            int              auxScomp,
                            const MFIter&    /*mfi*/)
{

   const Real* StateDataPtr  = State.dataPtr();
   const Real* auxDataPtr = Aux.dataPtr(auxScomp);

   const Real  grav     = gravity;

   if (ParallelDescriptor::IOProcessor() && getForceVerbose)
   {
       const int*  f_lo     = force.loVect();
       const int*  f_hi     = force.hiVect();
       const int*  v_lo     = State.loVect();
       const int*  v_hi     = State.hiVect();
       const int*  s_lo     = Aux.loVect();
       const int*  s_hi     = Aux.hiVect();
       
       amrex::Print() << "NavierStokesBase::getForce(): Entered..." << std::endl 
		      << "time      = " << time << std::endl
		      << "scomp     = " << scomp << std::endl
		      << "ncomp     = " << ncomp << std::endl
		      << "auxScomp = " << auxScomp << std::endl;
       
       if  (ncomp==1) amrex::Print() << "Doing only component " << scomp << std::endl;
       else if (scomp==0 && ncomp==AMREX_SPACEDIM) amrex::Print() << "Doing velocities only" << std::endl;
       else if (scomp>=AMREX_SPACEDIM) amrex::Print() << "Doing " << ncomp << " component(s) starting with component " << scomp << std::endl;

       amrex::Print() << "NavierStokesBase::getForce(): Filling Force on box:"
		      << bx << std::endl;
#if (AMREX_SPACEDIM == 3)
       amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
       amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << "," << f_lo[2] << ") - "
		      << "(" << f_hi[0] << "," << f_hi[1] << "," << f_hi[2] << ")" << std::endl;
       amrex::Print() << "NavierStokesBase::getForce(): Vel Domain:" << std::endl;
       amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << "," << v_lo[2] << ") - "
		      << "(" << v_hi[0] << "," << v_hi[1] << "," << v_hi[2] << ")" << std::endl;
       amrex::Print() << "NavierStokesBase::getForce(): Scal Domain:" << std::endl;
       amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << "," << s_lo[2] << ") - "
		      << "(" << s_hi[0] << "," << s_hi[1] << "," << s_hi[2] << ")" << std::endl;
#else
       amrex::Print() << "NavierStokesBase::getForce(): Force Domain:" << std::endl;
       amrex::Print() << "(" << f_lo[0] << "," << f_lo[1] << ") - "
		      << "(" << f_hi[0] << "," << f_hi[1] << ")" << std::endl;
       amrex::Print() << "NavierStokesBase::getForce(): State Domain:" << std::endl;
       amrex::Print() << "(" << v_lo[0] << "," << v_lo[1] << ") - "
		      << "(" << v_hi[0] << "," << v_hi[1] << ")" << std::endl;
       amrex::Print() << "NavierStokesBase::getForce(): Aux Domain:" << std::endl;
       amrex::Print() << "(" << s_lo[0] << "," << s_lo[1] << ") - "
		      << "(" << s_hi[0] << "," << s_hi[1] << ")" << std::endl;
#endif

       // Compute min/max
       Real min, max;
       for (int n=0; n<ncomp; n++) {
	   min= 1.e234;
	   max=-1.e234;
	   amrex::Print() << "State comp " << scomp+n << " min/max " 
			  << State.min<RunOn::Gpu>(scomp+n) << " / "
			  << State.max<RunOn::Gpu>(scomp+n) << std::endl;
       }
       for (int n=auxScomp; n<Aux.nComp(); n++) {
	   min= 1.e234;
	   max=-1.e234;
	   amrex::Print() << "aux comp " << n << " min/max " 
			  << Aux.min<RunOn::Gpu>(n) << " / "
			  << Aux.max<RunOn::Gpu>(n) << std::endl;
       }
   } //end if(getForceVerbose)

   //
   // Here's the meat
   //
   // Velocity forcing
   //
   if ( scomp<AMREX_SPACEDIM ){
       AMREX_ALWAYS_ASSERT(scomp==Xvel);
       AMREX_ALWAYS_ASSERT(ncomp>=AMREX_SPACEDIM);
   }

   if ( scomp==Xvel ){
     //
     // TODO: add some switch for user-supplied/problem-dependent forcing
     //
     auto const& frc = force.array(scomp);
     auto const& aux = Aux.array(auxScomp);
     auto const& state = State.array(auxScomp);

     if ( std::abs(grav) > 0.0001) {
       amrex::ParallelFor(bx, [frc, aux, grav]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
	 frc(i,j,k,0) = Real(0.0);
#if ( AMREX_SPACEDIM == 2 )
         frc(i,j,k,1) = grav*aux(i,j,k,0);
#elif ( AMREX_SPACEDIM == 3 )
         frc(i,j,k,1) = Real(0.0);
         frc(i,j,k,2) = grav*aux(i,j,k,0);
#endif
       });
     }
     else {
       const Real* dom_lo = geom.ProbLo();
       const Real* dx = geom.CellSize();
       NavierStokes::RayleighBenard rb = NavierStokes::getRayleighBenard();

       amrex::ParallelFor(bx, [frc, state, aux, rb, dom_lo, dx]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
#if ( AMREX_SPACEDIM == 2 )
         frc(i,j,k,0) = 0.0;
         Real y = dom_lo[1] + (j + 0.5_rt) * dx[1];
         Real m = aux(i,j,k,2) + rb.M0 + rb.dM*y;
         Real d = aux(i,j,k,1) + rb.D0 + rb.dD*y;
         frc(i,j,k,1) = std::max(m, d - rb.N2*y);
#elif ( AMREX_SPACEDIM == 3 )
	 frc(i,j,k,0) = aux(i,j,k,0) * rb.omega * state(i,j,k,1);
         frc(i,j,k,1) = -aux(i,j,k,0) * rb.omega * state(i,j,k,0);
         Real z = dom_lo[2] + (k + 0.5_rt) * dx[2];
         Real m = aux(i,j,k,2) + rb.M0 + rb.dM*z;
         Real d = aux(i,j,k,1) + rb.D0 + rb.dD*z;
         frc(i,j,k,2) = std::max(m, d - rb.N2*z);
#endif
	 // define dD = (DH-D0)/H and dM = (MH-M0)/H
	 // with this from, DBC = 0 in the buoyancy equation
       });
       // force.setVal<RunOn::Gpu>(0.0, bx, Xvel, AMREX_SPACEDIM);
     }
   }

   //
   // Scalar forcing
   //
   if ( scomp >= AMREX_SPACEDIM || scomp+ncomp >= AMREX_SPACEDIM) {
     // Doing only scalars
     // force.setVal<RunOn::Gpu>(0.0, bx, 0, ncomp);
     // auto const& frc  = force.array();
     // amrex::ParallelFor(bx, ncomp, [frc]
     // AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
     // {
     // 	 frc(i,j,k,n) = 0.0_rt;
     //});
     NavierStokes::RayleighBenard rb = NavierStokes::getRayleighBenard();
     const Real* dom_lo = geom.ProbLo();
     const Real* dom_hi = geom.ProbHi();
     const Real* dx = geom.CellSize();
     const Real Pi = 3.141592653589793238462643383279502884197;

#if ( AMREX_SPACEDIM == 2)

     const Real H = dom_hi[1] - dom_lo[1];

     // We are filling these all at once
     if ( scomp == 0 && scomp+ncomp >= AMREX_SPACEDIM+3 ) 
     {
         auto const& frc = force.array();
         auto const& vel = State.array();
         amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
         AMREX_GPU_DEVICE(int i, int j, int k) noexcept
         {
             Real y = dom_lo[1] + (k + 0.5_rt) * dx[1];
             frc(i,j,k,2) = 0.0_rt;
             frc(i,j,k,3) = -vel(i,j,k,1)*rb.dD - rb.qrad * sin(Pi*y/H);
             frc(i,j,k,4) = -vel(i,j,k,1)*rb.dM - 0.5 * rb.qrad * sin(Pi*y/H);
         });
     }

     // We are filling density, trac and trac2
     if ( scomp == 2 && ncomp >= 3) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real y = dom_lo[1] + (k + 0.5_rt) * dx[1];
         frc(i,j,k,0) = 0.0_rt;
         frc(i,j,k,1) = -vel(i,j,k,1)*rb.dD - rb.qrad * sin(Pi*y/H);
         frc(i,j,k,2) = -vel(i,j,k,1)*rb.dM - 0.5 * rb.qrad * sin(Pi*y/H);
     });
     }

     // We are filling only density 
     if ( scomp == AMREX_SPACEDIM && ncomp == 1 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         frc(i,j,k,0) = 0.0_rt;
     });
     }

     // We are filling only trac 
     if ( scomp == AMREX_SPACEDIM+1 && ncomp == 1 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real y = dom_lo[1] + (k + 0.5_rt) * dx[1];
         frc(i,j,k,0) = -vel(i,j,k,1)*rb.dD - rb.qrad * sin(Pi*y/H);
     });
     }
     // We are filling trac and trac2 
     if ( scomp == AMREX_SPACEDIM+1 && ncomp == 2 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real y = dom_lo[1] + (k + 0.5_rt) * dx[1];
         frc(i,j,k,0) = -vel(i,j,k,1)*rb.dD - rb.qrad * sin(Pi*y/H);
	 frc(i,j,k,1) = -vel(i,j,k,1)*rb.dM - 0.5 * rb.qrad * sin(Pi*y/H);
     });
     }

     // We are filling only trac2
     if ( scomp == AMREX_SPACEDIM+2 && ncomp == 1 ) {
     auto const& frc = force.array();
     auto const& vel = State.array(); 
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real y = dom_lo[1] + (k + 0.5_rt) * dx[1];
         frc(i,j,k,0) = -vel(i,j,k,1)*rb.dM - 0.5 * rb.qrad * sin(Pi*y/H);
     });
     }
#elif ( AMREX_SPACEDIM == 3)
     const Real H = dom_hi[2] - dom_lo[2];

     // We are filling these all at once
     if ( scomp == 0 && scomp+ncomp >= AMREX_SPACEDIM+3 )
     {
         auto const& frc = force.array();
         auto const& vel = State.array();
         amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
         AMREX_GPU_DEVICE(int i, int j, int k) noexcept
         {
             Real z = dom_lo[2] + (k + 0.5_rt) * dx[2];
             frc(i,j,k,3) = 0.0_rt;
             frc(i,j,k,4) = -vel(i,j,k,2)*rb.dD - rb.qrad * sin(Pi*z/H);
             frc(i,j,k,5) = -vel(i,j,k,2)*rb.dM - 0.5 * rb.qrad * sin(Pi*z/H);
         });
     }
     // We are filling scalers at once
     if ( scomp == AMREX_SPACEDIM && ncomp == 3 ) {
     auto const& frc = force.array();
     auto const& vel = State.array(); 
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real z = dom_lo[2] + (k + 0.5_rt) * dx[2];
         frc(i,j,k,0) = 0.0_rt;
         frc(i,j,k,1) = -vel(i,j,k,2)*rb.dD - rb.qrad * sin(Pi*z/H);
         frc(i,j,k,2) = -vel(i,j,k,2)*rb.dM - 0.5 *rb.qrad * sin(Pi*z/H);
     });
     }

     // We are filling only density 
     if ( scomp == AMREX_SPACEDIM && ncomp == 1 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         frc(i,j,k,0) = 0.0_rt;
     });
     }
     // We are filling only trac 
     if ( scomp == AMREX_SPACEDIM+1 && ncomp == 1 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real z = dom_lo[2] + (k + 0.5_rt) * dx[2];
         frc(i,j,k,0) = -vel(i,j,k,2)*rb.dD - rb.qrad * sin(Pi*z/H);
     });
     }
     // We are filling trac and trac2
     if ( scomp == AMREX_SPACEDIM+1 && ncomp == 2 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real z = dom_lo[2] + (k + 0.5_rt) * dx[2];
         frc(i,j,k,0) = -vel(i,j,k,2)*rb.dD - rb.qrad * sin(Pi*z/H);
	 frc(i,j,k,1) = -vel(i,j,k,2)*rb.dM - 0.5 * rb.qrad * sin(Pi*z/H);
     });
     }
     // We are filling trac2
     if ( scomp == AMREX_SPACEDIM+2 && ncomp == 1 ) {
     auto const& frc = force.array();
     auto const& vel = State.array();
     amrex::ParallelFor(bx, [frc, vel, rb, dx, H, dom_lo, Pi]
     AMREX_GPU_DEVICE(int i, int j, int k) noexcept
     {
         Real z = dom_lo[2] + (k + 0.5_rt) * dx[2];
         frc(i,j,k,0) = -vel(i,j,k,2)*rb.dM - 0.5 * rb.qrad * sin(Pi*z/H);
     });
     }
#endif     

   }

   if (ParallelDescriptor::IOProcessor() && getForceVerbose) {
       // Compute min/max
       Real min, max;
       for (int n=0; n<ncomp; n++) {
	   min= 1.e234;
	   max=-1.e234;
	   amrex::Print() << "Force comp " << scomp+n << " min/max " 
			  << force.min<RunOn::Gpu>(scomp+n) << " / "
			  << force.max<RunOn::Gpu>(scomp+n) << std::endl;
       }

      amrex::Print() << "NavierStokesBase::getForce(): Leaving..." 
                     << std::endl << "---" << std::endl;
   }
}

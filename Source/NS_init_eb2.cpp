#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>

using namespace amrex;

// called in main before Amr->init(start,stop)
void
initialize_EB2 (const Geometry& geom, const int required_coarsening_level,
		const int max_coarsening_level)
{
    // read in EB parameters
    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

    Print()<<"geom "<<geom<<"\n\n";
    //Build geometry -- WIP!!!
    // seems like this is a call that can create the implicit function,
    // GeometryShop, and Geometry all at "once"... geom might already be init
    //AMReX_EB2::Build (const Geometry& geom, int required_coarsening_level,
    //   int max_coarsening_level, int ngrow)
    //fixme -- not sure that 4 is the right number here, just taken from CNS
    EB2::Build(geom, required_coarsening_level, max_coarsening_level, 4);
}
#endif
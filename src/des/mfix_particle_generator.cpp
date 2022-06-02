#include <mfix_particle_generator.H>
#include <mfix_calc_cell.H>
#include <mfix_ic_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_des_parts_gen_K.H>

#include <AMReX_AmrParGDB.H>

#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
//                                                                      !
//  Purpose: Generate particle configuration based on maximum particle  !
//           radius and filling from top to bottom within specified     !
//           bounds                                                     !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

using namespace amrex;


ParticlesGenerator::ParticlesGenerator (const amrex::IntVect& bx_lo,
                                        const amrex::IntVect& bx_hi,
                                        const amrex::GpuArray<amrex::Real,3>& plo,
                                        const amrex::GpuArray<amrex::Real,3>& dx,
                                        const int id,
                                        const int cpu,
                                        const int icv,
                                        const int type)
  : m_bx_lo(bx_lo)
  , m_bx_hi(bx_hi)
  , m_plo(plo)
  , m_dx(dx)
  , m_id(id)
  , m_cpu(cpu)
  , m_icv(icv)
  , m_type(type)
{}


void
ParticlesGenerator::generate (int& particles_count,
                              ParticleTileType& particles)
{
  std::string ic_pack_type_str = IC::ic[m_icv].packing;

  int cube_base(-1);

  const size_t cube_pos = ic_pack_type_str.find("-cube");

  // Get the X-ncube integer input
  if(cube_pos != std::string::npos) {
    const size_t cube_len = ic_pack_type_str.length();
    std::string str_base = ic_pack_type_str.erase(cube_pos, cube_len-cube_pos);
    cube_base = std::stoi(str_base);
  }

  // Handle oneper and eightper cases
  if (ic_pack_type_str.compare("oneper") == 0)
    cube_base = 1;
  else if (ic_pack_type_str.compare("eightper") == 0)
    cube_base = 2;

  // Call generate with specific positions generator
  if (DEM::solve && ic_pack_type_str.compare("hcp") == 0) {

    Hex_ClosePack hex_close_pack(m_plo, m_dx);
    hex_close_pack.setup(m_bx_lo, m_bx_hi, m_icv, m_type);
    generate(particles_count, particles, hex_close_pack);

  } else if (DEM::solve && ic_pack_type_str.compare("random") == 0) {

    m_h_data.clear();
    m_d_data.clear();

    RandomFill_DEM random_fill_dem(m_plo, m_dx);
    random_fill_dem.setup(m_bx_lo, m_bx_hi, m_icv, m_type, m_h_data, m_d_data, false);
    generate(particles_count, particles, random_fill_dem);

  } else if (PIC::solve && ic_pack_type_str.compare("random") == 0) {

    RandomFill_PIC random_fill_pic(m_plo, m_dx);
    random_fill_pic.setup(m_bx_lo, m_bx_hi, m_icv, m_type, false);
    generate(particles_count, particles, random_fill_pic);

  } else if (DEM::solve && ic_pack_type_str.compare("pseudo_random") == 0) {

    m_h_data.clear();
    m_d_data.clear();

    RandomFill_DEM random_fill_dem(m_plo, m_dx);
    random_fill_dem.setup(m_bx_lo, m_bx_hi, m_icv, m_type, m_h_data, m_d_data, true);

    generate(particles_count, particles, random_fill_dem);

  } else if (PIC::solve && ic_pack_type_str.compare("pseudo_random") == 0) {

    RandomFill_PIC random_fill_pic(m_plo, m_dx);
    random_fill_pic.setup(m_bx_lo, m_bx_hi, m_icv, m_type, true);
    generate(particles_count, particles, random_fill_pic);

  } else if (cube_base > 0) {

    nCubePer_Fill n_cube_per_fill(cube_base, m_plo, m_dx);
    n_cube_per_fill.setup(m_bx_lo, m_bx_hi, m_icv, m_type);
    generate(particles_count, particles, n_cube_per_fill);

  } else {

    amrex::Abort("Unknown particle generator fill type");
  }
}


template <typename F1>
void ParticlesGenerator::generate (int& particles_count,
                                   ParticleTileType& particles,
                                   F1 positions_generator)
{
  particles_count = positions_generator.get_particles_number();

  const int current_size = particles.numParticles();

  particles.resize(current_size + particles_count);

  const SOLIDS_t& solid = IC::ic[m_icv].solids[m_type];

  // Setup particle diameters parameters

  int diameter_distr_uniform(solid.diameter.is_uniform());
  int diameter_distr_normal(solid.diameter.is_normal());
  Real diameter_mean = solid.diameter.get_mean();
  Real diameter_stddev = solid.diameter.get_stddev();
  Real diameter_min = solid.diameter.get_min();
  Real diameter_max = solid.diameter.get_max();

  // Setup particle densities parameters

  int density_distr_uniform(solid.density.is_uniform());
  int density_distr_normal(solid.density.is_normal());
  Real density_mean(solid.density.get_mean());
  Real density_stddev(solid.density.get_stddev());
  Real density_min(solid.density.get_min());
  Real density_max(solid.density.get_max());

  // Get particles initial velocity
  const Real ic_u_s = solid.velocity[0];
  const Real ic_v_s = solid.velocity[1];
  const Real ic_w_s = solid.velocity[2];

  const Real statwt = solid.statwt;

  auto& aos = particles.GetArrayOfStructs();
  ParticleType* pstruct = aos().dataPtr();

  auto& soa = particles.GetStructOfArrays();
  auto p_realarray = soa.realarray();
  auto p_intarray = soa.intarray();

  const Real picmulti = (DEM::solve) ? 1.0 : 0.0;

  const int type = m_type;
  const int local_cg_dem=DEM::cg_dem;
  const int id = m_id;
  const int cpu = m_cpu;

  amrex::ParallelForRNG(particles_count, [pstruct,p_realarray,p_intarray,
      picmulti,ic_u_s,ic_v_s,ic_w_s,statwt,type,id,cpu,local_cg_dem,
      diameter_distr_uniform,diameter_distr_normal,diameter_mean,current_size,
      diameter_stddev,diameter_min,diameter_max,density_distr_uniform,
      density_distr_normal,density_mean,density_stddev,density_min,density_max,
      positions_generator]
    AMREX_GPU_DEVICE (int p, RandomEngine const& engine) noexcept
  {
    const int p_tot = current_size + p;

    ParticleType& part = pstruct[p_tot];

    part.id() = id+p;
    part.cpu() = cpu;

    RealVect position = positions_generator.template operator()<run_on>(p, engine);
    part.pos(0) = position[0];
    part.pos(1) = position[1];
    part.pos(2) = position[2];

    Real diameter = 0;

    if (diameter_distr_uniform) {
      diameter = diameter_min + (diameter_max-diameter_min)*amrex::Random(engine);
    } else if (diameter_distr_normal) {
      diameter = amrex::RandomNormal(diameter_mean, diameter_stddev, engine);
      while (diameter < diameter_min || diameter_max < diameter){
        diameter = amrex::RandomNormal(diameter_mean, diameter_stddev, engine);
      }
    } else {
      diameter = diameter_mean;
    }

    if (local_cg_dem) {
      diameter *= std::cbrt(statwt);
    }

    Real rad = .5*diameter;

    Real rho = 0;

    if (density_distr_uniform) {
      rho = density_min + (density_max-density_min)*amrex::Random(engine);
    } else if (density_distr_normal) {
      rho = amrex::RandomNormal(density_mean, density_stddev, engine);
      while(rho < density_min || density_max < rho) {
        rho = amrex::RandomNormal(density_mean, density_stddev, engine);
      }
    } else {
      rho = density_mean;
    }

    Real vol  = (4.0/3.0)*M_PI*rad*rad*rad;
    Real mass = vol * rho;
    Real omoi = 2.5/(mass * rad*rad);

    p_realarray[SoArealData::velx][p_tot] = ic_u_s;
    p_realarray[SoArealData::vely][p_tot] = ic_v_s;
    p_realarray[SoArealData::velz][p_tot] = ic_w_s;

    p_realarray[SoArealData::statwt][p_tot] = statwt;

    p_realarray[SoArealData::radius][p_tot] = rad;
    p_realarray[SoArealData::density][p_tot] = rho;

    p_realarray[SoArealData::volume][p_tot] = vol;
    p_realarray[SoArealData::mass][p_tot] = mass;
    p_realarray[SoArealData::oneOverI][p_tot] = omoi*picmulti;

    p_realarray[SoArealData::omegax][p_tot] = 0.0;
    p_realarray[SoArealData::omegay][p_tot] = 0.0;
    p_realarray[SoArealData::omegaz][p_tot] = 0.0;

    p_realarray[SoArealData::dragcoeff][p_tot] = 0.0;

    p_realarray[SoArealData::dragx][p_tot] = 0.0;
    p_realarray[SoArealData::dragy][p_tot] = 0.0;
    p_realarray[SoArealData::dragz][p_tot] = 0.0;

    p_intarray[SoAintData::phase][p_tot] = type+1;
    p_intarray[SoAintData::state][p_tot] = 1;
  });

  return;
}


//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
//                                                                     !
//                                                                     !
//                                                                     !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void
ParticlesGenerator::write (const int nrp,
                           ParticleTileType& particles,
                           const int nstep) const
{
  const auto& aos = particles.GetArrayOfStructs();
  const ParticleType* pstruct = aos().dataPtr();

  auto& soa = particles.GetStructOfArrays();
  auto p_realarray = soa.realarray();

  std::ostringstream nstep_stream;
  nstep_stream << nstep;

  std::string nstep_string(nstep_stream.str());
  std::string filename = "test-" + nstep_string + ".vtp";

  std::ofstream output_file;
  output_file.open(filename.c_str());

  // Write the necessary header information for a PolyData file type
  output_file << "<?xml version=\"1.0\"?>" << std::endl;
  output_file << "<VTKFile type=\"PolyData\"" << " "
              << "version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  output_file << "   " << "<PolyData>" << std::endl;

  // Write Piece tag and identify the number of parts in the system.
  output_file << "      " << "<Piece NumberOfPoints=\""
              << std::fixed << std::setw(10) << std::setprecision(10)
              << std::setfill('0') << nrp << "\" NumberOfVerts=\"0\" "
              << "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl;

  output_file << "         " << "<PointData>" << std::endl;

  output_file << "            "
    << "<DataArray type=\"Float32\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">"
    << std::endl;

  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << std::setfill(' ') << Real(p_realarray[SoArealData::radius][lc1]) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "            "
    << "<DataArray type=\"Float32\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">"
    << std::endl;

  for(int lc1 = 0; lc1 < nrp; ++lc1)
  {
    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << Real(p_realarray[SoArealData::density][lc1]) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;

  output_file << "         " << "</PointData>" << std::endl;

  output_file << "         " << "<Points>" << std::endl;

  output_file << "            " << "<DataArray type=\"Float32\" "
              << "Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">"
              << std::endl;

  for(int lc1 = 0; lc1 < nrp; lc1++)
  {
    const ParticleType& part = pstruct[lc1];

    output_file << "               "
                << std::scientific << std::setw(13) << std::setprecision(6)
                << Real(part.pos(0)) << "    "
                << Real(part.pos(1)) << "    "
                << Real(part.pos(2)) << std::endl;
  }

  output_file << "            " << "</DataArray>" << std::endl;
  output_file << "         " << "</Points>" << std::endl;

  // Write tags for data not included (vtp format style)
  output_file << "         " << "<Verts></Verts>" << std::endl;
  output_file << "         " << "<Lines></Lines>" << std::endl;
  output_file << "         " << "<Strips></Strips>" << std::endl;
  output_file << "         " << "<Polys></Polys>" << std::endl;

  output_file << "      " << "</Piece>" << std::endl;
  output_file << "   " << "</PolyData>" << std::endl;
  output_file << "</VTKFile>" << std::endl;

  output_file.close();

  return;
}

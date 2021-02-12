target_include_directories(mfixcore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( mfixcore
   PRIVATE
   mfix_mol_predict.cpp
   mfix_mol_predict_vels_on_faces_K.H
   mfix_mol_predict_vels_on_faces_eb_K.H
   mfix_mol_convective_fluxes.cpp
   mfix_mol_convective_fluxes_eb.cpp
   mfix_mol_divergence.cpp
   mfix_mol_flux_redistribute.cpp
   )

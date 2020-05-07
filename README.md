# Hyperspectral_Simulation_Using_HTFRTC
Using HTFRTC to Simulate Hyper-Spectral Satellite Radiance

# How to Compile
Tested on UW S4 Cluster

mpiifort -O3 -traceback Calculate_Radiance.f90 -I/data/users/qzhang/local/intel_18.0.3/netcdf-4.1.2/include -I/data/users/qzhang/local/intel_18.0.3/rttov-12.3/mod -I/data/users/qzhang/local/intel_18.0.3/rttov-12.3/include -L/data/users/qzhang/local/intel_18.0.3/netcdf-4.1.2/lib -lnetcdff -lnetcdf -L/data/users/qzhang/local/intel_18.0.3/rttov-12.3/lib -lrttov_emis_atlas -lrttov_brdf_atlas -lrttov_other -lrttov_mw_scatt -lrttov_other -lrttov_coef_io -lrttov_hdf -lrttov_main -L/data/users/qzhang/local/intel_18.0.3/hdf5-1.8.19/lib -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz -L/data/users/qzhang/local/intel_18.0.3/lapack-3.5.0 -llapack -lrefblas -o HTFRTC_Cal.exe

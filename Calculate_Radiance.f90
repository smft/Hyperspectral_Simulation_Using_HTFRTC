! alculate HTFRTC Needed Variables From WRF Output
! Pressure in hPa
! Potential Temperature in K
! Water Vapor in kg/kg
! 2m AGL Temperature in K
! 2m Water Vapor in K
! Surface Pressure in hPa
! 10m Zonal Wind in m/s
! 10m Meridinal Wind in m/s
! Skin Temperature in K
! Terrain Height in m
! Landsea Cover
! Cloud Fraction
! Satellite Zenth angle

program test
  use mpi
  use netcdf
  use rttov_const, only: errorstatus_success,errorstatus_fatal
  use rttov_types, only: rttov_options,rttov_coefs,rttov_profile,&
                          rttov_transmission,rttov_radiance,rttov_chanprof,&
                          rttov_emissivity,rttov_pccomp
  use mod_rttov_emis_atlas, only: rttov_emis_atlas_data,atlas_type_ir
  use parkind1, only: jpim,jplm,jprb
  use rttov_unix_env, only: rttov_exit
  implicit none
  include "rttov_direct.interface"
  include "rttov_parallel_direct.interface"
  include "rttov_read_coefs_htfrtc.interface"
  include "rttov_dealloc_coefs.interface"
  include "rttov_alloc_prof.interface"
  include "rttov_alloc_pccomp.interface"
  include "rttov_print_opts.interface"
  include "rttov_print_profile.interface"
  include "rttov_skipcommentline.interface"
  include "rttov_setup_emis_atlas.interface"
  include "rttov_get_emis.interface"
  include "rttov_deallocate_emis_atlas.interface"

  character(len=999) :: nc_name,fname_coef,fname_sensor,emis_path,save_name
  character(len=999),parameter :: var_name_1="T2",var_name_2="TSK",var_name_3="Q2",&
                                  var_name_4="PSFC",var_name_5="U10",var_name_6="V10",&
                                  var_name_7="HGT",var_name_8="P",var_name_9="PB",&
                                  var_name_10="QVAPOR",var_name_11="T",var_name_12="CLDFRA",&
                                  var_name_13="XLAT",var_name_14="XLONG",var_name_15="XLAND"
  integer :: nlat,nlon,nlevels,month,npcscores,i,j,k,total_count=1,rank,size,ierr
  integer :: cell_chunk,cell_mod,cell_profiles,nchannels,y_dimid,x_dimid,dimids(2)
  integer :: ncid,var_name_1_id,var_name_2_id,var_name_3_id,var_name_4_id,var_name_5_id,&
             var_name_6_id,var_name_7_id,var_name_8_id,var_name_9_id,var_name_10_id,&
             var_name_11_id,var_name_12_id,var_name_13_id,var_name_14_id,var_name_15_id,&
             var_tbb_id,var_lat_id,var_lon_id
  real :: base_temperature,sat_lat,sat_lon,start_sec,end_sec
  real,dimension(:,:,:),allocatable :: pressure,p,pb,watervapor,theta,cloudfra,&
                                       temperature
  real,dimension(:,:),allocatable :: t2,tsk,q2,psfc,u10,v10,&
                                     terrain,lat,lon,surface_type,&
                                     pressure_2d,temperature_2d,watervapor_2d
  real,dimension(:,:),allocatable :: pressure_recv,temperature_recv,watervapor_recv,&
                                     tbb,tbb_all
  real,dimension(:),allocatable :: t2_1d,tsk_1d,q2_1d,psfc_1d,u10_1d,v10_1d,&
                                     cloudtop_1d,terrain_1d,lat_1d,lon_1d,z_angle_1d,&
                                     cloudfra_1d,surface_type_1d
  real,dimension(:),allocatable :: surface_type_recv,t2_recv,tsk_recv,q2_recv,&
                                      psfc_recv,u10_recv,v10_recv,cloudtop_recv,&
                                      terrain_recv,lat_recv,lon_recv,z_angle_recv,&
                                      cloudfra_recv
  namelist /HTFRTC/ nc_name,fname_coef,fname_sensor,emis_path,save_name,nlat,nlon,nlevels,&
                    base_temperature,sat_lat,sat_lon,month,npcscores,nchannels
  open(15,file='htfrtc.namelist')
  read(15,HTFRTC)
  close(15)
  allocate(t2(nlon,nlat),tsk(nlon,nlat),q2(nlon,nlat),&
           psfc(nlon,nlat),u10(nlon,nlat),v10(nlon,nlat),&
           terrain(nlon,nlat),lat(nlon,nlat),lon(nlon,nlat),&
           surface_type(nlon,nlat))
  allocate(pressure(nlon,nlat,nlevels),p(nlon,nlat,nlevels),pb(nlon,nlat,nlevels),&
           watervapor(nlon,nlat,nlevels),theta(nlon,nlat,nlevels),cloudfra(nlon,nlat,nlevels),&
           temperature(nlon,nlat,nlevels))
  allocate(pressure_2d(nlon*nlat,nlevels),temperature_2d(nlon*nlat,nlevels),&
           watervapor_2d(nlon*nlat,nlevels))
  allocate(t2_1d(nlon*nlat),tsk_1d(nlon*nlat),q2_1d(nlon*nlat),psfc_1d(nlon*nlat),&
           u10_1d(nlon*nlat),v10_1d(nlon*nlat),cloudtop_1d(nlon*nlat),terrain_1d(nlon*nlat),&
           lat_1d(nlon*nlat),lon_1d(nlon*nlat),z_angle_1d(nlon*nlat),surface_type_1d(nlon*nlat),&
           cloudfra_1d(nlon*nlat))
  call MPI_INIT(ierr) 
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  cell_chunk=(nlon*nlat)/(size-1)
  cell_mod=mod(nlon*nlat,(size-1))
if (rank .eq. 0) then
  call check(nf90_open(trim(nc_name),nf90_nowrite,ncid))
  call check(nf90_inq_varid(ncid,trim(var_name_1),var_name_1_id))
  call check(nf90_inq_varid(ncid,trim(var_name_2),var_name_2_id))
  call check(nf90_inq_varid(ncid,trim(var_name_3),var_name_3_id))
  call check(nf90_inq_varid(ncid,trim(var_name_4),var_name_4_id))
  call check(nf90_inq_varid(ncid,trim(var_name_5),var_name_5_id))
  call check(nf90_inq_varid(ncid,trim(var_name_6),var_name_6_id))
  call check(nf90_inq_varid(ncid,trim(var_name_7),var_name_7_id))
  call check(nf90_inq_varid(ncid,trim(var_name_8),var_name_8_id))
  call check(nf90_inq_varid(ncid,trim(var_name_9),var_name_9_id))
  call check(nf90_inq_varid(ncid,trim(var_name_10),var_name_10_id))
  call check(nf90_inq_varid(ncid,trim(var_name_11),var_name_11_id))
  call check(nf90_inq_varid(ncid,trim(var_name_12),var_name_12_id))
  call check(nf90_inq_varid(ncid,trim(var_name_13),var_name_13_id))
  call check(nf90_inq_varid(ncid,trim(var_name_14),var_name_14_id))
  call check(nf90_inq_varid(ncid,trim(var_name_15),var_name_15_id))
  call check(nf90_get_var(ncid,var_name_1_id,t2,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_2_id,tsk,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_3_id,q2,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_4_id,psfc,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_5_id,u10,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_6_id,v10,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_7_id,terrain,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_13_id,lat,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_14_id,lon,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_15_id,surface_type,(/1,1,1/),(/nlon,nlat,1/)))
  call check(nf90_get_var(ncid,var_name_8_id,p,(/1,1,1,1/),(/nlon,nlat,nlevels,1/)))
  call check(nf90_get_var(ncid,var_name_9_id,pb,(/1,1,1,1/),(/nlon,nlat,nlevels,1/)))
  call check(nf90_get_var(ncid,var_name_10_id,watervapor,(/1,1,1,1/),(/nlon,nlat,nlevels,1/)))
  call check(nf90_get_var(ncid,var_name_11_id,theta,(/1,1,1,1/),(/nlon,nlat,nlevels,1/)))
  call check(nf90_get_var(ncid,var_name_12_id,cloudfra,(/1,1,1,1/),(/nlon,nlat,nlevels,1/)))
  call check(nf90_close(ncid))
  pressure=(p+pb)/100.0
  temperature=(theta+base_temperature)*((pressure/1000.0)**.286)
  surface_type=surface_type-1
  do i=1,nlon
    do j=1,nlat
      t2_1d(total_count)=t2(i,j)
      tsk_1d(total_count)=tsk(i,j)
      q2_1d(total_count)=q2(i,j)
      psfc_1d(total_count)=psfc(i,j)
      u10_1d(total_count)=u10(i,j)
      v10_1d(total_count)=v10(i,j)
      cloudtop_1d(total_count)=cloud_top_pressure(cloudfra(i,j,:),pressure(i,j,:),nlevels)
      terrain_1d(total_count)=terrain(i,j)
      lat_1d(total_count)=lat(i,j)
      lon_1d(total_count)=lon(i,j)
      z_angle_1d(total_count)=zenith_angle(lat(i,j),lon(i,j),sat_lat,sat_lon)
      surface_type_1d(total_count)=surface_type(i,j)
      cloudfra_1d(total_count)=min(sum(cloudfra(i,j,:)),1.0)
      pressure_2d(total_count,:)=pressure(i,j,:)
      temperature_2d(total_count,:)=temperature(i,j,:)
      watervapor_2d(total_count,:)=watervapor(i,j,:)
      total_count=total_count+1
    end do
  end do
  do i=1,(size-1)
    if (i .lt. (size-1)) then
      call MPI_SEND(pressure_2d(1+cell_chunk*(i-1):i*cell_chunk,:nlevels),cell_chunk*nlevels,&
                    MPI_REAL,i,1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(temperature_2d(1+cell_chunk*(i-1):i*cell_chunk,:nlevels),cell_chunk*nlevels,&
                    MPI_REAL,i,2,MPI_COMM_WORLD,ierr)
      call MPI_SEND(watervapor_2d(1+cell_chunk*(i-1):i*cell_chunk,:nlevels),cell_chunk*nlevels,&
                    MPI_REAL,i,3,MPI_COMM_WORLD,ierr)
      call MPI_SEND(t2_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,4,MPI_COMM_WORLD,ierr)
      call MPI_SEND(q2_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,5,MPI_COMM_WORLD,ierr)
      call MPI_SEND(tsk_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,6,MPI_COMM_WORLD,ierr)
      call MPI_SEND(psfc_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,7,MPI_COMM_WORLD,ierr)
      call MPI_SEND(cloudtop_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,8,MPI_COMM_WORLD,ierr)
      call MPI_SEND(cloudfra_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,9,MPI_COMM_WORLD,ierr)
      call MPI_SEND(terrain_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,10,MPI_COMM_WORLD,ierr)
      call MPI_SEND(lat_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,11,MPI_COMM_WORLD,ierr)
      call MPI_SEND(lon_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,12,MPI_COMM_WORLD,ierr)
      call MPI_SEND(surface_type_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,13,MPI_COMM_WORLD,ierr)
      call MPI_SEND(u10_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,14,MPI_COMM_WORLD,ierr)
      call MPI_SEND(v10_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,15,MPI_COMM_WORLD,ierr)
      call MPI_SEND(z_angle_1d(1+cell_chunk*(i-1):i*cell_chunk),cell_chunk,&
                    MPI_REAL,i,16,MPI_COMM_WORLD,ierr)
    else
      call MPI_SEND(pressure_2d(1+cell_chunk*(i-1):,:nlevels),(nlon*nlat-cell_chunk*(i-1))*nlevels,&
                    MPI_REAL,i,1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(temperature_2d(1+cell_chunk*(i-1):,:nlevels),(nlon*nlat-cell_chunk*(i-1))*nlevels,&
                    MPI_REAL,i,2,MPI_COMM_WORLD,ierr)
      call MPI_SEND(watervapor_2d(1+cell_chunk*(i-1):,:nlevels),(nlon*nlat-cell_chunk*(i-1))*nlevels,&
                    MPI_REAL,i,3,MPI_COMM_WORLD,ierr)
      call MPI_SEND(t2_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,4,MPI_COMM_WORLD,ierr)
      call MPI_SEND(q2_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,5,MPI_COMM_WORLD,ierr)
      call MPI_SEND(tsk_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,6,MPI_COMM_WORLD,ierr)
      call MPI_SEND(psfc_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,7,MPI_COMM_WORLD,ierr)
      call MPI_SEND(cloudtop_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,8,MPI_COMM_WORLD,ierr)
      call MPI_SEND(cloudfra_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,9,MPI_COMM_WORLD,ierr)
      call MPI_SEND(terrain_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,10,MPI_COMM_WORLD,ierr)
      call MPI_SEND(lat_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,11,MPI_COMM_WORLD,ierr)
      call MPI_SEND(lon_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,12,MPI_COMM_WORLD,ierr)
      call MPI_SEND(surface_type_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,13,MPI_COMM_WORLD,ierr)
      call MPI_SEND(u10_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,14,MPI_COMM_WORLD,ierr)
      call MPI_SEND(v10_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,15,MPI_COMM_WORLD,ierr)
      call MPI_SEND(z_angle_1d(1+cell_chunk*(i-1):),nlon*nlat-cell_chunk*(i-1),&
                    MPI_REAL,i,16,MPI_COMM_WORLD,ierr)
    end if
  end do
  allocate(tbb_all(nlon*nlat,nchannels+2))
else
  if (rank .lt. (size-1)) then
    allocate(pressure_recv(cell_chunk,nlevels))
    allocate(temperature_recv(cell_chunk,nlevels))
    allocate(watervapor_recv(cell_chunk,nlevels))
    allocate(surface_type_recv(cell_chunk))
    allocate(t2_recv(cell_chunk))
    allocate(tsk_recv(cell_chunk))
    allocate(q2_recv(cell_chunk))
    allocate(psfc_recv(cell_chunk))
    allocate(u10_recv(cell_chunk))
    allocate(v10_recv(cell_chunk))
    allocate(cloudtop_recv(cell_chunk))
    allocate(terrain_recv(cell_chunk))
    allocate(lat_recv(cell_chunk))
    allocate(lon_recv(cell_chunk))
    allocate(z_angle_recv(cell_chunk))
    allocate(cloudfra_recv(cell_chunk))
    call MPI_RECV(pressure_recv,cell_chunk*nlevels,MPI_REAL,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(temperature_recv,cell_chunk*nlevels,MPI_REAL,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(watervapor_recv,cell_chunk*nlevels,MPI_REAL,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(t2_recv,cell_chunk*1,MPI_REAL,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(q2_recv,cell_chunk,MPI_REAL,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(tsk_recv,cell_chunk,MPI_REAL,0,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(psfc_recv,cell_chunk,MPI_REAL,0,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(cloudtop_recv,cell_chunk,MPI_REAL,0,8,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(cloudfra_recv,cell_chunk,MPI_REAL,0,9,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(terrain_recv,cell_chunk,MPI_REAL,0,10,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(lat_recv,cell_chunk,MPI_REAL,0,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(lon_recv,cell_chunk,MPI_REAL,0,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(surface_type_recv,cell_chunk,MPI_REAL,0,13,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(u10_recv,cell_chunk,MPI_REAL,0,14,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(v10_recv,cell_chunk,MPI_REAL,0,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(z_angle_recv,cell_chunk,MPI_REAL,0,16,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    cell_profiles=cell_chunk
  else
    allocate(pressure_recv(nlon*nlat-cell_chunk*(rank-1),nlevels))
    allocate(temperature_recv(nlon*nlat-cell_chunk*(rank-1),nlevels))
    allocate(watervapor_recv(nlon*nlat-cell_chunk*(rank-1),nlevels))
    allocate(surface_type_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(t2_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(tsk_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(q2_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(psfc_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(u10_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(v10_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(cloudtop_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(terrain_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(lat_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(lon_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(z_angle_recv(nlon*nlat-cell_chunk*(rank-1)))
    allocate(cloudfra_recv(nlon*nlat-cell_chunk*(rank-1)))
    call MPI_RECV(pressure_recv,(cell_chunk+cell_mod)*nlevels,MPI_REAL,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(temperature_recv,(cell_chunk+cell_mod)*nlevels,MPI_REAL,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(watervapor_recv,(cell_chunk+cell_mod)*nlevels,MPI_REAL,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(t2_recv,(nlon*nlat-cell_chunk*(rank-1))*1,MPI_REAL,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(q2_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(tsk_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(psfc_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(cloudtop_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,8,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(cloudfra_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,9,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(terrain_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,10,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(lat_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(lon_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(surface_type_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,13,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(u10_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,14,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(v10_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(z_angle_recv,nlon*nlat-cell_chunk*(rank-1),MPI_REAL,0,16,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    cell_profiles=cell_chunk+cell_mod
    allocate(tbb(cell_profiles,nchannels+2))
  end if
  call cpu_time(start_sec)
  tbb=calculate_htfrtc(temperature_recv,&
                       pressure_recv,&
                       watervapor_recv,&
                       t2_recv,&
                       q2_recv,&
                       psfc_recv,&
                       u10_1d,&
                       v10_1d,&
                       tsk_recv,&
                       terrain_recv,&
                       lat_recv,&
                       lon_recv,&
                       z_angle_recv,&
                       surface_type_recv,&
                       cloudtop_recv,&
                       cloudfra_recv,&
                       nlevels,&
                       npcscores,&
                       month,&
                       cell_profiles,&
                       fname_coef,&
                       fname_sensor,&
                       emis_path)
  call cpu_time(end_sec)
  print '("HTFRTC Calculation on Rank ",i3,", Time Use = ",f8.3," seconds.")',rank,end_sec-start_sec
  call MPI_SEND(tbb,cell_profiles*nchannels,MPI_REAL,0,rank,MPI_COMM_WORLD,ierr)
end if

if (rank .eq. 0) then
  do i=1,(size-1)
    if (i .lt. (size-1)) then
      call MPI_RECV(tbb_all(1+cell_chunk*(i-1):i*cell_chunk,:),&
                    cell_chunk*nchannels,MPI_REAL,&
                    i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    else
      call MPI_RECV(tbb_all(1+cell_chunk*(i-1):,:),&
                    (nlon*nlat-cell_chunk*(i-1))*nchannels,MPI_REAL,&
                    i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end if
  end do
  !do i=i,nchannels
  !  write(*,*) sum(tbb_all(:,i))/(nlon*nlat)
  !end do
  call cpu_time(start_sec)
  call check(nf90_create(save_name,NF90_CLOBBER,ncid)) 
  call check(nf90_def_dim(ncid,"nlocs",nlon*nlat,y_dimid))
  call check(nf90_def_dim(ncid,"nchannels",nchannels,x_dimid))
  dimids=(/y_dimid,x_dimid/)
  call check(nf90_def_var(ncid,"IASI_TBB",NF90_DOUBLE,dimids,var_tbb_id))
  call check(nf90_def_var(ncid,"LAT",NF90_DOUBLE,(/y_dimid/),var_lat_id))
  call check(nf90_def_var(ncid,"LON",NF90_DOUBLE,(/y_dimid/),var_lon_id))
  call check(nf90_enddef(ncid))
  call check(nf90_put_var(ncid,var_tbb_id,tbb_all(:,:nchannels)))
  call check(nf90_put_var(ncid,var_lat_id,tbb_all(:,nchannels+1)))
  call check(nf90_put_var(ncid,var_lon_id,tbb_all(:,nchannels+2)))
  call check(nf90_close(ncid))
  call cpu_time(end_sec)
  print '("Write NetCDF File, Time Use = ",f8.3,"seconds.")',end_sec-start_sec
end if

contains
  subroutine check(status)
    integer,intent(in) :: status
    if(status/=nf90_noerr) then 
      write(*,*), trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  function zenith_angle(lat,lon,sat_lat,sat_lon) result(za)
    real,parameter :: DTR=0.017453293   ! degrees to radians
    real,parameter :: RE=6370.0         ! radius of earth
    real,parameter :: PI=3.1415926      ! pi
    real,parameter :: RTD=57.29577951   ! radians to degrees
    real,parameter :: RES=40589641.0    ! (earth radius)**2
    real,parameter :: C1=1826835337.0   ! RES + RZS
    real,parameter :: C2=538527888.0    ! 2*RE*RES
    real,parameter :: RZS=1786245696.0  ! (geosynchronous height)**2 for GOES
    real,intent(in) :: sat_lat,sat_lon
    real,intent(in) :: lat,lon
    real :: za
    real :: A,B,C,COSD,COSE,E,P,PP
    A=(90.0-lat)*DTR  ! colatitude of point
    B=(90.0-sat_lat)*DTR   ! colatitude of satellite
    if (lon>180.) then
      C=abs(lon-360.-sat_lon)*DTR
    else
      C=abs(lon-sat_lon)*DTR
    end if    
    COSD=cos(A)*cos(B)+sin(A)*sin(B)*sin(C) ! great circle arc
    PP=C1-C2*COSD
    P=sqrt(PP)
    COSE=(PP+RES-RZS)/(2.*RE*P)
    COSE=max(min(COSE,1.),-1.)
    E=acos(COSE)
    za=PI-E
    za=za*RTD
    za=max(za,0.)
  end function zenith_angle

  function cloud_top_pressure(cloudfra,pressure,nlevels) result(cloudtop)
    integer :: k
    integer,intent(in) :: nlevels
    real,intent(in) :: cloudfra(nlevels),pressure(nlevels)
    real :: cloudtop,trans(nlevels)
    trans=trans+999999999.9
    do k=1,nlevels
      if (cloudfra(k)>0) then
        trans(k)=pressure(k)
      end if
    end do
    cloudtop=minval(trans)
  end function cloud_top_pressure

  function calculate_htfrtc(temperature,pressure,watervapor,&
                  t2,q2,psfc,u10,v10,tsk,height,lat,lon,&
                  z_angle,surface_type,cloudtop,cloudfra,nlevels,&
                  npcscores,month,profs,&
                  fname_coef,fname_sensor,emis_path) result(transmistance)

    integer(kind=jpim),intent(in) :: nlevels
    integer(kind=jpim),intent(in) :: npcscores,month ! number of pcscores (1-300)
    integer(kind=jpim),intent(in) :: profs
    character(len=*),intent(in) :: fname_coef,fname_sensor,emis_path
    real,intent(in) :: surface_type(profs)
    real,intent(in) :: temperature(profs,nlevels),pressure(profs,nlevels),watervapor(profs,nlevels)
    real,intent(in) :: t2(profs),q2(profs),psfc(profs),u10(profs),v10(profs),tsk(profs),height(profs),&
                       lat(profs),lon(profs),z_angle(profs),cloudtop(profs),cloudfra(profs)
    real,dimension(:,:),allocatable :: transmistance

    integer(kind=jpim) :: iup ! unit for input profile file
    integer(kind=jpim) :: ioout,i
    integer(kind=jpim) :: nchannels_rec
    type(rttov_options) :: opts ! Options structure
    type(rttov_coefs) :: coefs ! Coefficients structure
    type(rttov_chanprof), pointer :: chanprof(:) => NULL() ! Required argument,unused by HTFRTC
    logical(kind=jplm), pointer :: calcemis(:) => NULL() ! Optional, flag to indicate calculation of emissivity within HTFRTC
    type(rttov_emissivity), pointer :: emissivity(:)  => NULL() ! Optional,input/output surface emissivity
    type(rttov_profile), pointer :: profiles(:) => NULL() ! Input profiles
    type(rttov_transmission) :: transmission ! Required argument, unused by HTFRTC
    type(rttov_radiance) :: radiance ! Required argument, unused by HTFRTC
    type(rttov_pccomp) :: pccomp ! Output PC structure
    type(rttov_emis_atlas_data) :: emis_atlas ! Data structure for emissivity atlas
    integer(kind=jpim) :: errorstatus ! Return error status of RTTOV subroutine calls
    integer(kind=jpim) :: nprof
    integer(kind=jpim) :: nchanprof
    integer(KIND=jpim) :: lo,hi
    integer :: ios
    errorstatus=0_jpim
    nprof=profs

    opts%htfrtc_opts%htfrtc=.true. ! Select HTFRTC
    opts%htfrtc_opts%n_pc_in=npcscores
    opts%htfrtc_opts%simple_cloud=.true. ! Set to true to enable simple cloud scheme
    opts%htfrtc_opts%overcast=.true. ! Set to true to enable overcast radiance calculations
    opts%rt_ir%ozone_data=.false.
    opts%rt_ir%co2_data=.false.
    opts%rt_ir%n2o_data=.false.
    opts%rt_ir%ch4_data=.false.
    opts%rt_ir%co_data=.false.
    opts%rt_ir%so2_data=.false.
    opts%htfrtc_opts%reconstruct=.true.

    call rttov_read_coefs_htfrtc(errorstatus,coefs,fname_coef,fname_sensor)
    if (errorstatus/=errorstatus_success) then
      write(*,*) 'fatal error reading coefficients'
      call rttov_exit(errorstatus)
    end if
    nchannels_rec=coefs%coef_htfrtc%n_ch
    nchanprof=nchannels_rec
    allocate(transmistance(profs,nchannels_rec+2))
    allocate(profiles(1),chanprof(nchanprof),calcemis(nchanprof),&
              emissivity(nchanprof),stat=errorstatus)
    if(errorstatus/=errorstatus_success) then
      write(*,*) 'allocation error'
      call rttov_exit(errorstatus)
    end if
    call rttov_alloc_prof(errorstatus,1,profiles,nlevels,opts,asw=1_jpim,init=.TRUE._jplm)
    if(errorstatus/=errorstatus_success) then
      write(*,*) 'allocation error for profiles structure'
      call rttov_exit(errorstatus)
    end if
    call rttov_alloc_pccomp(errorstatus,pccomp=pccomp,npcscores=npcscores,asw=1_jpim,&
                          init=.TRUE._jplm,nchannels_rec=nchannels_rec,&
                          opts=opts,nlevels=nlevels)
    if(errorstatus/=errorstatus_success) then
      write(*,*) 'allocation error for pccomp structure'
      call rttov_exit(errorstatus)
    end if

    if(month>=1 .and. month<=12) then
      call rttov_setup_emis_atlas(errorstatus,opts,month,atlas_type_ir,emis_atlas,path=emis_path,coefs=coefs)
      if(errorstatus/=errorstatus_success) then
        write(*,*) 'error initialising emissivity atlas'
        call rttov_exit(errorstatus)
      end if
    end if
    if(month >= 1 .and. month<= 12) then
      call rttov_get_emis(errorstatus,opts,chanprof,profiles,coefs,emis_atlas,emissivity(:)%emis_in)
      if(errorstatus/=errorstatus_success) then
        write(*,*) 'error reading emissivity atlas'
        call rttov_exit(errorstatus)
      end if
      calcemis(:)=(emissivity(:)%emis_in<=0._jprb)
    else
      emissivity(:)%emis_in=0._jprb
      calcemis(:)=.true.
    end if
    do i=1,nprof
      profiles(1)%gas_units=1
      profiles(1)%p(:)=pressure(i,nlevels:1:-1)
      profiles(1)%t(:)=temperature(i,nlevels:1:-1)
      profiles(1)%q(:)=watervapor(i,nlevels:1:-1)
      profiles(1)%s2m%t=t2(i)
      profiles(1)%s2m%q=q2(i)
      profiles(1)%s2m%p=psfc(i)
      profiles(1)%s2m%u=u10(i)
      profiles(1)%s2m%v=v10(i)
      profiles(1)%skin%t=tsk(i)
      profiles(1)%skin%surftype=int(surface_type(i),kind=jpim)
      profiles(1)%elevation=height(i)
      profiles(1)%latitude=lat(i)
      profiles(1)%longitude=lon(i)
      profiles(1)%zenangle=z_angle(i)
      if(cloudtop(i) .eq. 999999999.9) then
        profiles(1)%ctp=cloudtop(i)
        profiles(1)%cfraction=cloudfra(i)
      else
        profiles(1)%ctp=0
        profiles(1)%cfraction=0
      end if
      call rttov_direct(errorstatus,chanprof,opts,profiles,coefs,transmission,radiance,&
                        calcemis=calcemis,emissivity=emissivity,pccomp=pccomp)
      lo=1
      hi=nchannels_rec
      transmistance(i,lo:hi)=pccomp%bt_pccomp(lo:hi)
      transmistance(i,hi+1)=lat(i)
      transmistance(i,hi+2)=lon(i)
    end do
    call rttov_dealloc_coefs(errorstatus,coefs)
    deallocate(profiles,chanprof,calcemis,emissivity)
    call rttov_deallocate_emis_atlas(emis_atlas)
  end function calculate_htfrtc

end program test

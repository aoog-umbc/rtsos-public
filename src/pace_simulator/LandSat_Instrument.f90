MODULE LandSat_INSTRUMENT_DESIGN
USE HDF5
INTEGER,PARAMETER :: NWV_BAND=26,NWV_BAND_SEG1=17, NWV_BAND_SEG2=4,&
                     NWV_BAND_SEG1P2=21,NWV_BAND_SEG3=5
! NWV_BAND_SEG1: Virtual Bands from 430 nm - 676 nm
! NWV_BAND_SEG2: NIR + SWIR
! NWV_BAND_SEG3: CoastalAerosol+Blue+Green+Red+Pan

!REAL*8,PARAMETER :: WV_START=340D0 !,WV_END=1000.D0

REAL*8,DIMENSION(NWV_BAND):: WV_BAND,LandSat_FWHM=0.0d0

INTEGER,PARAMETER :: NILS_LandSat_VIRT=319
INTEGER,DIMENSION(9) :: NILS_LandSat

REAL*8,DIMENSION(9,NILS_LandSat_VIRT) :: ILS_WaveLength
REAL*8,DIMENSION(NWV_BAND,NILS_LandSat_VIRT) :: ILS_DeltaWaveLength,ILS_LandSat

CONTAINS

SUBROUTINE WV_LandSat_SETUP(aux_dir)
IMPLICIT NONE
INTEGER IWV,ILS_LINE
REAL*8 Delta_Wavelength_Max,Delta_Wavelength_Sub,Factor

! HDF 5 DEFINITION
! This should map to REAL*8 on most modern processors
!INTEGER, PARAMETER :: real_kind_15 = SELECTED_REAL_KIND(Fortran_REAL_8)
CHARACTER(LEN=160),intent(in) :: aux_dir
CHARACTER(LEN=180) :: LandSat_RSR_Filename

INTEGER(HID_T)  :: file, space, dset, attr ! Handles
INTEGER :: hdferr
INTEGER(hsize_t),   DIMENSION(1:2) :: dims
INTEGER(hsize_t),DIMENSION(1) :: dimscl
INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
TYPE(C_PTR) :: f_ptr
REAL*4,DIMENSION(1),TARGET :: HDF5RTMP
INTEGER,DIMENSION(1),TARGET :: HDF5ITMP
!REAL*4,DIMENSION(:),TARGET,ALLOCATABLE :: HDF5RARRD1
REAL*8,DIMENSION(:,:),TARGET,ALLOCATABLE :: HDF5RARR2DIM

IF(NWV_BAND /= NWV_BAND_SEG1+NWV_BAND_SEG2+NWV_BAND_SEG3) STOP 'CHECK NWV_BAND'
IF(NWV_BAND_SEG1P2 /= NWV_BAND_SEG1+NWV_BAND_SEG2) STOP 'CHECK NWV_BAND_SEG1P2'


LandSat_FWHM(1:NWV_BAND_SEG1)=16.0D0
LandSat_FWHM(NWV_BAND_SEG1+1)=28.25D0
LandSat_FWHM(NWV_BAND_SEG1+2)=20.39D0
LandSat_FWHM(NWV_BAND_SEG1+3)=84.72D0
LandSat_FWHM(NWV_BAND_SEG1+4)=186.66D0

LandSat_FWHM(NWV_BAND_SEG1P2+1)=15.98D0
LandSat_FWHM(NWV_BAND_SEG1P2+2)=60.04D0
LandSat_FWHM(NWV_BAND_SEG1P2+3)=57.33D0
LandSat_FWHM(NWV_BAND_SEG1P2+4)=37.47D0
LandSat_FWHM(NWV_BAND_SEG1P2+5)=172.40D0

WV_BAND(1)=434.0D0
DO IWV=2,NWV_BAND_SEG1
   WV_BAND(IWV)=WV_BAND(IWV-1)+LandSat_FWHM(IWV)
ENDDO
WV_BAND(NWV_BAND_SEG1+1)=864.67D0
WV_BAND(NWV_BAND_SEG1+2)=1373.43D0
WV_BAND(NWV_BAND_SEG1+3)=1608.86D0
WV_BAND(NWV_BAND_SEG1+4)=2200.73D0

WV_BAND(NWV_BAND_SEG1P2+1)=442.96D0
WV_BAND(NWV_BAND_SEG1P2+2)=482.04D0
WV_BAND(NWV_BAND_SEG1P2+3)=561.41D0
WV_BAND(NWV_BAND_SEG1P2+4)=654.59D0
WV_BAND(NWV_BAND_SEG1P2+5)=589.50D0

NILS_LandSat(1)=72 ! FOR 864.67 NM
NILS_LandSat(2)=70 ! FOR 1373.43 NM
NILS_LandSat(3)=183 ! FOR 1608 NM
NILS_LandSat(4)=319 ! FOR 2200 NM

NILS_LandSat(5)=33 ! FOR 442.96 NM
NILS_LandSat(6)=93 ! FOR 482.04 NM
NILS_LandSat(7)=99 ! FOR 561.41 NM
NILS_LandSat(8)=67 ! FOR 654.59 NM
NILS_LandSat(9)=205 ! FOR 589.50 NM

LandSat_RSR_Filename=trim(aux_dir)//'/Ball_BA_RSR.v1.2.h5'

ILS_WaveLength=0.0D0

CALL h5open_f(hdferr)
CALL h5fopen_f(LandSat_RSR_Filename, H5F_ACC_RDONLY_F, file, hdferr)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(1)))
CALL h5dopen_f (file, '/RSR_NIR', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(1,1:NILS_LandSat(1))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+1,1:NILS_LandSat(1))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(2)))
CALL h5dopen_f (file, '/RSR_Cirrus', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(2,1:NILS_LandSat(2))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+2,1:NILS_LandSat(2))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(3)))
CALL h5dopen_f (file, '/RSR_SWIR1', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(3,1:NILS_LandSat(3))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+3,1:NILS_LandSat(3))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(4)))
CALL h5dopen_f (file, '/RSR_SWIR2', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(4,1:NILS_LandSat(4))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+4,1:NILS_LandSat(4))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(5)))
CALL h5dopen_f (file, '/RSR_CA', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(5,1:NILS_LandSat(5))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+5,1:NILS_LandSat(5))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(6)))
CALL h5dopen_f (file, '/RSR_Blue', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(6,1:NILS_LandSat(6))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+6,1:NILS_LandSat(6))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(7)))
CALL h5dopen_f (file, '/RSR_Green', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(7,1:NILS_LandSat(7))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+7,1:NILS_LandSat(7))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(8)))
CALL h5dopen_f (file, '/RSR_Red', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(8,1:NILS_LandSat(8))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+8,1:NILS_LandSat(8))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

ALLOCATE(HDF5RARR2DIM(3,NILS_LandSat(9)))
CALL h5dopen_f (file, '/RSR_Pan', dset, hdferr)
CALL h5dget_space_f(dset,space, hdferr)
f_ptr = C_LOC(HDF5RARR2DIM(1,1))
CALL h5dread_f( dset, H5T_NATIVE_DOUBLE,f_ptr, hdferr)
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
ILS_WaveLength(9,1:NILS_LandSat(9))=HDF5RARR2DIM(1,:)
ILS_LandSat(NWV_BAND_SEG1+9,1:NILS_LandSat(9))=HDF5RARR2DIM(2,:)
DEALLOCATE(HDF5RARR2DIM)

CALL h5fclose_f(file , hdferr)
CALL h5close_f(hdferr)

ILS_DeltaWaveLength=0.0d0
DO IWV=NWV_BAND_SEG1+1,NWV_BAND
DO ILS_LINE=1,NILS_LandSat(IWV-NWV_BAND_SEG1)
   ILS_DeltaWaveLength(IWV,ILS_LINE)= &
         ILS_WaveLength(IWV-NWV_BAND_SEG1,ILS_LINE)-WV_BAND(IWV)
ENDDO
ENDDO

DO IWV=1,NWV_BAND_SEG1

Delta_Wavelength_Max=2.0d0*LandSat_FWHM(IWV)
Delta_Wavelength_Sub=Delta_Wavelength_Max/(1.0d0*NILS_LandSat_VIRT)
Factor=4.0D0*log(2.0D0)/(LandSat_FWHM(IWV)**2)

DO ILS_LINE=1,NILS_LandSat_VIRT
   ILS_DeltaWaveLength(IWV,ILS_LINE)=(ILS_LINE-NILS_LandSat_VIRT/2)*Delta_Wavelength_Sub
   ILS_LandSat(IWV,ILS_LINE)=exp(-Factor*(ILS_DeltaWaveLength(IWV,ILS_LINE)**2));
ENDDO
ENDDO

! RESET RSR IF SMALLER THAN 0
DO IWV=NWV_BAND_SEG1+1,NWV_BAND
DO ILS_LINE=1,NILS_LandSat_VIRT
  IF(ILS_LandSat(IWV,ILS_LINE)<0.0D0)ILS_LandSat(IWV,ILS_LINE)=0.0D0
ENDDO
ENDDO

END SUBROUTINE WV_LandSat_SETUP

END MODULE LandSat_INSTRUMENT_DESIGN

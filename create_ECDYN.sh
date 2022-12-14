# This is a script in Bash.

# Author: Lionel GUEZ

# This script creates the file "ECDYN.nc" at the first day, 0h, of a
# given year and month, from ERA interim data at IDRIS or on
# Ciclad. You need to be in the UNIX group subipsl at IDRIS, or in the
# UNIX group ecmwf on Ciclad, to access the data. "ECDYN.nc" also
# contains ST and CDSW (which could be in "ECPHY.nc").CDSW is not
# modified because we do not care about that variable (and we do not
# know what to put into it).

# This script calls the accompanying Python script "lnsp_to_SP.py",
# which should be in the same directory.

set -xe

year=2005
month=03

IDRIS=y

if [[ $IDRIS == y ]]
then
    ERAI_dir=ergon:~rpsl376/ERAI/NETCDF/GLOBAL_075/4xdaily
    module load python/3.3.2

    export PYTHONPATH=$PYTHONPATH:/workgpfs/rech/lmd/rlmd542/lib/python3.3/site-packages
    # (for the netCDF4 python module at IDRIS)
else
    # Ciclad:
    ERAI_dir=/bdd/ERAI/NETCDF/GLOBAL_075/4xdaily
    module load python/3.4-anaconda3
fi

for var in u v ta r
do
    rcp $ERAI_dir/AN_PL/$year/$var.$year$month.aphei.GLOBAL_075.nc \
	${var}_${year}_$month.nc
done

for var in geopt stl1
do
    rcp $ERAI_dir/AN_SF/$year/$var.$year$month.ashei.GLOBAL_075.nc \
	${var}_${year}_$month.nc
done

rcp $ERAI_dir/AN_ML/$year/lnsp.$year$month.amhei.GLOBAL_075.nc \
    lnsp_${year}_$month.nc

for var in u v ta r geopt stl1 lnsp
do
    chmod 600 ${var}_${year}_$month.nc
    ncpdq -U --dimension=time,0 --overwrite ${var}_${year}_$month.nc \
	${var}_${year}_${month}_01_00.nc
done

# Remove degenerate dimension "level" and variable "level" in lnsp
# file in preparation of paste:
ncwa -a level --overwrite lnsp_${year}_${month}_01_00.nc plouf.nc
ncks --exclude --variable=level --overwrite plouf.nc \
    lnsp_${year}_${month}_01_00.nc

# Start from CDSW so that the time value is overwritten:

if [[ $IDRIS == y ]]
then
    rcp ergon:~rpsl035/IGCM/INIT/ATM/LMDZ/ECDYN.nc.20020101 ECDYN_high_res.nc
else
    # Ciclad:
    wget http://www.lmd.jussieu.fr/~lmdz/LMDZ_Init/ECDYN_high_res.nc
fi

ncks --variable=CDSW --overwrite ECDYN_high_res.nc \
    ECDYN_${year}_${month}_01_00.nc

for var in u v ta r geopt stl1 lnsp
do
    ncks --append ${var}_${year}_${month}_01_00.nc \
	ECDYN_${year}_${month}_01_00.nc
done

ncrename --variable=u,U --variable=v,V --variable=ta,TEMP --variable=r,R \
    --variable=geopt,Z --variable=stl1,ST ECDYN_${year}_${month}_01_00.nc

script_dir=`dirname $0`
$script_dir/lnsp_to_SP.py ECDYN_${year}_${month}_01_00.nc

# Clean up:

for var in u v ta r geopt stl1 lnsp
do
    rm ${var}_${year}_$month.nc ${var}_${year}_${month}_01_00.nc
done

rm -f plouf.nc ECDYN_high_res.nc

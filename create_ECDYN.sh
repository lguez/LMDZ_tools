# This is a script in Bash.

# Author: Lionel GUEZ

# This script creates the file "ECDYN.nc" at the first day, 0h, of a
# given year and month, from ERA interim data at IDRIS or on
# Spirit. On Spirit, you need to be in the UNIX group ecmwf to access
# the data. "ECDYN.nc" also contains ST and CDSW (which could be in
# "ECPHY.nc"). CDSW is not modified because we do not care about that
# variable (and we do not know what to put into it).

# This script calls the accompanying Python script "lnsp_to_SP.py",
# which should be in the same directory.

set -e

year=2005
month=03
IDRIS=y
module purge

if [[ $IDRIS == y ]]
then
    ERAI_dir=/gpfsstore/rech/psl/rpsl376/ergon/ERAI/NETCDF/GLOBAL_075/4xdaily
    module load nco
else
    # Spirit:
    ERAI_dir=/bdd/ERAI/NETCDF/GLOBAL_075/4xdaily
fi

set -x

for var in u v ta r
do
    cp $ERAI_dir/AN_PL/$year/$var.$year$month.aphei.GLOBAL_075.nc \
	${var}_${year}_$month.nc
done

for var in geopt stl1
do
    cp $ERAI_dir/AN_SF/$year/$var.$year$month.ashei.GLOBAL_075.nc \
	${var}_${year}_$month.nc
done

cp $ERAI_dir/AN_ML/$year/lnsp.$year$month.amhei.GLOBAL_075.nc \
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
cp $R_IN/ATM/ECDYN.nc.20020101 ECDYN_high_res.nc
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

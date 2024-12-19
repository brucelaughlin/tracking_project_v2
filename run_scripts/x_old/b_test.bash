cmd = '''testFile='/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed/wc15_fiechter_era5_glorys_avg_1995.nc'
testLine="$(ncdump -h $testFile | head -n 3 | tail -1)"
echo "$testLine"'''

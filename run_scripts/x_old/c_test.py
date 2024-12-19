import re
import subprocess

pattern = re.compile('\\s*ocean_time = UNLIMITED ; \\/\\/ \\(([0-9]+) currently\\)')

testFile = '/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed/wc15_fiechter_era5_glorys_avg_1997.nc'

cmd = f'''testLine="$(ncdump -h {testFile})"
echo "$testLine"'''

# Apparently this is "horrible" according to the SO author I copied
bash_output = subprocess.run(cmd,shell=True,text=True,capture_output=True,check=True)
bash_output = bash_output.stdout

re_data = re.search(pattern, bash_output)
